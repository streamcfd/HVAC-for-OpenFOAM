/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    ISO7243Foam

Description
    Post-processing utility for ISO 7243 heat stress assessment via WBGT.

    Output fields (all in degC, stored as dimensionless scalars):
    - TnwbISO7243_C       : Natural wet-bulb temperature
    - TgISO7243_C         : Globe temperature
    - WBGTIndoorISO7243_C : WBGT for indoor/no-solar conditions
    - WBGTOutdoorISO7243_C: WBGT for outdoor/solar conditions
    - WBGTISO7243_C       : Selected WBGT formula (indoor or outdoor)
    - ISO7243RiskLevel    : 0..3 based on configurable WBGT thresholds

    Input handling:
    - Required field: T [K]
    - RH field is searched automatically (RH, relativeHumidity, thermoRelHum, relHum)
    - If no RH field is available, RH from ISO7243Dict is used
    - If Tnwb field is absent, Stull approximation from Ta/RH is used
    - If Tg field is absent, estimate from Tmrt/G/IDefault/qr/wall-temperature fallback

    Notes:
    - For strict ISO 7243 compliance, provide physically based Tnwb and Tg fields.
    - Built-in Tnwb/Tg estimation is intended for practical CFD post-processing.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "cellSet.H"
#include "wallFvPatch.H"
#include "mathematicalConstants.H"

namespace Foam
{
namespace ISO7243Constants
{
    const scalar stefanBoltzmann = 5.670374419e-8;  // W/(m^2 K^4)
    const scalar kelvinOffset = 273.15;
}
}

Foam::scalar clampRH(const Foam::scalar rh)
{
    return Foam::max(0.0, Foam::min(100.0, rh));
}

Foam::scalar toCelsiusAuto(const Foam::scalar value)
{
    // Auto-detect common OpenFOAM temperature conventions:
    // values above ~170 are treated as Kelvin, otherwise degC.
    return (value > 170.0) ? value - ISO7243Constants::kelvinOffset : value;
}

Foam::scalar estimateNaturalWetBulbStull
(
    const Foam::scalar taC,
    const Foam::scalar rhPercent
)
{
    // Stull (2011) approximation, valid for typical HVAC/ambient ranges.
    const scalar rh = Foam::max(1.0, Foam::min(100.0, rhPercent));

    return
        taC*Foam::atan(0.151977*Foam::sqrt(rh + 8.313659))
      + Foam::atan(taC + rh)
      - Foam::atan(rh - 1.676331)
      + 0.00391838*Foam::pow(rh, 1.5)*Foam::atan(0.023101*rh)
      - 4.686035;
}

Foam::scalar estimateRadiantTemperatureFromG(const Foam::scalar gWm2)
{
    const scalar g = Foam::max(SMALL, gWm2);
    return
        Foam::pow(g/(4.0*ISO7243Constants::stefanBoltzmann), 0.25)
      - ISO7243Constants::kelvinOffset;
}

Foam::scalar estimateRadiantTemperatureFromQr
(
    const Foam::scalar qrWm2,
    const Foam::scalar taK,
    const Foam::scalar emissivity
)
{
    const scalar eps = Foam::max(0.01, Foam::min(1.0, emissivity));
    const scalar denom = Foam::max(SMALL, eps*ISO7243Constants::stefanBoltzmann);
    scalar trad4 = Foam::pow(taK, 4) + qrWm2/denom;
    trad4 = Foam::max(SMALL, trad4);
    return Foam::pow(trad4, 0.25) - ISO7243Constants::kelvinOffset;
}

Foam::scalar calculateWallAverageTemperatureC(const Foam::fvMesh& mesh)
{
    const volScalarField& T = mesh.lookupObject<volScalarField>("T");

    scalar weightedWallTemp(0.0);
    scalar totalWallArea(0.0);

    const fvPatchList& patches = mesh.boundary();
    forAll(patches, patchI)
    {
        if (!isType<wallFvPatch>(patches[patchI]))
        {
            continue;
        }

        const label pI = patches[patchI].index();
        const scalar area = gSum(mesh.magSf().boundaryField()[pI]);

        if (area > SMALL)
        {
            weightedWallTemp += gSum
            (
                mesh.magSf().boundaryField()[pI] * T.boundaryField()[pI]
            );
            totalWallArea += area;
        }
    }

    if (totalWallArea > SMALL)
    {
        return (weightedWallTemp/totalWallArea) - ISO7243Constants::kelvinOffset;
    }

    // Fallback when no wall patches exist in region
    scalar weightedCellTemp(0.0);
    scalar totalVolume(0.0);
    forAll(mesh.C(), cellI)
    {
        weightedCellTemp += mesh.V()[cellI]*T[cellI];
        totalVolume += mesh.V()[cellI];
    }

    reduce(weightedCellTemp, sumOp<scalar>());
    reduce(totalVolume, sumOp<scalar>());

    return (weightedCellTemp/Foam::max(SMALL, totalVolume))
         - ISO7243Constants::kelvinOffset;
}

Foam::labelList getCellsToAnalyze(const Foam::fvMesh& mesh, const Foam::argList& args)
{
    labelList cellsToAnalyze;

    if (args.found("cellSet"))
    {
        const word cellSetName = args.get<word>("cellSet");
        Info<< "Loading cellSet: " << cellSetName << endl;

        cellSet selectedCells(mesh, cellSetName);
        if (selectedCells.empty())
        {
            FatalErrorInFunction
                << "cellSet " << cellSetName << " is empty or does not exist"
                << abort(FatalError);
        }

        cellsToAnalyze = selectedCells.toc();
    }
    else if (args.found("cellZone"))
    {
        const word cellZoneName = args.get<word>("cellZone");
        Info<< "Loading cellZone: " << cellZoneName << endl;

        const cellZoneMesh& cellZones = mesh.cellZones();
        const label zoneId = cellZones.findZoneID(cellZoneName);

        if (zoneId == -1)
        {
            FatalErrorInFunction
                << "cellZone " << cellZoneName << " does not exist" << nl
                << "Available cellZones: " << cellZones.names()
                << abort(FatalError);
        }

        cellsToAnalyze = cellZones[zoneId];
    }
    else
    {
        cellsToAnalyze.setSize(mesh.nCells());
        forAll(cellsToAnalyze, i)
        {
            cellsToAnalyze[i] = i;
        }
    }

    const label nCells = returnReduce(cellsToAnalyze.size(), sumOp<label>());
    Info<< "Analyzing " << nCells << " cells" << endl;

    return cellsToAnalyze;
}

Foam::scalar computeRiskLevel
(
    const Foam::scalar wbgtC,
    const Foam::scalar lowLimit,
    const Foam::scalar moderateLimit,
    const Foam::scalar highLimit
)
{
    if (wbgtC > highLimit)
    {
        return 3.0;
    }

    if (wbgtC > moderateLimit)
    {
        return 2.0;
    }

    if (wbgtC > lowLimit)
    {
        return 1.0;
    }

    return 0.0;
}

int main(int argc, char *argv[])
{
    Foam::timeSelector::addOptions();

    argList::addBoolOption
    (
        "outdoorSolar",
        "Use ISO 7243 outdoor (with solar) formula for WBGT output"
    );

    argList::addBoolOption
    (
        "indoorNoSolar",
        "Use ISO 7243 indoor/no-solar formula for WBGT output"
    );

    argList::addOption
    (
        "cellSet",
        "name",
        "Analyze only cells in the specified cellSet"
    );

    argList::addOption
    (
        "cellZone",
        "name",
        "Analyze only cells in the specified cellZone"
    );

    #include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"

    if (args.found("outdoorSolar") && args.found("indoorNoSolar"))
    {
        FatalErrorInFunction
            << "Use either -outdoorSolar or -indoorNoSolar, not both."
            << abort(FatalError);
    }

    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);
        Info<< "\nTime = " << runTime.timeName() << endl;

        IOobject dictHeader
        (
            "ISO7243Dict",
            runTime.time().constant(),
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        );

        autoPtr<IOdictionary> dictPtr;
        if (dictHeader.typeHeaderOk<IOdictionary>())
        {
            dictPtr.reset(new IOdictionary(dictHeader));
            Info<< "Using settings from constant/ISO7243Dict" << endl;
        }
        else
        {
            Info<< "constant/ISO7243Dict not found: using internal defaults" << endl;
        }

        const scalar defaultRH = dictPtr.valid()
            ? dictPtr().lookupOrDefault<scalar>("RH", 50.0)
            : 50.0;

        const bool useOutdoorFromDict = dictPtr.valid()
            ? dictPtr().lookupOrDefault<bool>("useOutdoorSolarFormula", false)
            : false;

        const scalar globeOffsetC = dictPtr.valid()
            ? dictPtr().lookupOrDefault<scalar>("globeTemperatureOffset", 0.0)
            : 0.0;

        const scalar globeEmissivity = dictPtr.valid()
            ? dictPtr().lookupOrDefault<scalar>("globeEmissivity", 0.95)
            : 0.95;

        const scalar lowRiskLimit = dictPtr.valid()
            ? dictPtr().lookupOrDefault<scalar>("lowRiskLimit", 25.0)
            : 25.0;
        const scalar moderateRiskLimit = dictPtr.valid()
            ? dictPtr().lookupOrDefault<scalar>("moderateRiskLimit", 28.0)
            : 28.0;
        const scalar highRiskLimit = dictPtr.valid()
            ? dictPtr().lookupOrDefault<scalar>("highRiskLimit", 31.0)
            : 31.0;

        if (!(lowRiskLimit < moderateRiskLimit && moderateRiskLimit < highRiskLimit))
        {
            FatalErrorInFunction
                << "WBGT thresholds must satisfy: lowRiskLimit < moderateRiskLimit < highRiskLimit"
                << abort(FatalError);
        }

        bool useOutdoorFormula = useOutdoorFromDict;
        if (args.found("outdoorSolar"))
        {
            useOutdoorFormula = true;
        }
        if (args.found("indoorNoSolar"))
        {
            useOutdoorFormula = false;
        }

        const word tNwbFieldName = dictPtr.valid()
            ? dictPtr().lookupOrDefault<word>("tNwbFieldName", "Tnwb")
            : "Tnwb";
        const word tgFieldName = dictPtr.valid()
            ? dictPtr().lookupOrDefault<word>("tgFieldName", "Tg")
            : "Tg";
        const word tmrtFieldName = dictPtr.valid()
            ? dictPtr().lookupOrDefault<word>("tmrtFieldName", "Tmrt")
            : "Tmrt";

        IOobject THeader
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );
        volScalarField T(THeader, mesh);

        // Optional humidity source fields (searched in priority order)
        autoPtr<volScalarField> rhField;
        word rhSource("ISO7243Dict RH");
        {
            // Prioritize humidity fields written by buoyantHumiditySimpleFoam.
            const wordList rhCandidates
            ({
                "thermo:relHum",
                "thermoRelHum",
                "relHum",
                "RH",
                "relativeHumidity"
            });

            forAll(rhCandidates, i)
            {
                IOobject rhHeader
                (
                    rhCandidates[i],
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                );

                if (rhHeader.typeHeaderOk<volScalarField>())
                {
                    rhField.reset(new volScalarField(rhHeader, mesh));
                    rhSource = rhCandidates[i];
                    break;
                }
            }
        }

        autoPtr<volScalarField> tNwbField;
        word tNwbSource("Stull(T,RH) approximation");
        {
            IOobject tNwbHeader
            (
                tNwbFieldName,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            );

            if (tNwbHeader.typeHeaderOk<volScalarField>())
            {
                tNwbField.reset(new volScalarField(tNwbHeader, mesh));
                tNwbSource = tNwbFieldName;
            }
        }

        autoPtr<volScalarField> tgField;
        autoPtr<volScalarField> tmrtField;
        autoPtr<volScalarField> gField;
        autoPtr<volScalarField> qrField;
        autoPtr<volScalarField> iDefaultField;

        word tgSource("wallTemperature fallback");
        {
            IOobject tgHeader
            (
                tgFieldName,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            );
            if (tgHeader.typeHeaderOk<volScalarField>())
            {
                tgField.reset(new volScalarField(tgHeader, mesh));
                tgSource = tgFieldName;
            }

            IOobject tmrtHeader
            (
                tmrtFieldName,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            );
            if (tmrtHeader.typeHeaderOk<volScalarField>())
            {
                tmrtField.reset(new volScalarField(tmrtHeader, mesh));
                if (!tgField.valid())
                {
                    tgSource = tmrtFieldName;
                }
            }

            IOobject gHeader
            (
                "G",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            );
            if (gHeader.typeHeaderOk<volScalarField>())
            {
                gField.reset(new volScalarField(gHeader, mesh));
                if (!tgField.valid() && !tmrtField.valid())
                {
                    tgSource = "G";
                }
            }

            IOobject qrHeader
            (
                "qr",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            );
            if (qrHeader.typeHeaderOk<volScalarField>())
            {
                qrField.reset(new volScalarField(qrHeader, mesh));
                if (!tgField.valid() && !tmrtField.valid() && !gField.valid())
                {
                    tgSource = "qr";
                }
            }

            IOobject iDefaultHeader
            (
                "IDefault",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            );
            if (iDefaultHeader.typeHeaderOk<volScalarField>())
            {
                iDefaultField.reset(new volScalarField(iDefaultHeader, mesh));
                if
                (
                    !tgField.valid()
                 && !tmrtField.valid()
                 && !gField.valid()
                 && !qrField.valid()
                )
                {
                    tgSource = "IDefault";
                }
            }
        }

        const scalar wallAverageTemperatureC = calculateWallAverageTemperatureC(mesh);
        const labelList cellsToAnalyze = getCellsToAnalyze(mesh, args);

        volScalarField tNwbOut
        (
            IOobject
            (
                "TnwbISO7243_C",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0),
            "zeroGradient"
        );

        volScalarField tgOut
        (
            IOobject
            (
                "TgISO7243_C",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0),
            "zeroGradient"
        );

        volScalarField wbgtIndoorOut
        (
            IOobject
            (
                "WBGTIndoorISO7243_C",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0),
            "zeroGradient"
        );

        volScalarField wbgtOutdoorOut
        (
            IOobject
            (
                "WBGTOutdoorISO7243_C",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0),
            "zeroGradient"
        );

        volScalarField wbgtSelectedOut
        (
            IOobject
            (
                "WBGTISO7243_C",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0),
            "zeroGradient"
        );

        volScalarField riskLevelOut
        (
            IOobject
            (
                "ISO7243RiskLevel",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0),
            "zeroGradient"
        );

        scalar volumeWeightedTaC(0.0);
        scalar volumeWeightedRh(0.0);
        scalar volumeWeightedTnwbC(0.0);
        scalar volumeWeightedTgC(0.0);
        scalar volumeWeightedWbgtC(0.0);
        scalar totalVolume(0.0);
        scalar minWbgt(GREAT);
        scalar maxWbgt(-GREAT);

        label riskCount0(0);
        label riskCount1(0);
        label riskCount2(0);
        label riskCount3(0);

        forAll(cellsToAnalyze, i)
        {
            const label cellI = cellsToAnalyze[i];

            const scalar taK = T[cellI];
            const scalar taC = taK - ISO7243Constants::kelvinOffset;

            scalar rh = defaultRH;
            if (rhField.valid())
            {
                const scalar rawRh = rhField()[cellI];
                rh = (rawRh <= 1.2) ? rawRh*100.0 : rawRh;
            }
            rh = clampRH(rh);

            scalar tNwbC = 0.0;
            if (tNwbField.valid())
            {
                tNwbC = toCelsiusAuto(tNwbField()[cellI]);
            }
            else
            {
                tNwbC = estimateNaturalWetBulbStull(taC, rh);
            }

            scalar tgC = 0.0;
            if (tgField.valid())
            {
                tgC = toCelsiusAuto(tgField()[cellI]);
            }
            else if (tmrtField.valid())
            {
                tgC = toCelsiusAuto(tmrtField()[cellI]) + globeOffsetC;
            }
            else if (gField.valid())
            {
                tgC = estimateRadiantTemperatureFromG
                (
                    Foam::max(0.0, gField()[cellI])
                ) + globeOffsetC;
            }
            else if (iDefaultField.valid())
            {
                const scalar gFromIDefault =
                    4.0*constant::mathematical::pi*Foam::max(0.0, iDefaultField()[cellI]);
                tgC = estimateRadiantTemperatureFromG(gFromIDefault) + globeOffsetC;
            }
            else if (qrField.valid())
            {
                tgC = estimateRadiantTemperatureFromQr
                (
                    qrField()[cellI],
                    taK,
                    globeEmissivity
                ) + globeOffsetC;
            }
            else
            {
                tgC = wallAverageTemperatureC + globeOffsetC;
            }

            const scalar wbgtIndoorC = 0.7*tNwbC + 0.3*tgC;
            const scalar wbgtOutdoorC = 0.7*tNwbC + 0.2*tgC + 0.1*taC;
            const scalar wbgtC = useOutdoorFormula ? wbgtOutdoorC : wbgtIndoorC;

            const scalar riskLevel = computeRiskLevel
            (
                wbgtC,
                lowRiskLimit,
                moderateRiskLimit,
                highRiskLimit
            );

            tNwbOut[cellI] = tNwbC;
            tgOut[cellI] = tgC;
            wbgtIndoorOut[cellI] = wbgtIndoorC;
            wbgtOutdoorOut[cellI] = wbgtOutdoorC;
            wbgtSelectedOut[cellI] = wbgtC;
            riskLevelOut[cellI] = riskLevel;

            if (riskLevel < 0.5)
            {
                riskCount0++;
            }
            else if (riskLevel < 1.5)
            {
                riskCount1++;
            }
            else if (riskLevel < 2.5)
            {
                riskCount2++;
            }
            else
            {
                riskCount3++;
            }

            minWbgt = Foam::min(minWbgt, wbgtC);
            maxWbgt = Foam::max(maxWbgt, wbgtC);

            const scalar cellVolume = mesh.V()[cellI];
            volumeWeightedTaC += taC*cellVolume;
            volumeWeightedRh += rh*cellVolume;
            volumeWeightedTnwbC += tNwbC*cellVolume;
            volumeWeightedTgC += tgC*cellVolume;
            volumeWeightedWbgtC += wbgtC*cellVolume;
            totalVolume += cellVolume;
        }

        tNwbOut.write();
        tgOut.write();
        wbgtIndoorOut.write();
        wbgtOutdoorOut.write();
        wbgtSelectedOut.write();
        riskLevelOut.write();

        reduce(volumeWeightedTaC, sumOp<scalar>());
        reduce(volumeWeightedRh, sumOp<scalar>());
        reduce(volumeWeightedTnwbC, sumOp<scalar>());
        reduce(volumeWeightedTgC, sumOp<scalar>());
        reduce(volumeWeightedWbgtC, sumOp<scalar>());
        reduce(totalVolume, sumOp<scalar>());
        reduce(minWbgt, minOp<scalar>());
        reduce(maxWbgt, maxOp<scalar>());
        reduce(riskCount0, sumOp<label>());
        reduce(riskCount1, sumOp<label>());
        reduce(riskCount2, sumOp<label>());
        reduce(riskCount3, sumOp<label>());

        const scalar avgTaC = volumeWeightedTaC/Foam::max(SMALL, totalVolume);
        const scalar avgRh = volumeWeightedRh/Foam::max(SMALL, totalVolume);
        const scalar avgTnwbC = volumeWeightedTnwbC/Foam::max(SMALL, totalVolume);
        const scalar avgTgC = volumeWeightedTgC/Foam::max(SMALL, totalVolume);
        const scalar avgWbgtC = volumeWeightedWbgtC/Foam::max(SMALL, totalVolume);

        const label analyzedCells = returnReduce(cellsToAnalyze.size(), sumOp<label>());

        Info<< nl
            << "=============== ISO 7243 WBGT ANALYSIS ===============" << nl
            << "Formula: "
            << (useOutdoorFormula
                ? "Outdoor/solar WBGT = 0.7*Tnwb + 0.2*Tg + 0.1*Ta"
                : "Indoor/no-solar WBGT = 0.7*Tnwb + 0.3*Tg")
            << nl
            << "Analyzed cells:                 " << analyzedCells << nl
            << "Analysis volume:                " << totalVolume << " m^3" << nl
            << "Humidity source:                " << rhSource << nl
            << "Tnwb source:                    " << tNwbSource << nl
            << "Tg source:                      " << tgSource << nl
            << nl
            << "Average Ta:                     " << avgTaC << " degC" << nl
            << "Average RH:                     " << avgRh << " %" << nl
            << "Average Tnwb:                   " << avgTnwbC << " degC" << nl
            << "Average Tg:                     " << avgTgC << " degC" << nl
            << "Average WBGT:                   " << avgWbgtC << " degC" << nl
            << "Min/Max WBGT:                   " << minWbgt
            << " / " << maxWbgt << " degC" << nl
            << nl
            << "Risk thresholds (degC):         "
            << lowRiskLimit << " / "
            << moderateRiskLimit << " / "
            << highRiskLimit << nl
            << "Risk level counts (0/1/2/3):   "
            << riskCount0 << " / "
            << riskCount1 << " / "
            << riskCount2 << " / "
            << riskCount3 << nl
            << "======================================================" << endl;
    }

    return 0;
}

// ************************************************************************* //
