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
    ISO7933Foam

Description
    Dynamic PHS-oriented post-processing utility for ISO 7933 workflows.

    The model uses transient integration of core/skin temperatures, sweat
    response and dehydration over exposure time, using CFD fields as boundary
    conditions. It is intended for engineering assessment and hotspot ranking.

    Output fields:
    - TaISO7933_C                Air temperature [degC]
    - TrISO7933_C                Mean radiant temperature [degC]
    - PaISO7933_kPa              Water vapor partial pressure [kPa]
    - ISO7933Ereq_Wm2            Peak required evaporative heat loss [W/m^2]
    - ISO7933Emax_Wm2            Minimum evaporative capacity [W/m^2]
    - ISO7933Wreq                Peak required skin wettedness [-]
    - ISO7933ReqSweatRate_gm2h   Peak required sweat rate [g/(m^2 h)]
    - ISO7933SweatRate_gm2h      Peak simulated sweat rate [g/(m^2 h)]
    - ISO7933Eact_Wm2            Peak effective evaporation [W/m^2]
    - ISO7933Tre_C               Final core temperature [degC]
    - ISO7933Tsk_C               Final skin temperature [degC]
    - ISO7933PredWaterLoss_g     Predicted cumulative water loss [g]
    - ISO7933DlimTre_min         Time to core temperature limit [min]
    - ISO7933DlimLoss95_min      Time to 95% water loss criterion [min]
    - ISO7933DlimLoss50_min      Time to 50% water loss criterion [min]
    - ISO7933HeatStrainIndex     Combined dynamic heat-strain index [-]
    - ISO7933RiskLevel           0..3 screening level

    Notes:
    - This is a dynamic PHS-oriented engineering implementation, not a
      certified normative ISO 7933 compliance engine.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "cellSet.H"
#include "wallFvPatch.H"
#include "mathematicalConstants.H"

namespace Foam
{
namespace ISO7933Constants
{
    const scalar stefanBoltzmann = 5.670374419e-8;  // W/(m^2 K^4)
    const scalar kelvinOffset = 273.15;
    const scalar metToWm2 = 58.15;
    const scalar latentHeatEvapJPerG = 2426.0;
    const scalar emissivityHuman = 0.95;
}
}

struct PHSParameters
{
    Foam::scalar M;
    Foam::scalar W;
    Foam::scalar clo;
    Foam::scalar im;

    Foam::scalar exposureDurationMin;
    Foam::scalar timeStepMin;

    Foam::scalar bodySurfaceArea;
    Foam::scalar bodyMassKg;
    Foam::scalar bodySpecificHeat;

    Foam::scalar initialCoreTempC;
    Foam::scalar initialSkinTempC;
    Foam::scalar coreNeutralTempC;
    Foam::scalar skinNeutralTempC;
    Foam::scalar coreTempLimitC;

    Foam::scalar coreCompartmentFraction;
    Foam::scalar coreSkinConductance;

    Foam::scalar maxSweatRate;
    Foam::scalar sweatControlGain;
    Foam::scalar sweatTimeConstantMin;

    Foam::scalar lossLimit95Fraction;
    Foam::scalar lossLimit50Fraction;

    Foam::scalar targetStorageWm2;
    Foam::scalar evaporationEfficiency;

    Foam::scalar radiantOffsetC;
    Foam::scalar globeEmissivity;
};

struct CellEnvironment
{
    Foam::scalar taC;
    Foam::scalar trC;
    Foam::scalar paKPa;
    Foam::scalar va;
};

struct HeatExchange
{
    Foam::scalar dryLossWm2;
    Foam::scalar hc;
    Foam::scalar hr;
};

struct PHSResult
{
    Foam::scalar taC;
    Foam::scalar trC;
    Foam::scalar paKPa;

    Foam::scalar ereqPeakWm2;
    Foam::scalar emaxMinWm2;
    Foam::scalar wreqPeak;
    Foam::scalar reqSweatRatePeak;
    Foam::scalar sweatRatePeak;
    Foam::scalar eactPeakWm2;

    Foam::scalar treFinalC;
    Foam::scalar tskFinalC;
    Foam::scalar waterLossG;

    Foam::scalar dlimTreMin;
    Foam::scalar dlimLoss95Min;
    Foam::scalar dlimLoss50Min;

    Foam::scalar heatStrainIndex;
    Foam::scalar riskLevel;
};

Foam::scalar clampRH(const Foam::scalar rh)
{
    return Foam::max(0.0, Foam::min(100.0, rh));
}

Foam::scalar clampScalar
(
    const Foam::scalar value,
    const Foam::scalar minValue,
    const Foam::scalar maxValue
)
{
    return Foam::max(minValue, Foam::min(maxValue, value));
}

Foam::scalar toCelsiusAuto(const Foam::scalar value)
{
    return (value > 170.0) ? value - ISO7933Constants::kelvinOffset : value;
}

Foam::scalar saturationVapourPressureKPa(const Foam::scalar tC)
{
    // Tetens approximation
    return 0.61078*Foam::exp((17.2694*tC)/(tC + 237.3));
}

Foam::scalar estimateRadiantTemperatureFromG(const Foam::scalar gWm2)
{
    const Foam::scalar g = Foam::max(SMALL, gWm2);

    return Foam::pow(g/(4.0*ISO7933Constants::stefanBoltzmann), 0.25)
         - ISO7933Constants::kelvinOffset;
}

Foam::scalar estimateRadiantTemperatureFromQr
(
    const Foam::scalar qrWm2,
    const Foam::scalar taK,
    const Foam::scalar emissivity
)
{
    const Foam::scalar eps = clampScalar(emissivity, 0.01, 1.0);
    const Foam::scalar denom = Foam::max(SMALL, eps*ISO7933Constants::stefanBoltzmann);

    Foam::scalar trad4 = Foam::pow(taK, 4) + qrWm2/denom;
    trad4 = Foam::max(SMALL, trad4);

    return Foam::pow(trad4, 0.25) - ISO7933Constants::kelvinOffset;
}

Foam::scalar calculateWallAverageTemperatureC(const Foam::fvMesh& mesh)
{
    const Foam::volScalarField& T = mesh.lookupObject<Foam::volScalarField>("T");

    Foam::scalar weightedWallTemp(0.0);
    Foam::scalar totalWallArea(0.0);

    const Foam::fvPatchList& patches = mesh.boundary();
    forAll(patches, patchI)
    {
        if (!Foam::isType<Foam::wallFvPatch>(patches[patchI]))
        {
            continue;
        }

        const Foam::label pI = patches[patchI].index();
        const Foam::scalar area = Foam::gSum(mesh.magSf().boundaryField()[pI]);

        if (area > SMALL)
        {
            weightedWallTemp += Foam::gSum
            (
                mesh.magSf().boundaryField()[pI] * T.boundaryField()[pI]
            );
            totalWallArea += area;
        }
    }

    if (totalWallArea > SMALL)
    {
        return (weightedWallTemp/totalWallArea) - ISO7933Constants::kelvinOffset;
    }

    Foam::scalar weightedCellTemp(0.0);
    Foam::scalar totalVolume(0.0);
    forAll(mesh.C(), cellI)
    {
        weightedCellTemp += mesh.V()[cellI]*T[cellI];
        totalVolume += mesh.V()[cellI];
    }

    Foam::reduce(weightedCellTemp, Foam::sumOp<Foam::scalar>());
    Foam::reduce(totalVolume, Foam::sumOp<Foam::scalar>());

    return (weightedCellTemp/Foam::max(SMALL, totalVolume))
         - ISO7933Constants::kelvinOffset;
}

Foam::labelList getCellsToAnalyze(const Foam::fvMesh& mesh, const Foam::argList& args)
{
    Foam::labelList cellsToAnalyze;

    if (args.found("cellSet"))
    {
        const Foam::word cellSetName = args.get<Foam::word>("cellSet");
        Foam::Info<< "Loading cellSet: " << cellSetName << Foam::endl;

        Foam::cellSet selectedCells(mesh, cellSetName);
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
        const Foam::word cellZoneName = args.get<Foam::word>("cellZone");
        Foam::Info<< "Loading cellZone: " << cellZoneName << Foam::endl;

        const Foam::cellZoneMesh& cellZones = mesh.cellZones();
        const Foam::label zoneId = cellZones.findZoneID(cellZoneName);

        if (zoneId == -1)
        {
            FatalErrorInFunction
                << "cellZone " << cellZoneName << " does not exist" << Foam::nl
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

    const Foam::label nCells = Foam::returnReduce(cellsToAnalyze.size(), Foam::sumOp<Foam::label>());
    Foam::Info<< "Analyzing " << nCells << " cells" << Foam::endl;

    return cellsToAnalyze;
}

HeatExchange computeDryHeatLoss
(
    const CellEnvironment& env,
    const PHSParameters& params,
    const Foam::scalar tskC,
    const Foam::scalar fcl
)
{
    HeatExchange hx;

    const Foam::scalar hcf = 12.1*Foam::sqrt(Foam::max(0.01, env.va));
    const Foam::scalar hcn = 2.38*Foam::pow(Foam::mag(tskC - env.taC), 0.25);
    hx.hc = Foam::max(2.0, Foam::max(hcf, hcn));

    const Foam::scalar tskK = tskC + ISO7933Constants::kelvinOffset;
    const Foam::scalar trK = env.trC + ISO7933Constants::kelvinOffset;
    const Foam::scalar tMeanK = 0.5*(tskK + trK);

    hx.hr = 4.0*ISO7933Constants::emissivityHuman
          * ISO7933Constants::stefanBoltzmann
          * Foam::pow(Foam::max(150.0, tMeanK), 3);

    const Foam::scalar combinedH = Foam::max(SMALL, hx.hc + hx.hr);
    const Foam::scalar toC = (hx.hc*env.taC + hx.hr*env.trC)/combinedH;

    hx.dryLossWm2 = fcl*combinedH*(tskC - toC);

    return hx;
}

void computeRespiratoryLosses
(
    const Foam::scalar M,
    const CellEnvironment& env,
    Foam::scalar& cres,
    Foam::scalar& eres
)
{
    cres = 0.0014*M*(34.0 - env.taC);
    eres = 0.0173*M*(5.87 - env.paKPa);
}

Foam::scalar computeEvaporativeCapacity
(
    const CellEnvironment& env,
    const PHSParameters& params,
    const Foam::scalar tskC,
    const Foam::scalar Icl,
    const Foam::scalar hc
)
{
    const Foam::scalar pskKPa = saturationVapourPressureKPa(tskC);

    const Foam::scalar rea = 1.0/(16.5*Foam::max(2.0, hc));
    const Foam::scalar recl = Icl/(Foam::max(0.05, params.im)*16.5);

    Foam::scalar emax = (pskKPa - env.paKPa)/Foam::max(SMALL, rea + recl);
    emax = Foam::max(0.0, emax);

    return params.evaporationEfficiency*emax;
}

Foam::scalar computeRiskLevel
(
    const Foam::scalar ereqPeak,
    const Foam::scalar heatStrainIndex,
    const Foam::scalar treRatio
)
{
    if (ereqPeak <= SMALL && treRatio < 0.2)
    {
        return 0.0;
    }

    if (heatStrainIndex < 0.75)
    {
        return 1.0;
    }

    if (heatStrainIndex < 1.0)
    {
        return 2.0;
    }

    return 3.0;
}

PHSResult simulatePHS
(
    const CellEnvironment& env,
    const PHSParameters& params
)
{
    PHSResult result;

    result.taC = env.taC;
    result.trC = env.trC;
    result.paKPa = env.paKPa;

    result.ereqPeakWm2 = 0.0;
    result.emaxMinWm2 = GREAT;
    result.wreqPeak = 0.0;
    result.reqSweatRatePeak = 0.0;
    result.sweatRatePeak = 0.0;
    result.eactPeakWm2 = 0.0;

    result.treFinalC = params.initialCoreTempC;
    result.tskFinalC = params.initialSkinTempC;
    result.waterLossG = 0.0;

    result.dlimTreMin = params.exposureDurationMin;
    result.dlimLoss95Min = params.exposureDurationMin;
    result.dlimLoss50Min = params.exposureDurationMin;

    result.heatStrainIndex = 0.0;
    result.riskLevel = 0.0;

    const Foam::scalar Icl = 0.155*params.clo;
    const Foam::scalar fcl = (Icl < 0.078)
        ? (1.0 + 1.29*Icl)
        : (1.05 + 0.645*Icl);

    const Foam::scalar dtMinRaw = Foam::max(0.1, params.timeStepMin);
    Foam::label nSteps = Foam::label(params.exposureDurationMin/dtMinRaw + 0.999999);
    nSteps = Foam::max(Foam::label(1), nSteps);

    const Foam::scalar dtMin = params.exposureDurationMin/static_cast<Foam::scalar>(nSteps);
    const Foam::scalar dtSec = 60.0*dtMin;
    const Foam::scalar dtHour = dtMin/60.0;

    const Foam::scalar coreHeatCapacity = Foam::max
    (
        1.0,
        params.bodyMassKg*params.bodySpecificHeat*params.coreCompartmentFraction
    );
    const Foam::scalar skinHeatCapacity = Foam::max
    (
        1.0,
        params.bodyMassKg*params.bodySpecificHeat*(1.0 - params.coreCompartmentFraction)
    );

    const Foam::scalar lossLimit95G =
        params.bodyMassKg*1000.0*params.lossLimit95Fraction;
    const Foam::scalar lossLimit50G =
        params.bodyMassKg*1000.0*params.lossLimit50Fraction;

    Foam::scalar treC = params.initialCoreTempC;
    Foam::scalar tskC = params.initialSkinTempC;
    Foam::scalar sweatRateGm2h = 0.0;
    Foam::scalar waterLossG = 0.0;

    bool treLimitReached = false;
    bool loss95Reached = false;
    bool loss50Reached = false;

    for (Foam::label stepI = 0; stepI < nSteps; ++stepI)
    {
        const HeatExchange hx = computeDryHeatLoss(env, params, tskC, fcl);

        Foam::scalar cres(0.0), eres(0.0);
        computeRespiratoryLosses(params.M, env, cres, eres);

        const Foam::scalar ereq =
            params.M - params.W - hx.dryLossWm2 - cres - eres - params.targetStorageWm2;
        const Foam::scalar ereqPos = Foam::max(0.0, ereq);

        const Foam::scalar emax = Foam::max
        (
            SMALL,
            computeEvaporativeCapacity(env, params, tskC, Icl, hx.hc)
        );

        const Foam::scalar wreq = ereqPos/Foam::max(SMALL, emax);
        const Foam::scalar reqSweatRate =
            ereqPos*3600.0/ISO7933Constants::latentHeatEvapJPerG;

        // Thermoregulatory drive combines required evaporation with core/skin rise.
        const Foam::scalar thermalDrive =
            Foam::max(0.0, treC - params.coreNeutralTempC)
          + 0.2*Foam::max(0.0, tskC - params.skinNeutralTempC);

        const Foam::scalar controlSweatRate = params.sweatControlGain*thermalDrive;

        const Foam::scalar desiredSweatRate = Foam::min
        (
            params.maxSweatRate,
            Foam::max(reqSweatRate, controlSweatRate)
        );

        const Foam::scalar alpha =
            1.0 - Foam::exp(-dtMin/Foam::max(SMALL, params.sweatTimeConstantMin));
        sweatRateGm2h += alpha*(desiredSweatRate - sweatRateGm2h);
        sweatRateGm2h = clampScalar(sweatRateGm2h, 0.0, params.maxSweatRate);

        const Foam::scalar evapPotential =
            sweatRateGm2h*ISO7933Constants::latentHeatEvapJPerG/3600.0;
        const Foam::scalar eact = Foam::min(emax, evapPotential);

        const Foam::scalar sweatEvapRate =
            eact*3600.0/ISO7933Constants::latentHeatEvapJPerG;
        const Foam::scalar respWaterRate =
            Foam::max(0.0, eres)*3600.0/ISO7933Constants::latentHeatEvapJPerG;

        waterLossG +=
            (sweatEvapRate + respWaterRate)*Foam::max(0.1, params.bodySurfaceArea)*dtHour;

        // Two-compartment transient energy balance.
        const Foam::scalar metabolicW =
            (params.M - params.W)*Foam::max(0.1, params.bodySurfaceArea);
        const Foam::scalar respW =
            (cres + eres)*Foam::max(0.1, params.bodySurfaceArea);
        const Foam::scalar coreToSkinW =
            params.coreSkinConductance*Foam::max(0.1, params.bodySurfaceArea)*(treC - tskC);
        const Foam::scalar envSkinLossW =
            (hx.dryLossWm2 + eact)*Foam::max(0.1, params.bodySurfaceArea);

        const Foam::scalar qCore = metabolicW - respW - coreToSkinW;
        const Foam::scalar qSkin = coreToSkinW - envSkinLossW;

        treC += (qCore/coreHeatCapacity)*dtSec;
        tskC += (qSkin/skinHeatCapacity)*dtSec;

        treC = clampScalar(treC, 35.0, 43.0);
        tskC = clampScalar(tskC, 20.0, 42.0);

        const Foam::scalar currentTimeMin = (stepI + 1)*dtMin;

        if (!treLimitReached && treC >= params.coreTempLimitC)
        {
            result.dlimTreMin = currentTimeMin;
            treLimitReached = true;
        }

        if (!loss95Reached && waterLossG >= lossLimit95G)
        {
            result.dlimLoss95Min = currentTimeMin;
            loss95Reached = true;
        }

        if (!loss50Reached && waterLossG >= lossLimit50G)
        {
            result.dlimLoss50Min = currentTimeMin;
            loss50Reached = true;
        }

        result.ereqPeakWm2 = Foam::max(result.ereqPeakWm2, ereqPos);
        result.emaxMinWm2 = Foam::min(result.emaxMinWm2, emax);
        result.wreqPeak = Foam::max(result.wreqPeak, wreq);
        result.reqSweatRatePeak = Foam::max(result.reqSweatRatePeak, reqSweatRate);
        result.sweatRatePeak = Foam::max(result.sweatRatePeak, sweatRateGm2h);
        result.eactPeakWm2 = Foam::max(result.eactPeakWm2, eact);
    }

    if (result.emaxMinWm2 == GREAT)
    {
        result.emaxMinWm2 = 0.0;
    }

    result.treFinalC = treC;
    result.tskFinalC = tskC;
    result.waterLossG = waterLossG;

    const Foam::scalar treDenom =
        Foam::max(SMALL, params.coreTempLimitC - params.initialCoreTempC);
    const Foam::scalar treRatio =
        Foam::max(0.0, (result.treFinalC - params.initialCoreTempC)/treDenom);

    const Foam::scalar loss95Ratio =
        result.waterLossG/Foam::max(SMALL, lossLimit95G);
    const Foam::scalar sweatRatio =
        result.sweatRatePeak/Foam::max(50.0, params.maxSweatRate);

    result.heatStrainIndex = Foam::max
    (
        result.wreqPeak,
        Foam::max(treRatio, Foam::max(loss95Ratio, sweatRatio))
    );

    result.riskLevel = computeRiskLevel
    (
        result.ereqPeakWm2,
        result.heatStrainIndex,
        treRatio
    );

    return result;
}

int main(int argc, char *argv[])
{
    Foam::timeSelector::addOptions();

    Foam::argList::addOption
    (
        "cellSet",
        "name",
        "Analyze only cells in the specified cellSet"
    );

    Foam::argList::addOption
    (
        "cellZone",
        "name",
        "Analyze only cells in the specified cellZone"
    );

    #include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"

    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);
        Foam::Info<< "\nTime = " << runTime.timeName() << Foam::endl;

        Foam::IOobject dictHeader
        (
            "ISO7933Dict",
            runTime.time().constant(),
            runTime,
            Foam::IOobject::READ_IF_PRESENT,
            Foam::IOobject::NO_WRITE
        );

        Foam::autoPtr<Foam::IOdictionary> dictPtr;
        if (dictHeader.typeHeaderOk<Foam::IOdictionary>())
        {
            dictPtr.reset(new Foam::IOdictionary(dictHeader));
            Foam::Info<< "Using settings from constant/ISO7933Dict" << Foam::endl;
        }
        else
        {
            Foam::Info<< "constant/ISO7933Dict not found: using internal defaults" << Foam::endl;
        }

        const Foam::scalar defaultRH = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("RH", 50.0)
            : 50.0;

        const Foam::scalar metInput = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("met", 1.6)
            : 1.6;
        const Foam::scalar wmeInput = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("wme", 0.0)
            : 0.0;
        const bool metIsMetUnit = dictPtr.valid()
            ? dictPtr().lookupOrDefault<bool>("metIsMetUnit", true)
            : true;

        PHSParameters params;
        params.M = metIsMetUnit ? metInput*ISO7933Constants::metToWm2 : metInput;
        params.W = metIsMetUnit ? wmeInput*ISO7933Constants::metToWm2 : wmeInput;

        params.clo = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("clo", 0.5)
            : 0.5;
        params.im = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("im", 0.38)
            : 0.38;

        params.exposureDurationMin = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("exposureDuration", 480.0)
            : 480.0;
        params.timeStepMin = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("timeStep", 1.0)
            : 1.0;

        params.bodySurfaceArea = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("bodySurfaceArea", 1.8)
            : 1.8;
        params.bodyMassKg = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("bodyMass", 75.0)
            : 75.0;
        params.bodySpecificHeat = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("bodySpecificHeat", 3470.0)
            : 3470.0;

        const Foam::scalar legacySkinTemp = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("skinTemperature", 34.1)
            : 34.1;
        const Foam::scalar legacyLossFraction = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("maxWaterLossFraction", 0.05)
            : 0.05;

        params.initialCoreTempC = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("initialCoreTemperature", 37.0)
            : 37.0;
        params.initialSkinTempC = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("initialSkinTemperature", legacySkinTemp)
            : legacySkinTemp;
        params.coreNeutralTempC = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("coreNeutralTemperature", 37.0)
            : 37.0;
        params.skinNeutralTempC = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("skinNeutralTemperature", 34.1)
            : 34.1;
        params.coreTempLimitC = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("coreTemperatureLimit", 38.0)
            : 38.0;

        params.coreCompartmentFraction = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("coreCompartmentFraction", 0.9)
            : 0.9;
        params.coreSkinConductance = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("coreSkinConductance", 5.3)
            : 5.3;

        params.maxSweatRate = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("maxSweatRate", 500.0)
            : 500.0;
        params.sweatControlGain = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("sweatControlGain", 350.0)
            : 350.0;
        params.sweatTimeConstantMin = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("sweatTimeConstant", 10.0)
            : 10.0;

        params.lossLimit95Fraction = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("lossLimit95Fraction", legacyLossFraction)
            : legacyLossFraction;
        params.lossLimit50Fraction = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("lossLimit50Fraction", 0.075)
            : 0.075;

        params.targetStorageWm2 = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("targetStorage", 0.0)
            : 0.0;
        params.evaporationEfficiency = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("evaporationEfficiency", 1.0)
            : 1.0;

        params.radiantOffsetC = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("radiantTemperatureOffset", 0.0)
            : 0.0;
        params.globeEmissivity = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::scalar>("globeEmissivity", 0.95)
            : 0.95;

        const Foam::word trFieldName = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::word>("trFieldName", "Tr")
            : "Tr";
        const Foam::word tmrtFieldName = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::word>("tmrtFieldName", "Tmrt")
            : "Tmrt";
        const Foam::word tgFieldName = dictPtr.valid()
            ? dictPtr().lookupOrDefault<Foam::word>("tgFieldName", "Tg")
            : "Tg";

        params.evaporationEfficiency = clampScalar(params.evaporationEfficiency, 0.0, 1.0);
        params.coreCompartmentFraction = clampScalar(params.coreCompartmentFraction, 0.6, 0.98);

        if (params.lossLimit50Fraction < params.lossLimit95Fraction)
        {
            Foam::Info
                << "Swapping lossLimit50Fraction/lossLimit95Fraction to keep 50% >= 95% criterion"
                << Foam::endl;
            const Foam::scalar tmp = params.lossLimit50Fraction;
            params.lossLimit50Fraction = params.lossLimit95Fraction;
            params.lossLimit95Fraction = tmp;
        }

        if
        (
            params.M <= 0.0
         || params.W < 0.0
         || params.clo < 0.0
         || params.im <= 0.0
         || params.exposureDurationMin <= 0.0
         || params.timeStepMin <= 0.0
         || params.bodySurfaceArea <= 0.0
         || params.bodyMassKg <= 0.0
        )
        {
            FatalErrorInFunction
                << "Invalid ISO7933Dict parameters: check met/wme/clo/im/time/body values"
                << abort(FatalError);
        }

        Foam::IOobject THeader
        (
            "T",
            runTime.timeName(),
            mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::NO_WRITE
        );
        Foam::volScalarField T(THeader, mesh);

        Foam::IOobject UHeader
        (
            "U",
            runTime.timeName(),
            mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::NO_WRITE
        );
        Foam::volVectorField U(UHeader, mesh);

        Foam::autoPtr<Foam::volScalarField> rhField;
        Foam::word rhSource("ISO7933Dict RH");
        {
            const Foam::wordList rhCandidates
            ({
                "thermo:relHum",
                "thermoRelHum",
                "relHum",
                "RH",
                "relativeHumidity"
            });

            forAll(rhCandidates, i)
            {
                Foam::IOobject rhHeader
                (
                    rhCandidates[i],
                    runTime.timeName(),
                    mesh,
                    Foam::IOobject::READ_IF_PRESENT,
                    Foam::IOobject::NO_WRITE
                );

                if (rhHeader.typeHeaderOk<Foam::volScalarField>())
                {
                    rhField.reset(new Foam::volScalarField(rhHeader, mesh));
                    rhSource = rhCandidates[i];
                    break;
                }
            }
        }

        Foam::autoPtr<Foam::volScalarField> trField;
        Foam::autoPtr<Foam::volScalarField> tmrtField;
        Foam::autoPtr<Foam::volScalarField> tgField;
        Foam::autoPtr<Foam::volScalarField> gField;
        Foam::autoPtr<Foam::volScalarField> qrField;
        Foam::autoPtr<Foam::volScalarField> iDefaultField;

        Foam::word trSource("wallTemperature fallback");
        {
            Foam::IOobject trHeader
            (
                trFieldName,
                runTime.timeName(),
                mesh,
                Foam::IOobject::READ_IF_PRESENT,
                Foam::IOobject::NO_WRITE
            );
            if (trHeader.typeHeaderOk<Foam::volScalarField>())
            {
                trField.reset(new Foam::volScalarField(trHeader, mesh));
                trSource = trFieldName;
            }

            Foam::IOobject tmrtHeader
            (
                tmrtFieldName,
                runTime.timeName(),
                mesh,
                Foam::IOobject::READ_IF_PRESENT,
                Foam::IOobject::NO_WRITE
            );
            if (tmrtHeader.typeHeaderOk<Foam::volScalarField>())
            {
                tmrtField.reset(new Foam::volScalarField(tmrtHeader, mesh));
                if (!trField.valid())
                {
                    trSource = tmrtFieldName;
                }
            }

            Foam::IOobject tgHeader
            (
                tgFieldName,
                runTime.timeName(),
                mesh,
                Foam::IOobject::READ_IF_PRESENT,
                Foam::IOobject::NO_WRITE
            );
            if (tgHeader.typeHeaderOk<Foam::volScalarField>())
            {
                tgField.reset(new Foam::volScalarField(tgHeader, mesh));
                if (!trField.valid() && !tmrtField.valid())
                {
                    trSource = tgFieldName;
                }
            }

            Foam::IOobject gHeader
            (
                "G",
                runTime.timeName(),
                mesh,
                Foam::IOobject::READ_IF_PRESENT,
                Foam::IOobject::NO_WRITE
            );
            if (gHeader.typeHeaderOk<Foam::volScalarField>())
            {
                gField.reset(new Foam::volScalarField(gHeader, mesh));
                if (!trField.valid() && !tmrtField.valid() && !tgField.valid())
                {
                    trSource = "G";
                }
            }

            Foam::IOobject qrHeader
            (
                "qr",
                runTime.timeName(),
                mesh,
                Foam::IOobject::READ_IF_PRESENT,
                Foam::IOobject::NO_WRITE
            );
            if (qrHeader.typeHeaderOk<Foam::volScalarField>())
            {
                qrField.reset(new Foam::volScalarField(qrHeader, mesh));
                if
                (
                    !trField.valid()
                 && !tmrtField.valid()
                 && !tgField.valid()
                 && !gField.valid()
                )
                {
                    trSource = "qr";
                }
            }

            Foam::IOobject iDefaultHeader
            (
                "IDefault",
                runTime.timeName(),
                mesh,
                Foam::IOobject::READ_IF_PRESENT,
                Foam::IOobject::NO_WRITE
            );
            if (iDefaultHeader.typeHeaderOk<Foam::volScalarField>())
            {
                iDefaultField.reset(new Foam::volScalarField(iDefaultHeader, mesh));
                if
                (
                    !trField.valid()
                 && !tmrtField.valid()
                 && !tgField.valid()
                 && !gField.valid()
                 && !qrField.valid()
                )
                {
                    trSource = "IDefault";
                }
            }
        }

        const Foam::scalar wallAverageTempC = calculateWallAverageTemperatureC(mesh);
        const Foam::labelList cellsToAnalyze = getCellsToAnalyze(mesh, args);

        Foam::volScalarField taOut
        (
            Foam::IOobject
            (
                "TaISO7933_C",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("zero", Foam::dimless, 0.0),
            "zeroGradient"
        );

        Foam::volScalarField trOut
        (
            Foam::IOobject
            (
                "TrISO7933_C",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("zero", Foam::dimless, 0.0),
            "zeroGradient"
        );

        Foam::volScalarField paOut
        (
            Foam::IOobject
            (
                "PaISO7933_kPa",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("zero", Foam::dimless, 0.0),
            "zeroGradient"
        );

        Foam::volScalarField ereqOut
        (
            Foam::IOobject
            (
                "ISO7933Ereq_Wm2",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("zero", Foam::dimless, 0.0),
            "zeroGradient"
        );

        Foam::volScalarField emaxOut
        (
            Foam::IOobject
            (
                "ISO7933Emax_Wm2",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("zero", Foam::dimless, 0.0),
            "zeroGradient"
        );

        Foam::volScalarField wreqOut
        (
            Foam::IOobject
            (
                "ISO7933Wreq",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("zero", Foam::dimless, 0.0),
            "zeroGradient"
        );

        Foam::volScalarField reqSweatRateOut
        (
            Foam::IOobject
            (
                "ISO7933ReqSweatRate_gm2h",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("zero", Foam::dimless, 0.0),
            "zeroGradient"
        );

        Foam::volScalarField sweatRateOut
        (
            Foam::IOobject
            (
                "ISO7933SweatRate_gm2h",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("zero", Foam::dimless, 0.0),
            "zeroGradient"
        );

        Foam::volScalarField eactOut
        (
            Foam::IOobject
            (
                "ISO7933Eact_Wm2",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("zero", Foam::dimless, 0.0),
            "zeroGradient"
        );

        Foam::volScalarField treOut
        (
            Foam::IOobject
            (
                "ISO7933Tre_C",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("zero", Foam::dimless, 0.0),
            "zeroGradient"
        );

        Foam::volScalarField tskOut
        (
            Foam::IOobject
            (
                "ISO7933Tsk_C",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("zero", Foam::dimless, 0.0),
            "zeroGradient"
        );

        Foam::volScalarField waterLossOut
        (
            Foam::IOobject
            (
                "ISO7933PredWaterLoss_g",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("zero", Foam::dimless, 0.0),
            "zeroGradient"
        );

        Foam::volScalarField dlimTreOut
        (
            Foam::IOobject
            (
                "ISO7933DlimTre_min",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("zero", Foam::dimless, 0.0),
            "zeroGradient"
        );

        Foam::volScalarField dlimLoss95Out
        (
            Foam::IOobject
            (
                "ISO7933DlimLoss95_min",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("zero", Foam::dimless, 0.0),
            "zeroGradient"
        );

        Foam::volScalarField dlimLoss50Out
        (
            Foam::IOobject
            (
                "ISO7933DlimLoss50_min",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("zero", Foam::dimless, 0.0),
            "zeroGradient"
        );

        Foam::volScalarField heatStrainOut
        (
            Foam::IOobject
            (
                "ISO7933HeatStrainIndex",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("zero", Foam::dimless, 0.0),
            "zeroGradient"
        );

        Foam::volScalarField riskOut
        (
            Foam::IOobject
            (
                "ISO7933RiskLevel",
                runTime.timeName(),
                mesh,
                Foam::IOobject::NO_READ,
                Foam::IOobject::AUTO_WRITE
            ),
            mesh,
            Foam::dimensionedScalar("zero", Foam::dimless, 0.0),
            "zeroGradient"
        );

        Foam::scalar weightedTa(0.0), weightedTr(0.0), weightedPa(0.0);
        Foam::scalar weightedEreq(0.0), weightedEmax(0.0), weightedWreq(0.0);
        Foam::scalar weightedReqSweat(0.0), weightedSweat(0.0), weightedEact(0.0);
        Foam::scalar weightedTre(0.0), weightedTsk(0.0), weightedWaterLoss(0.0);
        Foam::scalar weightedDlimTre(0.0), weightedDlimLoss95(0.0), weightedDlimLoss50(0.0);
        Foam::scalar weightedHSI(0.0);
        Foam::scalar totalVolume(0.0);
        Foam::scalar maxHSI(-GREAT);
        Foam::label risk0(0), risk1(0), risk2(0), risk3(0);

        forAll(cellsToAnalyze, i)
        {
            const Foam::label cellI = cellsToAnalyze[i];

            const Foam::scalar taK = T[cellI];
            const Foam::scalar taC = taK - ISO7933Constants::kelvinOffset;

            Foam::scalar rh = defaultRH;
            if (rhField.valid())
            {
                const Foam::scalar rawRh = rhField()[cellI];
                rh = (rawRh <= 1.2) ? rawRh*100.0 : rawRh;
            }
            rh = clampRH(rh);

            Foam::scalar trC = wallAverageTempC + params.radiantOffsetC;
            if (trField.valid())
            {
                trC = toCelsiusAuto(trField()[cellI]);
            }
            else if (tmrtField.valid())
            {
                trC = toCelsiusAuto(tmrtField()[cellI]) + params.radiantOffsetC;
            }
            else if (tgField.valid())
            {
                trC = toCelsiusAuto(tgField()[cellI]) + params.radiantOffsetC;
            }
            else if (gField.valid())
            {
                trC = estimateRadiantTemperatureFromG(Foam::max(0.0, gField()[cellI]))
                    + params.radiantOffsetC;
            }
            else if (iDefaultField.valid())
            {
                const Foam::scalar gFromIDefault =
                    4.0*Foam::constant::mathematical::pi*Foam::max(0.0, iDefaultField()[cellI]);
                trC = estimateRadiantTemperatureFromG(gFromIDefault) + params.radiantOffsetC;
            }
            else if (qrField.valid())
            {
                trC = estimateRadiantTemperatureFromQr
                (
                    qrField()[cellI],
                    taK,
                    params.globeEmissivity
                ) + params.radiantOffsetC;
            }

            CellEnvironment env;
            env.taC = taC;
            env.trC = trC;
            env.paKPa = (rh/100.0)*saturationVapourPressureKPa(taC);
            env.va = Foam::max(0.01, Foam::mag(U[cellI]));

            const PHSResult result = simulatePHS(env, params);

            taOut[cellI] = result.taC;
            trOut[cellI] = result.trC;
            paOut[cellI] = result.paKPa;
            ereqOut[cellI] = result.ereqPeakWm2;
            emaxOut[cellI] = result.emaxMinWm2;
            wreqOut[cellI] = result.wreqPeak;
            reqSweatRateOut[cellI] = result.reqSweatRatePeak;
            sweatRateOut[cellI] = result.sweatRatePeak;
            eactOut[cellI] = result.eactPeakWm2;
            treOut[cellI] = result.treFinalC;
            tskOut[cellI] = result.tskFinalC;
            waterLossOut[cellI] = result.waterLossG;
            dlimTreOut[cellI] = result.dlimTreMin;
            dlimLoss95Out[cellI] = result.dlimLoss95Min;
            dlimLoss50Out[cellI] = result.dlimLoss50Min;
            heatStrainOut[cellI] = result.heatStrainIndex;
            riskOut[cellI] = result.riskLevel;

            if (result.riskLevel < 0.5) risk0++;
            else if (result.riskLevel < 1.5) risk1++;
            else if (result.riskLevel < 2.5) risk2++;
            else risk3++;

            maxHSI = Foam::max(maxHSI, result.heatStrainIndex);

            const Foam::scalar cellVolume = mesh.V()[cellI];
            weightedTa += result.taC*cellVolume;
            weightedTr += result.trC*cellVolume;
            weightedPa += result.paKPa*cellVolume;
            weightedEreq += result.ereqPeakWm2*cellVolume;
            weightedEmax += result.emaxMinWm2*cellVolume;
            weightedWreq += result.wreqPeak*cellVolume;
            weightedReqSweat += result.reqSweatRatePeak*cellVolume;
            weightedSweat += result.sweatRatePeak*cellVolume;
            weightedEact += result.eactPeakWm2*cellVolume;
            weightedTre += result.treFinalC*cellVolume;
            weightedTsk += result.tskFinalC*cellVolume;
            weightedWaterLoss += result.waterLossG*cellVolume;
            weightedDlimTre += result.dlimTreMin*cellVolume;
            weightedDlimLoss95 += result.dlimLoss95Min*cellVolume;
            weightedDlimLoss50 += result.dlimLoss50Min*cellVolume;
            weightedHSI += result.heatStrainIndex*cellVolume;
            totalVolume += cellVolume;
        }

        taOut.write();
        trOut.write();
        paOut.write();
        ereqOut.write();
        emaxOut.write();
        wreqOut.write();
        reqSweatRateOut.write();
        sweatRateOut.write();
        eactOut.write();
        treOut.write();
        tskOut.write();
        waterLossOut.write();
        dlimTreOut.write();
        dlimLoss95Out.write();
        dlimLoss50Out.write();
        heatStrainOut.write();
        riskOut.write();

        Foam::reduce(weightedTa, Foam::sumOp<Foam::scalar>());
        Foam::reduce(weightedTr, Foam::sumOp<Foam::scalar>());
        Foam::reduce(weightedPa, Foam::sumOp<Foam::scalar>());
        Foam::reduce(weightedEreq, Foam::sumOp<Foam::scalar>());
        Foam::reduce(weightedEmax, Foam::sumOp<Foam::scalar>());
        Foam::reduce(weightedWreq, Foam::sumOp<Foam::scalar>());
        Foam::reduce(weightedReqSweat, Foam::sumOp<Foam::scalar>());
        Foam::reduce(weightedSweat, Foam::sumOp<Foam::scalar>());
        Foam::reduce(weightedEact, Foam::sumOp<Foam::scalar>());
        Foam::reduce(weightedTre, Foam::sumOp<Foam::scalar>());
        Foam::reduce(weightedTsk, Foam::sumOp<Foam::scalar>());
        Foam::reduce(weightedWaterLoss, Foam::sumOp<Foam::scalar>());
        Foam::reduce(weightedDlimTre, Foam::sumOp<Foam::scalar>());
        Foam::reduce(weightedDlimLoss95, Foam::sumOp<Foam::scalar>());
        Foam::reduce(weightedDlimLoss50, Foam::sumOp<Foam::scalar>());
        Foam::reduce(weightedHSI, Foam::sumOp<Foam::scalar>());
        Foam::reduce(totalVolume, Foam::sumOp<Foam::scalar>());
        Foam::reduce(maxHSI, Foam::maxOp<Foam::scalar>());
        Foam::reduce(risk0, Foam::sumOp<Foam::label>());
        Foam::reduce(risk1, Foam::sumOp<Foam::label>());
        Foam::reduce(risk2, Foam::sumOp<Foam::label>());
        Foam::reduce(risk3, Foam::sumOp<Foam::label>());

        const Foam::scalar avgTa = weightedTa/Foam::max(SMALL, totalVolume);
        const Foam::scalar avgTr = weightedTr/Foam::max(SMALL, totalVolume);
        const Foam::scalar avgPa = weightedPa/Foam::max(SMALL, totalVolume);
        const Foam::scalar avgEreq = weightedEreq/Foam::max(SMALL, totalVolume);
        const Foam::scalar avgEmax = weightedEmax/Foam::max(SMALL, totalVolume);
        const Foam::scalar avgWreq = weightedWreq/Foam::max(SMALL, totalVolume);
        const Foam::scalar avgReqSweat = weightedReqSweat/Foam::max(SMALL, totalVolume);
        const Foam::scalar avgSweat = weightedSweat/Foam::max(SMALL, totalVolume);
        const Foam::scalar avgEact = weightedEact/Foam::max(SMALL, totalVolume);
        const Foam::scalar avgTre = weightedTre/Foam::max(SMALL, totalVolume);
        const Foam::scalar avgTsk = weightedTsk/Foam::max(SMALL, totalVolume);
        const Foam::scalar avgWaterLoss = weightedWaterLoss/Foam::max(SMALL, totalVolume);
        const Foam::scalar avgDlimTre = weightedDlimTre/Foam::max(SMALL, totalVolume);
        const Foam::scalar avgDlimLoss95 = weightedDlimLoss95/Foam::max(SMALL, totalVolume);
        const Foam::scalar avgDlimLoss50 = weightedDlimLoss50/Foam::max(SMALL, totalVolume);
        const Foam::scalar avgHSI = weightedHSI/Foam::max(SMALL, totalVolume);

        const Foam::label analyzedCells =
            Foam::returnReduce(cellsToAnalyze.size(), Foam::sumOp<Foam::label>());

        Foam::Info<< Foam::nl
            << "============= ISO 7933 DYNAMIC PHS ANALYSIS =============" << Foam::nl
            << "Analyzed cells:                  " << analyzedCells << Foam::nl
            << "Analysis volume:                 " << totalVolume << " m^3" << Foam::nl
            << "Humidity source:                 " << rhSource << Foam::nl
            << "Radiant temperature source:      " << trSource << Foam::nl
            << Foam::nl
            << "Metabolic rate M:               " << params.M << " W/m^2" << Foam::nl
            << "External work W:                " << params.W << " W/m^2" << Foam::nl
            << "Clothing insulation:            " << params.clo << " clo" << Foam::nl
            << "Permeability index im:          " << params.im << Foam::nl
            << "Exposure duration:              " << params.exposureDurationMin << " min" << Foam::nl
            << "Time step:                      " << params.timeStepMin << " min" << Foam::nl
            << Foam::nl
            << "Average Ta:                     " << avgTa << " degC" << Foam::nl
            << "Average Tr:                     " << avgTr << " degC" << Foam::nl
            << "Average Pa:                     " << avgPa << " kPa" << Foam::nl
            << "Average Ereq(peak):             " << avgEreq << " W/m^2" << Foam::nl
            << "Average Emax(min):              " << avgEmax << " W/m^2" << Foam::nl
            << "Average Wreq(peak):             " << avgWreq << Foam::nl
            << "Average req sweat(peak):        " << avgReqSweat << " g/(m^2 h)" << Foam::nl
            << "Average sweat(peak):            " << avgSweat << " g/(m^2 h)" << Foam::nl
            << "Average Eact(peak):             " << avgEact << " W/m^2" << Foam::nl
            << "Average Tre(final):             " << avgTre << " degC" << Foam::nl
            << "Average Tsk(final):             " << avgTsk << " degC" << Foam::nl
            << "Average water loss:             " << avgWaterLoss << " g" << Foam::nl
            << "Average DlimTre:                " << avgDlimTre << " min" << Foam::nl
            << "Average DlimLoss95:             " << avgDlimLoss95 << " min" << Foam::nl
            << "Average DlimLoss50:             " << avgDlimLoss50 << " min" << Foam::nl
            << "Average HeatStrainIndex:        " << avgHSI << Foam::nl
            << "Max HeatStrainIndex:            " << maxHSI << Foam::nl
            << Foam::nl
            << "Risk level counts (0/1/2/3):    "
            << risk0 << " / " << risk1 << " / " << risk2 << " / " << risk3 << Foam::nl
            << "=========================================================" << Foam::endl;
    }

    return 0;
}

// ************************************************************************* //
