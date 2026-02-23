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
    comfortFoam

Description
    This tool calculates thermal comfort parameters according to DIN EN ISO 7730:
    - PMV (Predicted Mean Vote), Range: -3 (cold) to 3 (hot)
    - PPD (Predicted Percentage of Dissatisfied), Range: 0 to 100%
    - DR  (Draught Rating), Range: 0 to 100%
    - TOp (Operative Temperature)
    - Mean radiant temperature
    - Turbulent intensity %

    Analysis can be limited to specific regions:
    - Use -cellSet <name> (or -setFields <name>) for a specific cellSet
    - Use -cellZone <name> for a specific cellZone
    - Without options: analyzes entire mesh

    For humidity calculations, run: buoyantHumiditySimpleFoam / buoyantHumidityPimpleFoam
    
    Supported humidity fields:
    - thermo:relHum (from buoyantHumiditySimpleFoam, range [0,1])
    - thermoRelHum (legacy naming, range [0,1])
    - relHum (legacy field, range [0,100] or [0,1])
    - Default constant value if no field available
    
    Supported radiation fields (in order of preference):
    - G: Incident radiation field (W/m^2) from radiation models
    - qr: Radiative heat flux field (W/m^2) from some radiation models
    - IDefault: Default radiation intensity (W/m^2/sr) from DOM models
    - Fallback: Area-weighted wall temperature calculation

    Usage examples:
    comfortFoam                           # Analyze entire mesh
    comfortFoam -cellSet occupancyZone    # Analyze only specified cellSet
    comfortFoam -setFields occupancyZone  # Alias for -cellSet
    comfortFoam -cellZone roomA           # Analyze only specified cellZone

Background
    DIN EN ISO 7730

Authors
    Thomas Tian
    Tobias Holzmann
    Manuel Scheu

Version
    3.3

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "wallFvPatch.H"
#include "externalWallHeatFluxTemperatureFvPatchScalarField.H"
#include "cellSet.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "radiationModel.H"
#include "atmBoundaryLayerInletVelocityFvPatchVectorField.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * Constants * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace ComfortConstants
    {
        // Physical constants
        const scalar stefanBoltzmannConstant = 5.67e-8;  // W/(m^2*K^4)
        const scalar baseMetabolicRate = 58.15;          // W/m^2
        
        // Heat loss coefficients (from ISO 7730)
        const scalar skinDiffusionCoeff = 3.05e-3;
        const scalar thermalCoeff = 6.99;
        const scalar basePressure = 5733;                // Pa
        const scalar respirationCoeff1 = 1.7e-5;
        const scalar respirationCoeff2 = 0.0014;
        const scalar respirationTemp = 34.0;             // degrees C
        const scalar respirationHumidity = 5867;         // Pa
        const scalar radiationCoeff = 3.96;
        const scalar thermalSensCoeff1 = 0.303;
        const scalar thermalSensCoeff2 = 0.036;  // ISO 7730 exponential coefficient
        const scalar thermalSensCoeff3 = 0.028;  // ISO 7730 additive constant
        
        // Clothing coefficients
        const scalar clothingFactor1 = 1.29;
        const scalar clothingFactor2 = 1.05;
        const scalar clothingFactor3 = 0.645;
        const scalar clothingThreshold = 0.078;
        
        // Convection coefficients
        const scalar forcedConvectionCoeff = 12.1;
        const scalar naturalConvectionCoeff = 2.38;
        const scalar naturalConvectionExp = 0.25;
        
        // Draft calculation coefficients
        const scalar draftBaseTemp = 34.0;               // degrees C
        const scalar draftVelThreshold = 0.05;           // m/s
        const scalar draftVelExp = 0.62;
        const scalar draftTurbCoeff = 0.37;
        const scalar draftConstant = 3.14;
        
        // Iteration parameters
        const label maxClothingIterations = 150;
        const scalar clothingConvergenceTol = 0.00015;  // ISO 7730 epsilon from BASIC line 340
        
        // PPD calculation coefficients
        const scalar ppdCoeff1 = 0.03353;
        const scalar ppdCoeff2 = 0.2179;
        const scalar ppdBase = 95.0;
        
        // Operative temperature correction factors
        const scalar opTempFactor1 = 0.5;  // v < 0.2 m/s
        const scalar opTempFactor2 = 0.6;  // 0.2 <= v <= 0.6 m/s
        const scalar opTempFactor3 = 0.7;  // v > 0.6 m/s
        const scalar velThreshold1 = 0.2;  // m/s
        const scalar velThreshold2 = 0.6;  // m/s
    }

// * * * * * * * * * * * * * * * Enums * * * * * * * * * * * * * * * * * * //

enum class ComfortCategory
{
    A,  // High comfort
    B,  // Medium comfort  
    C,  // Moderate comfort
    None // No category
};

enum class AnalysisRegionType
{
    EntireMesh,
    CellSet,
    CellZone
};

// Convert comfort category to string
Foam::word comfortCategoryToString(ComfortCategory category)
{
    switch (category)
    {
        case ComfortCategory::A: return "Category A (High comfort)";
        case ComfortCategory::B: return "Category B (Medium comfort)";
        case ComfortCategory::C: return "Category C (Moderate comfort)";
        default: return "No category (Poor comfort)";
    }
}

// Calculate volume-weighted averages for specific cells
void calculateVolumeWeightedAverages
(
    const fvMesh& mesh,
    const labelList& cellsToAnalyze,
    scalar& avgVelocity,
    scalar& avgTemperature,
    scalar& totalVolume
)
{
    const volVectorField& U = mesh.lookupObject<volVectorField>("U");
    const volScalarField& T = mesh.lookupObject<volScalarField>("T");
    
    vector weightedVelocity(Zero);
    scalar weightedTemperature(0);
    totalVolume = 0;
    
    forAll(cellsToAnalyze, i)
    {
        label cellI = cellsToAnalyze[i];
        scalar cellVolume = mesh.V()[cellI];
        
        weightedVelocity += U[cellI] * cellVolume;
        weightedTemperature += T[cellI] * cellVolume;
        totalVolume += cellVolume;
    }
    
    // Reduce for parallel operation
    reduce(weightedVelocity, sumOp<vector>());
    reduce(weightedTemperature, sumOp<scalar>());
    reduce(totalVolume, sumOp<scalar>());
    
    avgVelocity = mag(weightedVelocity / totalVolume);
    avgTemperature = weightedTemperature / totalVolume;
}
}

// * * * * * * * * * * * * * * * Functions * * * * * * * * * * * * * * * * //

// Validate input parameters
void validateInputParameters(scalar met, scalar clo, scalar wme, scalar rh)
{
    if (met < 0.8 || met > 4.0)
    {
        FatalErrorInFunction
            << "Metabolic rate out of valid range (0.8-4.0 met): " << met
            << abort(FatalError);
    }
    
    if (clo < 0 || clo > 2.0)
    {
        FatalErrorInFunction 
            << "Clothing insulation out of valid range (0-2.0 clo): " << clo
            << abort(FatalError);
    }
    
    if (wme < 0 || wme > met)
    {
        FatalErrorInFunction
            << "External work out of valid range (0-" << met << " met): " << wme
            << abort(FatalError);
    }
    
    if (rh < 0 || rh > 100)
    {
        FatalErrorInFunction
            << "Relative humidity out of valid range (0-100%): " << rh
            << abort(FatalError);
    }
}

// Calculate volume-weighted average velocity
Foam::vector calculateAverageVelocity(const fvMesh& mesh)
{
    const volVectorField& U = mesh.lookupObject<volVectorField>("U");
    vector avgU(Zero);
    scalar totalVolume(0);

    forAll(mesh.cells(), cellI)
    {
        avgU += U[cellI] * mesh.V()[cellI];
        totalVolume += mesh.V()[cellI];
    }

    return avgU / totalVolume;
}

// Calculate volume-weighted average temperature
Foam::scalar calculateAverageTemperature(const fvMesh& mesh)
{
    const volScalarField& T = mesh.lookupObject<volScalarField>("T");
    scalar avgT(0);
    scalar totalVolume(0);

    forAll(mesh.cells(), cellI)
    {
        avgT += T[cellI] * mesh.V()[cellI];
        totalVolume += mesh.V()[cellI];
    }

    return avgT / totalVolume;
}

// Calculate area-weighted radiation temperature from wall surfaces
Foam::scalar calculateRadiationTemperature
(
    const fvMesh& mesh,
    const fvPatchList& patches
)
{
    const volScalarField& T = mesh.lookupObject<volScalarField>("T");
    scalar weightedTemp(0);
    scalar totalArea(0);

    forAll(patches, patchI)
    {
        const label curPatch = patches[patchI].index();

        if (isType<wallFvPatch>(patches[patchI]))
        {
            scalar patchArea = gSum(mesh.magSf().boundaryField()[curPatch]);

            if (patchArea > SMALL)
            {
                weightedTemp += gSum
                (
                    mesh.magSf().boundaryField()[curPatch]
                  * T.boundaryField()[curPatch]
                );

                totalArea += patchArea;
            }
        }
    }

    return (weightedTemp / totalArea) - 273.0;  // Convert to Celsius (ISO 7730 uses 273.0)
}

// Calculate water vapour pressure
Foam::scalar calculateWaterVapourPressure(scalar temperature, scalar relativeHumidity)
{
    // Temperature should be in Kelvin, convert to Celsius for the formula
    // Note: Use 273.0 to match ISO 7730/C# implementation
    scalar tempCelsius = temperature - 273.0;
    return relativeHumidity * 10.0 
         * Foam::exp(16.6536 - (4030.183 / (tempCelsius + 235.0)));
}

// Calculate turbulent intensity
Foam::scalar calculateTurbulentIntensity
(
    scalar velocity,
    scalar k,  // turbulent kinetic energy
    scalar epsilon,  // turbulent dissipation rate
    scalar omega,    // specific dissipation rate
    const word& turbulenceModel
)
{
    if (velocity > SMALL)
    {
        if (turbulenceModel == "k-epsilon" || turbulenceModel == "k-omega" || turbulenceModel == "k-only")
        {
            // Standard formula: Tu = sqrt(2/3 * k) / U
            return Foam::sqrt(2.0/3.0 * k) / velocity;
        }
        else
        {
            // Default assumption for mixed convection
            return 0.4;  // 40%
        }
    }
    return 0.4;  // Default 40% for zero velocity
}

// Calculate draught rating (DR)
Foam::scalar calculateDraughtRating
(
    scalar temperature,
    scalar velocity,
    scalar turbulentIntensity
)
{
    using namespace ComfortConstants;
    
    if (velocity >= draftVelThreshold)
    {
        scalar dr = (draftBaseTemp - (temperature - 273.0))
                  * Foam::pow(velocity - draftVelThreshold, draftVelExp)
                  * ((draftTurbCoeff * velocity * turbulentIntensity) + draftConstant);
        
        return Foam::max(0.0, Foam::min(100.0, dr));
    }
    
    return 0.0;
}

// Solve for clothing surface temperature iteratively
// Returns Tuple2<scalar> with first = tcl (Celsius), second = xn (Kelvin/100)
Foam::Tuple2<scalar> calculateClothingSurfaceTemperature
(
    scalar airTemp,
    scalar velocity,
    scalar icl,
    scalar fcl,
    scalar radiationTemp,
    scalar metabolicRate,
    scalar externalWork
)
{
    using namespace ComfortConstants;
    
    // Convert airTemp from Kelvin to Celsius for the calculation
    // Note: Use 273.0 to match ISO 7730/C# implementation
    scalar airTempC = airTemp - 273.0;
    // ISO 7730 BASIC line 250: TCLA in Kelvin
    scalar tcla = airTemp + (35.5 - airTempC) / (3.5 * 6.45 * (icl + 0.1));
    // ISO 7730: iteration uses Kelvin divided by 100
    scalar xn = tcla / 100.0;
    scalar xf = xn;
    
    scalar p1 = icl * fcl;
    scalar p2 = p1 * radiationCoeff;
    scalar p3 = p1 * 100.0;
    scalar p4 = p1 * airTemp;  // Use Kelvin for p4
    scalar p5 = 308.7 - 0.028 * (metabolicRate - externalWork)
              + p2 * Foam::pow((radiationTemp + 273.0) / 100.0, 4);
    
    label iterationCount = 0;
    
    
    do
    {
        iterationCount++;
        xf = (xf + xn) / 2.0;  // ISO 7730 BASIC line 350
        
        scalar hcf = forcedConvectionCoeff * Foam::sqrt(velocity);
        scalar hcn = naturalConvectionCoeff * Foam::pow(mag(100.0 * xf - airTemp), naturalConvectionExp);
        scalar hc = Foam::max(hcf, hcn);
        
        xn = (p5 + p4 * hc - p2 * Foam::pow(xf, 4.0)) / (100.0 + p3 * hc);
        
        
        if (iterationCount > maxClothingIterations)
        {
            WarningInFunction
                << "Clothing temperature iteration did not converge after "
                << maxClothingIterations << " iterations" << endl;
            break;
        }
        
    } while (mag(xn - xf) > clothingConvergenceTol);
    
    // Return both tcl (Celsius) and xn (Kelvin/100)
    return Tuple2<scalar>(100.0 * xn - 273.0, xn);
}

// Calculate PMV (Predicted Mean Vote)
Foam::scalar calculatePMV
(
    scalar airTemp,
    scalar velocity,
    scalar relativeHumidity,
    scalar radiationTemp,
    scalar metabolicRate,
    scalar clothingInsulation,
    scalar externalWork
)
{
    using namespace ComfortConstants;
    
    // Convert airTemp from Kelvin to Celsius for consistency
    scalar airTempC = airTemp - 273.0;
    
    scalar icl = 0.155 * clothingInsulation;
    scalar fcl = (icl < clothingThreshold) ? 
                 (1.0 + clothingFactor1 * icl) : 
                 (clothingFactor2 + clothingFactor3 * icl);
    
    scalar pa = calculateWaterVapourPressure(airTemp, relativeHumidity);
    
    Tuple2<scalar> tclResult = calculateClothingSurfaceTemperature
    (
        airTemp, velocity, icl, fcl, radiationTemp, metabolicRate, externalWork
    );
    scalar tcl = tclResult.first();  // Temperature in Celsius
    scalar xn = tclResult.second();   // Kelvin/100
    
    // Calculate heat losses
    scalar hl1 = skinDiffusionCoeff * (basePressure - (thermalCoeff * (metabolicRate - externalWork)) - pa);
    
    scalar hl2 = 0.0;
    if ((metabolicRate - externalWork) > baseMetabolicRate)
    {
        hl2 = 0.42 * ((metabolicRate - externalWork) - baseMetabolicRate);
    }
    
    scalar hl3 = respirationCoeff1 * metabolicRate * (respirationHumidity - pa);
    scalar hl4 = respirationCoeff2 * metabolicRate * (respirationTemp - airTempC);
    
    // Use xn directly for HL5 calculation (like ISO 7730 BASIC line 480)
    scalar hl5 = radiationCoeff * fcl * 
                (Foam::pow(xn, 4) - Foam::pow((radiationTemp + 273.0) / 100.0, 4));
    
    scalar hcf = forcedConvectionCoeff * Foam::sqrt(velocity);
    scalar hcn = naturalConvectionCoeff * Foam::pow(mag(tcl - airTempC), naturalConvectionExp);
    scalar hc = Foam::max(hcf, hcn);
    
    scalar hl6 = fcl * hc * (tcl - airTempC);
    
    // ISO 7730 formula: TS = 0.303 * exp(-0.036 * M) + 0.028
    scalar ts = thermalSensCoeff1 * Foam::exp(-thermalSensCoeff2 * metabolicRate) + thermalSensCoeff3;
    
    return ts * ((metabolicRate - externalWork) - hl1 - hl2 - hl3 - hl4 - hl5 - hl6);
}

// Calculate PPD (Predicted Percentage of Dissatisfied)
Foam::scalar calculatePPD(scalar pmv)
{
    using namespace ComfortConstants;
    
    return 100.0 - ppdBase * Foam::exp(-ppdCoeff1 * Foam::pow(pmv, 4) - ppdCoeff2 * Foam::pow(pmv, 2));
}

// Calculate operative temperature
Foam::scalar calculateOperativeTemperature
(
    scalar airTemp,
    scalar radiationTemp,
    scalar velocity
)
{
    using namespace ComfortConstants;
    
    scalar correctionFactor;
    
    if (velocity < velThreshold1)
    {
        correctionFactor = opTempFactor1;
    }
    else if (velocity <= velThreshold2)  // Fixed logical error: AND instead of OR
    {
        correctionFactor = opTempFactor2;
    }
    else
    {
        correctionFactor = opTempFactor3;
    }
    
    return correctionFactor * airTemp + (1.0 - correctionFactor) * (radiationTemp + 273.0);
}

// Analyze comfort category according to ISO 7730
ComfortCategory analyzeComfortCategory(scalar pmv, scalar dr, scalar ppd)
{
    if (mag(pmv) <= 0.2 && dr < 10.0 && ppd < 6.0)
    {
        return ComfortCategory::A;
    }
    else if (mag(pmv) <= 0.5 && dr < 20.0 && ppd < 10.0)
    {
        return ComfortCategory::B;
    }
    else if (mag(pmv) <= 0.7 && dr < 30.0 && ppd < 15.0)
    {
        return ComfortCategory::C;
    }
    else
    {
        return ComfortCategory::None;
    }
}

// Normalize selection labels to valid local cell labels.
// In parallel, this also handles the case where a set/zone stores global labels.
Foam::labelList normalizeSelectedCells
(
    const fvMesh& mesh,
    const labelUList& rawCells,
    const word& selectionKind,
    const word& selectionName
)
{
    const label nLocalCells = mesh.nCells();
    bool hasOutOfLocalRange = false;

    forAll(rawCells, i)
    {
        const label cellI = rawCells[i];
        if (cellI < 0 || cellI >= nLocalCells)
        {
            hasOutOfLocalRange = true;
            break;
        }
    }

    if (Pstream::parRun())
    {
        reduce(hasOutOfLocalRange, orOp<bool>());
    }

    DynamicList<label> validCells(rawCells.size());
    label nDiscarded = 0;

    if (Pstream::parRun() && hasOutOfLocalRange)
    {
        const globalIndex& globalCellAddr = mesh.globalData().globalMeshCellAddr();

        forAll(rawCells, i)
        {
            const label globalCellI = rawCells[i];
            if (globalCellAddr.isLocal(globalCellI))
            {
                validCells.append(globalCellAddr.toLocal(globalCellI));
            }
            else
            {
                ++nDiscarded;
            }
        }

        const label totalDiscarded = returnReduce(nDiscarded, sumOp<label>());
        if (totalDiscarded > 0)
        {
            Info<< "Converted global labels for " << selectionKind << " "
                << selectionName << " and ignored " << totalDiscarded
                << " non-local entries across all processors." << endl;
        }
    }
    else
    {
        forAll(rawCells, i)
        {
            const label localCellI = rawCells[i];
            if (localCellI >= 0 && localCellI < nLocalCells)
            {
                validCells.append(localCellI);
            }
            else
            {
                ++nDiscarded;
            }
        }

        const label totalDiscarded = returnReduce(nDiscarded, sumOp<label>());
        if (totalDiscarded > 0)
        {
            WarningInFunction
                << "Ignoring " << totalDiscarded << " invalid labels in "
                << selectionKind << " " << selectionName << endl;
        }
    }

    labelList cellsToAnalyze(validCells.size());
    forAll(validCells, i)
    {
        cellsToAnalyze[i] = validCells[i];
    }

    // Defensive duplicate removal to avoid double-counting volume.
    if (!cellsToAnalyze.empty())
    {
        sort(cellsToAnalyze);
        label writeI = 0;
        for (label readI = 1; readI < cellsToAnalyze.size(); ++readI)
        {
            if (cellsToAnalyze[readI] != cellsToAnalyze[writeI])
            {
                cellsToAnalyze[++writeI] = cellsToAnalyze[readI];
            }
        }
        cellsToAnalyze.setSize(writeI + 1);
    }

    return cellsToAnalyze;
}

// Get list of cells to analyze (from CLI/dictionary cellSet/cellZone or entire mesh)
Foam::labelList getCellsToAnalyze
(
    const fvMesh& mesh,
    const argList& args,
    const dictionary& comfortFoamDict,
    AnalysisRegionType& selectionType,
    word& selectionName
)
{
    labelList cellsToAnalyze;
    selectionType = AnalysisRegionType::EntireMesh;
    selectionName.clear();

    const bool cliCellSetFound = args.found("cellSet");
    const bool cliSetFieldsFound = args.found("setFields");
    const bool cliCellZoneFound = args.found("cellZone");

    word cliCellSetName;
    if (cliCellSetFound)
    {
        cliCellSetName = args.get<word>("cellSet");
    }
    if (cliSetFieldsFound)
    {
        const word setFieldsName = args.get<word>("setFields");
        if (cliCellSetFound && setFieldsName != cliCellSetName)
        {
            FatalErrorInFunction
                << "Both -cellSet and -setFields were specified with different names: "
                << cliCellSetName << " vs " << setFieldsName << nl
                << "Use only one option (or use the same name)."
                << abort(FatalError);
        }
        cliCellSetName = setFieldsName;
    }

    const bool dictCellSetFound = comfortFoamDict.found("cellSet");
    const bool dictSetFieldsFound = comfortFoamDict.found("setFields");
    const bool dictCellZoneFound = comfortFoamDict.found("cellZone");

    word dictCellSetName;
    if (dictCellSetFound)
    {
        dictCellSetName = word(comfortFoamDict.lookup("cellSet"));
    }
    if (dictSetFieldsFound)
    {
        const word setFieldsName(word(comfortFoamDict.lookup("setFields")));
        if (dictCellSetFound && setFieldsName != dictCellSetName)
        {
            FatalErrorInFunction
                << "Both dictionary entries cellSet and setFields are defined with "
                << "different names: " << dictCellSetName << " vs " << setFieldsName << nl
                << "Use only one of these entries (or use the same name)."
                << abort(FatalError);
        }
        dictCellSetName = setFieldsName;
    }

    const bool cliAnySelection =
        cliCellZoneFound || !cliCellSetName.empty();
    const bool dictAnySelection =
        dictCellZoneFound || !dictCellSetName.empty();

    if (cliAnySelection && dictAnySelection)
    {
        Info<< "Both CLI and comfortFoamDict region selections are specified. "
            << "Using CLI selection." << endl;
    }

    // CLI selection has priority over dictionary selection
    bool useCellSet = false;
    bool useCellZone = false;
    word cellSetName;
    word cellZoneName;

    if (cliAnySelection)
    {
        useCellSet = !cliCellSetName.empty();
        useCellZone = cliCellZoneFound;
        cellSetName = cliCellSetName;
        if (cliCellZoneFound)
        {
            cellZoneName = args.get<word>("cellZone");
        }
    }
    else
    {
        useCellSet = !dictCellSetName.empty();
        useCellZone = dictCellZoneFound;
        cellSetName = dictCellSetName;
        if (dictCellZoneFound)
        {
            cellZoneName = word(comfortFoamDict.lookup("cellZone"));
        }
    }

    if (useCellSet && useCellZone)
    {
        FatalErrorInFunction
            << "Both cellSet and cellZone region selections are specified. "
            << "Use only one."
            << abort(FatalError);
    }

    // Check for cellSet option
    if (useCellSet)
    {
        Info<< "Loading cellSet: " << cellSetName << endl;
        
        cellSet selectedCells(mesh, cellSetName);

        const labelList rawCells = selectedCells.toc();
        cellsToAnalyze =
            normalizeSelectedCells(mesh, rawCells, "cellSet", cellSetName);

        const label nCells = returnReduce(cellsToAnalyze.size(), sumOp<label>());
        if (nCells == 0)
        {
            FatalErrorInFunction
                << "cellSet " << cellSetName
                << " is empty, invalid, or not available on this case."
                << abort(FatalError);
        }

        Info<< "Analyzing " << nCells 
            << " cells from cellSet " << cellSetName << endl;
        selectionType = AnalysisRegionType::CellSet;
        selectionName = cellSetName;
    }
    // Check for cellZone option
    else if (useCellZone)
    {
        Info<< "Loading cellZone: " << cellZoneName << endl;
        
        const cellZoneMesh& cellZones = mesh.cellZones();
        label zoneID = cellZones.findZoneID(cellZoneName);
        
        if (zoneID == -1)
        {
            FatalErrorInFunction
                << "cellZone " << cellZoneName << " does not exist" << nl
                << "Available cellZones: " << cellZones.names()
                << abort(FatalError);
        }
        
        const cellZone& cz = cellZones[zoneID];
        const labelList rawCells = cz;
        cellsToAnalyze =
            normalizeSelectedCells(mesh, rawCells, "cellZone", cellZoneName);

        const label nCells = returnReduce(cellsToAnalyze.size(), sumOp<label>());
        if (nCells == 0)
        {
            FatalErrorInFunction
                << "cellZone " << cellZoneName
                << " is empty after parallel label mapping."
                << abort(FatalError);
        }

        Info<< "Analyzing " << nCells 
            << " cells from cellZone " << cellZoneName << endl;
        selectionType = AnalysisRegionType::CellZone;
        selectionName = cellZoneName;
    }
    // Use entire mesh
    else
    {
        cellsToAnalyze.setSize(mesh.nCells());
        forAll(cellsToAnalyze, i)
        {
            cellsToAnalyze[i] = i;
        }
        
        label nCells = returnReduce(cellsToAnalyze.size(), sumOp<label>());
        Info<< "Analyzing entire mesh: " << nCells << " cells" << endl;
        selectionType = AnalysisRegionType::EntireMesh;
        selectionName.clear();
    }
    
    return cellsToAnalyze;
}

// * * * * * * * * * * * * * * * * Main Program * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::timeSelector::addOptions();
    
    argList::addBoolOption
    (
        "validate",
        "Run validation test with ISO 7730 standard conditions"
    );
    
    argList::addOption
    (
        "validateAirTemp",
        "value",
        "Air temperature in Celsius for validation (default: 22)"
    );
    
    argList::addOption
    (
        "validateRadTemp", 
        "value",
        "Radiation temperature in Celsius for validation (default: 22)"
    );
    
    argList::addOption
    (
        "validateVelocity",
        "value", 
        "Air velocity in m/s for validation (default: 0.1)"
    );
    
    argList::addOption
    (
        "validateRH",
        "value",
        "Relative humidity in % for validation (default: 60)"
    );
    
    argList::addOption
    (
        "validateMet",
        "value",
        "Metabolic rate in met for validation (default: 1.2)"
    );
    
    argList::addOption
    (
        "validateClo",
        "value",
        "Clothing insulation in clo for validation (default: 0.5)"
    );
    
    argList::addOption
    (
        "cellSet",
        "name",
        "Analyze only cells in the specified cellSet instead of entire mesh"
    );

    argList::addOption
    (
        "setFields",
        "name",
        "Alias for -cellSet (analyze only cells in the specified cellSet)"
    );
    
    argList::addOption
    (
        "cellZone", 
        "name",
        "Analyze only cells in the specified cellZone instead of entire mesh"
    );
    
    #include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"
    
    // Check if validation mode is requested
    if (args.found("validate"))
    {
        Info<< "\n========== VALIDATION MODE ==========\n" << endl;
        
        // Get validation parameters
        scalar valAirTemp = args.getOrDefault<scalar>("validateAirTemp", 22.0);
        scalar valRadTemp = args.getOrDefault<scalar>("validateRadTemp", 22.0);
        scalar valVelocity = args.getOrDefault<scalar>("validateVelocity", 0.1);
        scalar valRH = args.getOrDefault<scalar>("validateRH", 60.0);
        scalar valMet = args.getOrDefault<scalar>("validateMet", 1.2);
        scalar valClo = args.getOrDefault<scalar>("validateClo", 0.5);
        scalar valWme = 0.0; // External work is typically 0
        
        Info<< "Validation parameters:" << nl
            << "  Air temperature:        " << valAirTemp << " degrees C" << nl
            << "  Radiation temperature:  " << valRadTemp << " degrees C" << nl
            << "  Air velocity:           " << valVelocity << " m/s" << nl
            << "  Relative humidity:      " << valRH << " %" << nl
            << "  Metabolic rate:         " << valMet << " met" << nl
            << "  Clothing insulation:    " << valClo << " clo" << nl
            << "  External work:          " << valWme << " met" << nl << endl;
        
        // Convert to SI units
        // Note: ISO 7730 uses 273.0 for K conversion
        scalar airTempK = valAirTemp + 273.0;  // Match C# code exactly
        scalar metRate = valMet * 58.15;
        scalar wmeRate = valWme * 58.15;
        
        // Calculate clothing parameters
        scalar icl = 0.155 * valClo;
        scalar fcl = (icl < 0.078) ? 
                     (1.0 + 1.29 * icl) : 
                     (1.05 + 0.645 * icl);
        
        Info<< "Intermediate values:" << nl
            << "  Icl:                    " << icl << " m^2*K/W" << nl
            << "  fcl:                    " << fcl << nl
            << "  M:                      " << metRate << " W/m^2" << nl
            << "  W:                      " << wmeRate << " W/m^2" << nl;
        
        // Calculate clothing surface temperature
        Tuple2<scalar> tclResult = calculateClothingSurfaceTemperature
        (
            airTempK, valVelocity, icl, fcl, valRadTemp, metRate, wmeRate
        );
        scalar tcl = tclResult.first();
        scalar xn = tclResult.second();
        
        // Recalculate initial values for debugging
        scalar tcla_debug = airTempK + (35.5 - valAirTemp) / (3.5 * 6.45 * (icl + 0.1));
        scalar xn_initial = tcla_debug / 100.0;
        
        Info<< "  TCLA initial:           " << tcla_debug << " K" << nl
            << "  XN initial:             " << xn_initial << nl
            << "  Clothing temp (tcl):    " << tcl << " degrees C" << nl
            << "  Iteration result (xn):  " << xn << nl;
            
        // Calculate heat losses for debugging
        scalar pa = calculateWaterVapourPressure(airTempK, valRH);
        Info<< "  Water vapor pressure:   " << pa << " Pa" << nl;
        
        scalar hl1 = 3.05e-3 * (5733 - 6.99 * (metRate - wmeRate) - pa);
        scalar hl2 = (metRate - wmeRate > 58.15) ? 0.42 * ((metRate - wmeRate) - 58.15) : 0.0;
        scalar hl3 = 1.7e-5 * metRate * (5867 - pa);
        scalar hl4 = 0.0014 * metRate * (34.0 - valAirTemp);
        scalar hl5 = 3.96 * fcl * (Foam::pow(xn, 4) - Foam::pow((valRadTemp + 273.0) / 100.0, 4));
        
        scalar hcf = 12.1 * Foam::sqrt(valVelocity);
        scalar hcn = 2.38 * Foam::pow(mag(tcl - valAirTemp), 0.25);
        scalar hc = Foam::max(hcf, hcn);
        scalar hl6 = fcl * hc * (tcl - valAirTemp);
        
        Info<< nl << "Heat losses:" << nl
            << "  HL1 (skin diffusion):   " << hl1 << " W/m^2" << nl
            << "  HL2 (sweat):            " << hl2 << " W/m^2" << nl
            << "  HL3 (resp. latent):     " << hl3 << " W/m^2" << nl
            << "  HL4 (resp. sensible):   " << hl4 << " W/m^2" << nl
            << "  HL5 (radiation):        " << hl5 << " W/m^2" << nl
            << "  HL6 (convection):       " << hl6 << " W/m^2" << nl
            << "  Total heat loss:        " << hl1+hl2+hl3+hl4+hl5+hl6 << " W/m^2" << nl;
            
        scalar ts = 0.303 * Foam::exp(-0.036 * metRate) + 0.028;
        scalar balance = (metRate - wmeRate) - hl1 - hl2 - hl3 - hl4 - hl5 - hl6;
        Info<< "  TS coefficient:         " << ts << nl
            << "  Heat balance:           " << balance << " W/m^2" << nl;
        
        // Calculate PMV directly from the already computed values
        scalar pmv = ts * balance;
        
        // Calculate PPD
        scalar ppd = calculatePPD(pmv);
        
        // Calculate DR
        scalar dr = calculateDraughtRating
        (
            valAirTemp + 273.0,  // Convert to Kelvin
            valVelocity,
            0.4  // Assume 40% turbulence intensity for validation (as fraction)
        );
        
        Info<< nl << "VALIDATION RESULTS:" << nl
            << "  PMV:                    " << pmv << nl
            << "  PPD:                    " << ppd << " %" << nl
            << "  DR:                     " << dr << " %" << nl;
        
        Info<< nl << "Expected values (ISO 7730):" << nl
            << "  PMV:                    -0.75" << nl
            << "  PPD:                    17 %" << nl;
        
        Info<< nl << "Deviation:" << nl
            << "  PMV error:              " << mag(pmv - (-0.75)) << " (" 
            << mag(pmv - (-0.75))/0.75*100 << "%)" << nl
            << "  PPD error:              " << mag(ppd - 17.0) << " %" << nl;
            
        Info<< "\n====================================\n" << endl;
        
        return 0;
    }

    // Get times list
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        Info<< "Time = " << runTime.timeName() << endl;

        #include "createFields.H"

        // Load temperature field
        volScalarField T(THeader, mesh);
        const fvPatchList& patches = T.mesh().boundary();

        // Check for humidity field - multiple possible sources
        bool humidityAvailable = false;
        autoPtr<volScalarField> humidityField;
        word humiditySource("ISO7730Dict RH");

        const wordList humidityCandidates
        ({
            "thermo:relHum",
            "thermoRelHum",
            "relHum",
            "RH",
            "relativeHumidity"
        });

        forAll(humidityCandidates, i)
        {
            IOobject humidityHeader
            (
                humidityCandidates[i],
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT
            );

            if (humidityHeader.typeHeaderOk<volScalarField>())
            {
                humidityField.reset(new volScalarField(humidityHeader, mesh));
                humidityAvailable = true;
                humiditySource = humidityCandidates[i];
                break;
            }
        }

        if (humidityAvailable)
        {
            Info<< "Using humidity field: " << humiditySource << endl;
        }
        else
        {
            Info<< "No humidity field found - using default value: " << RH1 << "%" << endl;
        }

        // Check for turbulence field
        // Check for turbulence fields
        bool turbulenceAvailable = false;
        autoPtr<volScalarField> kField;
        autoPtr<volScalarField> epsilonField;
        autoPtr<volScalarField> omegaField;
        
        word turbulenceModel = "none";

        if (kHeader.typeHeaderOk<volScalarField>())
        {
            kField.reset(new volScalarField(kHeader, mesh));
            turbulenceAvailable = true;
            
            // Check which turbulence model is used
            if (epsilonHeader.typeHeaderOk<volScalarField>())
            {
                epsilonField.reset(new volScalarField(epsilonHeader, mesh));
                turbulenceModel = "k-epsilon";
                Info<< "Turbulence model: k-epsilon" << endl;
            }
            else if (omegaHeader.typeHeaderOk<volScalarField>())
            {
                omegaField.reset(new volScalarField(omegaHeader, mesh));
                turbulenceModel = "k-omega";
                Info<< "Turbulence model: k-omega SST" << endl;
            }
            else
            {
                turbulenceModel = "k-only";
                Info<< "Turbulence field k available (no epsilon/omega found)" << endl;
            }
        }
        else
        {
            Info<< "No turbulence fields available - assuming Tu = 40%" << endl;
        }

        // Load velocity field
        volVectorField U(UHeader, mesh);
        
        // Get cells to analyze (cellSet/cellZone from CLI or comfortFoamDict)
        AnalysisRegionType selectionType;
        word selectionName;
        labelList cellsToAnalyze =
            getCellsToAnalyze(mesh, args, comfortFoamDict, selectionType, selectionName);
        
        // Calculate averages for analysis region
        scalar avgVelocityMag, avgTemperature, totalAnalysisVolume;
        calculateVolumeWeightedAverages(mesh, cellsToAnalyze, avgVelocityMag, avgTemperature, totalAnalysisVolume);

        // Validate input parameters
        validateInputParameters(met, clo, wme, RH1);

        // Pre-calculate constants
        const scalar metRate = met * ComfortConstants::baseMetabolicRate;
        const scalar wmeRate = wme * ComfortConstants::baseMetabolicRate;

        // Volume-weighted averages for output
        scalar volumeWeightedPMV(0);
        scalar volumeWeightedPPD(0);
        scalar volumeWeightedDR(0);
        scalar volumeWeightedTOp(0);
        scalar volumeWeightedRadTemp(0);
        scalar volumeWeightedRH(0);
        scalar volumeWeightedTu(0);
        scalar totalVolume(0);

        // Radiation source selection must be rank-consistent in parallel.
        const bool useGField = G.headerOk();
        const bool useQrField = (!useGField && qrHeader.typeHeaderOk<volScalarField>());
        const bool useIDefaultField =
            (!useGField && !useQrField && IDefaultHeader.typeHeaderOk<volScalarField>());

        autoPtr<volScalarField> qrField;
        autoPtr<volScalarField> IDefaultField;
        scalar fallbackRadTemp = 0.0;

        if (useQrField)
        {
            qrField.reset(new volScalarField(qrHeader, mesh));
        }
        else if (useIDefaultField)
        {
            IDefaultField.reset(new volScalarField(IDefaultHeader, mesh));
        }
        else
        {
            // Uses gSum internally, so evaluate once per time-step (not per cell).
            fallbackRadTemp = calculateRadiationTemperature(mesh, patches);
        }

        // Main calculation loop - only over selected cells
        forAll(cellsToAnalyze, i)
        {
            label cellI = cellsToAnalyze[i];
            // Clamp temperature to reasonable range
            scalar cellTemp = Foam::min(400.0, T[cellI]);
            
            // Determine radiation temperature from available radiation fields
            scalar cellRadTemp;
            if (useGField)
            {
                // Use incident radiation field G (primary choice)
                scalar gValue = Foam::max(0.0, Foam::min(50000.0, G[cellI]));
                cellRadTemp = Foam::pow(gValue / (4.0 * ComfortConstants::stefanBoltzmannConstant), 0.25) - 273.0;
            }
            else if (qrField.valid())
            {
                // Use radiative heat flux field qr
                // qr is typically the net radiative heat flux in W/m^2
                // For a gray surface: qr = epsilon * sigma * (T^4 - T_rad^4)
                // Assuming epsilon = 0.9 and using local temperature to estimate T_rad
                scalar epsilon = 0.9;
                scalar localTemp = cellTemp;
                scalar qrValue = qrField()[cellI];
                // Solve for T_rad from: qr = epsilon * sigma * (T^4 - T_rad^4)
                scalar T4_rad = Foam::pow(localTemp, 4) - qrValue / (epsilon * ComfortConstants::stefanBoltzmannConstant);
                if (T4_rad > 0)
                {
                    cellRadTemp = Foam::pow(T4_rad, 0.25) - 273.0;
                }
                else
                {
                    // Fallback to local temperature if calculation gives invalid result
                    cellRadTemp = localTemp - 273.0;
                }
            }
            else if (IDefaultField.valid())
            {
                // Use default radiation intensity field IDefault
                // IDefault is radiation intensity in W/m^2/sr
                // Total irradiance G = 4*pi*I for isotropic radiation
                scalar gValue = 4.0 * constant::mathematical::pi * IDefaultField()[cellI];
                gValue = Foam::max(0.0, Foam::min(50000.0, gValue));
                cellRadTemp = Foam::pow(gValue / (4.0 * ComfortConstants::stefanBoltzmannConstant), 0.25) - 273.0;
            }
            else
            {
                // Fallback: Use area-weighted wall temperature
                cellRadTemp = fallbackRadTemp;
            }
            
            // Calculate relative humidity
            scalar cellRH;
            if (humidityAvailable)
            {
                // For thermo::relHum field, values are typically already in [0,1] range
                // For relHum field, values might be in [0,100] range
                scalar rawHumidity = humidityField()[cellI];
                
                // Auto-detect range and convert to percentage
                if (rawHumidity <= 1.0)
                {
                    cellRH = rawHumidity * 100.0;  // Convert from fraction to percentage
                }
                else
                {
                    cellRH = rawHumidity;  // Already in percentage
                }
                
                // Clamp to valid range
                cellRH = Foam::max(0.0, Foam::min(100.0, cellRH));
            }
            else
            {
                cellRH = RH1;  // Use default value
            }
            
            // Calculate turbulent intensity
            scalar Tu = calculateTurbulentIntensity
            (
                mag(U[cellI]),
                turbulenceAvailable && kField.valid() ? kField()[cellI] : 0.0,
                turbulenceAvailable && epsilonField.valid() ? epsilonField()[cellI] : 0.0,
                turbulenceAvailable && omegaField.valid() ? omegaField()[cellI] : 0.0,
                turbulenceModel
            );
            
            // Calculate comfort parameters
            scalar pmv = calculatePMV
            (
                cellTemp,
                mag(U[cellI]),
                cellRH,
                cellRadTemp,
                metRate,
                clo,
                wmeRate
            );
            
            scalar ppd = calculatePPD(pmv);
            
            scalar dr = calculateDraughtRating(cellTemp, mag(U[cellI]), Tu);
            
            scalar tOp = calculateOperativeTemperature(cellTemp, cellRadTemp, mag(U[cellI]));
            
            // Store results in fields
            PMV[cellI] = pmv;
            PPD[cellI] = ppd;
            DR[cellI] = dr;
            TOp[cellI] = tOp;
            
            // Accumulate volume-weighted averages
            scalar cellVolume = mesh.V()[cellI];
            volumeWeightedPMV += pmv * cellVolume;
            volumeWeightedPPD += ppd * cellVolume;
            volumeWeightedDR += dr * cellVolume;
            volumeWeightedTOp += tOp * cellVolume;
            volumeWeightedRadTemp += cellRadTemp * cellVolume;
            volumeWeightedRH += cellRH * cellVolume;
            volumeWeightedTu += Tu * cellVolume;
            totalVolume += cellVolume;
        } // End of forAll(cellsToAnalyze, i)

        // Write fields
        DR.write();
        PMV.write();
        PPD.write();
        TOp.write();

        // Reduce parallel values
        reduce(volumeWeightedPMV, sumOp<scalar>());
        reduce(volumeWeightedPPD, sumOp<scalar>());
        reduce(volumeWeightedDR, sumOp<scalar>());
        reduce(volumeWeightedTOp, sumOp<scalar>());
        reduce(volumeWeightedRadTemp, sumOp<scalar>());
        reduce(volumeWeightedRH, sumOp<scalar>());
        reduce(volumeWeightedTu, sumOp<scalar>());
        reduce(totalVolume, sumOp<scalar>());
        
        // Calculate final averages
        scalar avgPMV = volumeWeightedPMV / totalVolume;
        scalar avgPPD = volumeWeightedPPD / totalVolume;
        scalar avgDR = volumeWeightedDR / totalVolume;
        scalar avgTOp = volumeWeightedTOp / totalVolume;
        scalar avgRadTemp = volumeWeightedRadTemp / totalVolume;
        scalar avgRH = volumeWeightedRH / totalVolume;
        scalar avgTu = volumeWeightedTu / totalVolume;

        // Analyze comfort category
        ComfortCategory category = analyzeComfortCategory(avgPMV, avgDR, avgPPD);

        // Output results
        Info<< nl << "============ THERMAL COMFORT ANALYSIS RESULTS ============" << nl;
        
        // Show analysis region info
        if (selectionType == AnalysisRegionType::CellSet)
        {
            Info<< "Analysis region: cellSet " << selectionName << nl;
        }
        else if (selectionType == AnalysisRegionType::CellZone)
        {
            Info<< "Analysis region: cellZone " << selectionName << nl;
        }
        else
        {
            Info<< "Analysis region: entire mesh" << nl;
        }
        
        label totalCells = returnReduce(cellsToAnalyze.size(), sumOp<label>());
        Info<< "Analyzed cells: " << totalCells << nl
            << "Analysis volume: " << totalVolume << " m^3" << nl
            << nl
            << "Mean radiation temperature:     " << avgRadTemp << " degrees C" << nl
            << "Average air temperature:        " << avgTemperature - 273.0 << " degrees C" << nl
            << "Average air velocity:           " << avgVelocityMag << " m/s" << nl
            << "Average relative humidity:      " << avgRH << " %" << nl
            << "Average operative temperature:  " << avgTOp - 273.0 << " degrees C" << nl
            << "Average turbulent intensity:    " << avgTu * 100.0 << " %" << nl
            << nl
            << "COMFORT INDICES:" << nl
            << "PMV (Predicted Mean Vote):      " << avgPMV << nl
            << "PPD (Predicted % Dissatisfied): " << avgPPD << " %" << nl
            << "DR (Draught Rating):            " << avgDR << " %" << nl
            << nl
            << "COMFORT CATEGORY: " << comfortCategoryToString(category) << nl
            << "=========================================================" << endl;
    
    } // End of forAll(timeDirs, timei)

    return 0;
}

// ************************************************************************* //
