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
    ASHRAE55

Description
    This tool calculates thermal comfort according to ASHRAE 55-2020:
    - Operative room air temperature (TOp)
    - ASHRAE 55 adaptive comfort compliance (80% and 90% acceptability)
    - Designed for OpenFOAM 2412+ with solar radiation models
    - Fallback support for older versions using EPW solar calculations

Background
    ANSI/ASHRAE Standard 55-2020 - Adaptive Comfort Model
    Optimized for OpenFOAM 2412+ solar radiation capabilities

Author
    Thomas Tian

Version
    1.2 (Enhanced with improved EPW parser and detailed statistics)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "wallFvPatch.H"
#include "externalWallHeatFluxTemperatureFvPatchScalarField.H"
#include "cellSet.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "radiationModel.H"
#include <fstream>
#include <sstream>
#include <vector>

// Structure to hold EPW weather data
struct EPWData
{
    Foam::List<Foam::scalar> temperature;     // Hourly temperature data
    Foam::List<Foam::scalar> globalRadiation; // Global horizontal radiation (Wh/m2)
    Foam::List<Foam::scalar> directRadiation; // Direct normal radiation (Wh/m2)
    Foam::List<Foam::scalar> diffuseRadiation; // Diffuse horizontal radiation (Wh/m2)
    Foam::scalar latitude;                    // Site latitude from EPW header
    Foam::scalar longitude;                   // Site longitude from EPW header
    Foam::scalar timeZone;                    // Time zone (GMT offset) from EPW header
    Foam::string locationName;                // Location name from EPW header
};

namespace Foam
{
namespace ASHRAE55Constants
{
    const scalar baseMetabolicRate = 58.15;
    const scalar skinDiffusionCoeff = 3.05e-3;
    const scalar thermalCoeff = 6.99;
    const scalar basePressure = 5733;
    const scalar respirationCoeff1 = 1.7e-5;
    const scalar respirationCoeff2 = 0.0014;
    const scalar respirationTemp = 34.0;
    const scalar respirationHumidity = 5867;
    const scalar radiationCoeff = 3.96;
    const scalar thermalSensCoeff1 = 0.303;
    const scalar thermalSensCoeff2 = 0.036;
    const scalar thermalSensCoeff3 = 0.028;
    const scalar clothingFactor1 = 1.29;
    const scalar clothingFactor2 = 1.05;
    const scalar clothingFactor3 = 0.645;
    const scalar clothingThreshold = 0.078;
    const scalar forcedConvectionCoeff = 12.1;
    const scalar naturalConvectionCoeff = 2.38;
    const scalar naturalConvectionExp = 0.25;
    const label maxClothingIterations = 150;
    const scalar clothingConvergenceTol = 0.00015;
    const scalar ppdCoeff1 = 0.03353;
    const scalar ppdCoeff2 = 0.2179;
    const scalar ppdBase = 95.0;
}
}

// * * * * * * * * * * * * * * * * * * * Functions  * * * * * * * * * * * * * //

// If no Radiation flux is found, calculate the surface temperature as average Ts
Foam::scalar radiationTemperature
(
    const Foam::fvMesh& mesh_,
    const Foam::fvPatchList& Patches_
)
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    scalar PatchSumTemp(0), area(0), sumArea(0);

    forAll (Patches_, patchI)
    {
        const label curPatch = Patches_[patchI].index();

        if (isType<wallFvPatch>( Patches_[patchI] ))
        {
            area = gSum(mesh_.magSf().boundaryField()[curPatch]);

            if (area > 0)
            {
                PatchSumTemp +=
                    gSum
                    (
                        mesh_.magSf().boundaryField()[curPatch]
                      * T.boundaryField()[curPatch]
                    );

                sumArea += area;
            }
        }
    }

    return (PatchSumTemp / sumArea) - 273.15;
}

// Enhanced solar calculator - ONLY use if no OpenFOAM solar model is active
// Calculate solar MRT using EPW radiation data
Foam::scalar calculateSolarMRT
(
    const EPWData& epwData,
    const Foam::scalar& dayOfYear,
    const Foam::scalar& hour,
    const Foam::scalar& latitude,
    const Foam::scalar& longitude,
    const Foam::scalar& airTemp,
    const bool& showDebug = false
)
{
    // Solar position calculation with proper longitude correction
    scalar dayAngle = 2.0 * M_PI * (dayOfYear - 1) / 365.0;
    scalar declination = 23.45 * Foam::sin(dayAngle + 2.0 * M_PI * (284.0/365.0)) * M_PI/180.0;
    
    // Equation of time correction (in minutes) - standard formula
    scalar eqTimeMinutes = 229.183 * (0.000075 + 0.001868 * Foam::cos(dayAngle) - 
                          0.032077 * Foam::sin(dayAngle) - 0.014615 * Foam::cos(2*dayAngle) - 
                          0.040849 * Foam::sin(2*dayAngle));
    
    // Use time zone from EPW data if available, otherwise estimate from longitude
    scalar timeZoneOffset = epwData.timeZone;  // GMT offset in hours
    scalar standardMeridian = timeZoneOffset * 15.0;  // Convert to degrees
    
    // Time correction from standard meridian to actual longitude
    // 4 minutes per degree of longitude
    scalar longitudeCorrection = 4.0 * (longitude - standardMeridian);
    
    // Solar time = local standard time + equation of time + longitude correction
    scalar solarTime = hour + (eqTimeMinutes + longitudeCorrection) / 60.0;
    
    // Hour angle from solar noon (negative = morning, positive = afternoon)
    scalar hourAngle = (solarTime - 12.0) * 15.0 * M_PI/180.0;
    scalar latRad = latitude * M_PI/180.0;
    
    // Solar elevation angle
    scalar solarElevation = Foam::asin(Foam::sin(latRad) * Foam::sin(declination) + 
                                Foam::cos(latRad) * Foam::cos(declination) * Foam::cos(hourAngle));
    
    if (showDebug)
    {
        Info << "Solar calculation debug:" << endl;
        Info << "  Location: Lat " << latitude << ", Lon " << longitude << endl;
        Info << "  Time zone: GMT" << (timeZoneOffset >= 0 ? "+" : "") << timeZoneOffset << endl;
        Info << "  Standard meridian: " << standardMeridian << " degrees" << endl;
        Info << "  Equation of time: " << eqTimeMinutes << " minutes" << endl;
        Info << "  Longitude correction: " << longitudeCorrection << " minutes" << endl;
        Info << "  Hour: " << hour << ", Solar time: " << solarTime << endl;
        Info << "  Hour angle: " << (hourAngle * 180.0/M_PI) << " degrees" << endl;
        Info << "  Solar elevation: " << (solarElevation * 180.0/M_PI) << " degrees" << endl;
    }
    
    if (solarElevation <= 0) 
    {
        if (showDebug) Info << "Solar elevation <= 0, sun below horizon. Returning air temp." << endl;
        return airTemp; // No sun
    }
    
    // Solar azimuth (needed for directional effects)
    // Note: solarAzimuth calculation commented out as it's not currently used
    // scalar solarAzimuth = Foam::atan2(Foam::sin(hourAngle), 
    //                            Foam::cos(hourAngle) * Foam::sin(latRad) - Foam::tan(declination) * Foam::cos(latRad));
    
    // Get radiation data from EPW for the specific hour
    scalar directNormalIrradiance = 0.0;
    scalar diffuseHorizontalIrradiance = 0.0;
    scalar globalHorizontalIrradiance = 0.0;
    
    // Calculate hour index in EPW data (0-8759)
    // Note: hour should be 0-23, not 1-24
    label hourOfYear = label((dayOfYear - 1) * 24 + hour);
    if (showDebug)
    {
        Info << "\n=== Solar MRT Calculation Debug ===" << endl;
        Info << "Hour of year index: " << hourOfYear << " (day " << dayOfYear << ", hour " << hour << ")" << endl;
    }
    if (hourOfYear >= 0 && hourOfYear < 8760)
    {
        // EPW files store radiation in Wh/m² (energy), not W/m² (power)
        // Since these are hourly integrated values, they are already in W/m² (1 hour average)
        globalHorizontalIrradiance = epwData.globalRadiation[hourOfYear];
        directNormalIrradiance = epwData.directRadiation[hourOfYear];
        diffuseHorizontalIrradiance = epwData.diffuseRadiation[hourOfYear];
        
        // Debug output
        if (showDebug)
        {
            Info << "EPW radiation data for hour " << hour << " on day " << dayOfYear << ":" << endl;
            Info << "  Global horizontal: " << globalHorizontalIrradiance << " W/m2" << endl;
            Info << "  Direct normal: " << directNormalIrradiance << " W/m2" << endl;
            Info << "  Diffuse horizontal: " << diffuseHorizontalIrradiance << " W/m2" << endl;
        }
    }
    else
    {
        WarningInFunction
            << "Hour index " << hourOfYear << " out of range for EPW data" << endl;
    }
    
    scalar skyViewFactor = 0.5;
    scalar groundViewFactor = 0.5;
    
    scalar skyTemp = airTemp - 10.0;
    scalar groundTemp = airTemp + (directNormalIrradiance > 0 ? 5.0 : 0.0);
    
    scalar solarAbsorptivity = 0.7;
    scalar effectiveArea = 0.77;
    
    // Calculate solar heat gain on the person
    scalar solarHeatGain = solarAbsorptivity * effectiveArea * 
                          (directNormalIrradiance * Foam::sin(solarElevation) + 
                           diffuseHorizontalIrradiance * skyViewFactor);
    
    // Debug output for solar heat gain
    if (showDebug)
    {
        Info << "  Solar heat gain: " << solarHeatGain << " W/m²" << endl;
        Info << "  Solar elevation: " << (solarElevation * 180.0/M_PI) << " degrees" << endl;
    }
    
    // Calculate equivalent temperature rise from solar radiation
    // Using empirical correlation for outdoor thermal comfort
    // Typical values: 100 W/m² solar gain ≈ 2-3C MRT increase
    scalar solarTempRise = 0.0;
    if (solarHeatGain > 0)
    {
        // Empirical correlation for outdoor thermal comfort
        // Based on Thorsson et al. (2007) and similar studies
        // MRT rise approximately 0.25°C per 10 W/m² absorbed solar radiation
        solarTempRise = 0.025 * solarHeatGain;
        
        // Limit solar temperature rise to realistic values
        solarTempRise = min(solarTempRise, 25.0);  // Max 25C rise from solar
    }
    
    if (showDebug) Info << "  Solar temperature rise: " << solarTempRise << " C" << endl;
    
    // Calculate mean radiant temperature including solar effects
    scalar meanRadiantTemp = skyViewFactor * skyTemp + 
                            groundViewFactor * groundTemp + 
                            solarTempRise;
    
    if (showDebug) Info << "  Final MRT: " << meanRadiantTemp << " C (Air temp: " << airTemp << " C)" << endl;
    
    return meanRadiantTemp;
}

// Read EPW file and return complete weather data
EPWData readEPWFile(const Foam::fileName& epwFile)
{
    EPWData data;
    data.temperature.setSize(8760);    // 365 days * 24 hours
    data.globalRadiation.setSize(8760);
    data.directRadiation.setSize(8760);
    data.diffuseRadiation.setSize(8760);
    
    std::ifstream file(epwFile);
    if (!file.is_open())
    {
        FatalErrorInFunction
            << "Cannot open EPW file: " << epwFile
            << exit(FatalError);
    }
    
    // Parse EPW header to extract location data
    std::string line;
    
    // Line 1: LOCATION header contains location data
    if (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string token;
        std::vector<std::string> tokens;
        
        // Split by comma
        while (std::getline(iss, token, ','))
        {
            tokens.push_back(token);
        }
        
        // EPW format: LOCATION,City,State,Country,Source,WMO,Latitude,Longitude,TimeZone,Elevation
        if (tokens.size() >= 10 && tokens[0] == "LOCATION")
        {
            data.locationName = tokens[1] + ", " + tokens[3];  // City, Country
            data.latitude = std::stod(tokens[6]);
            data.longitude = std::stod(tokens[7]);
            data.timeZone = std::stod(tokens[8]);
            
            Info << "EPW Location Data:" << endl;
            Info << "  Location: " << data.locationName << endl;
            Info << "  Latitude: " << data.latitude << " degrees" << endl;
            Info << "  Longitude: " << data.longitude << " degrees" << endl;
            Info << "  Time Zone: GMT" << (data.timeZone >= 0 ? "+" : "") << data.timeZone << endl;
        }
        else
        {
            WarningInFunction
                << "Could not parse location data from EPW header" << endl;
            // Set default values
            data.latitude = 50.0;
            data.longitude = 8.0;
            data.timeZone = 1.0;
            data.locationName = "Unknown Location";
        }
    }
    
    // Skip remaining header lines (2-8)
    for (int i = 1; i < 8; i++)
    {
        std::getline(file, line);
    }
    
    label hourIndex = 0;
    while (std::getline(file, line) && hourIndex < 8760)
    {
        std::istringstream iss(line);
        std::string token;
        
        // Skip to temperature (column 7)
        for (int i = 0; i < 6; i++)
        {
            std::getline(iss, token, ',');
        }
        
        // Column 7: Dry bulb temperature
        std::getline(iss, token, ',');
        data.temperature[hourIndex] = std::stod(token);
        
        // Skip to radiation data (columns 13-15)
        for (int i = 0; i < 5; i++)
        {
            std::getline(iss, token, ',');
        }
        
        // Column 13: Global horizontal radiation
        std::getline(iss, token, ',');
        data.globalRadiation[hourIndex] = std::stod(token);
        
        // Column 14: Direct normal radiation
        std::getline(iss, token, ',');
        data.directRadiation[hourIndex] = std::stod(token);
        
        // Column 15: Diffuse horizontal radiation
        std::getline(iss, token, ',');
        data.diffuseRadiation[hourIndex] = std::stod(token);
        
        hourIndex++;
    }
    
    file.close();
    return data;
}

// Improved EPW parser with robust day-of-year calculation
Foam::scalar calculateRunningMeanFromEPW
(
    const Foam::fileName& epwFile,
    const Foam::scalar& dayOfYear
)
{
    std::ifstream file(epwFile);
    if (!file.is_open())
    {
        FatalErrorInFunction
            << "Cannot open EPW file: " << epwFile
            << exit(FatalError);
    }

    // Skip EPW header (8 lines)
    std::string line;
    for (int i = 0; i < 8; i++)
    {
        std::getline(file, line);
    }

    // Read hourly data and calculate daily averages
    List<scalar> dailyTemps(365, 0.0);
    List<label> hourCount(365, 0);
    
    // Days in each month (non-leap year)
    int daysInMonth[] = {31,28,31,30,31,30,31,31,30,31,30,31};

    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string token;
        
        // Parse CSV: year,month,day,hour,minute,datasource,temp,dewpoint,...
        std::getline(iss, token, ','); // year
        label year = std::stoi(token);
        std::getline(iss, token, ','); // month
        label month = std::stoi(token);
        std::getline(iss, token, ','); // day  
        label day = std::stoi(token);
        std::getline(iss, token, ','); // hour
        std::getline(iss, token, ','); // minute
        std::getline(iss, token, ','); // datasource
        std::getline(iss, token, ','); // dry bulb temperature
        scalar temp = std::stod(token);
        
        // Calculate day of year properly
        label doy = 0;
        
        // Check for leap year
        bool isLeapYear = (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
        if (isLeapYear && month > 2) doy += 1; // Add leap day if after February
        
        // Sum days of previous months
        for (int m = 0; m < month-1; m++) {
            doy += daysInMonth[m];
        }
        doy += day - 1; // Add current day (0-indexed)
        
        if (doy < 365)
        {
            dailyTemps[doy] += temp;
            hourCount[doy]++;
        }
    }

    // Average daily temperatures
    label validDays = 0;
    forAll(dailyTemps, dayI)
    {
        if (hourCount[dayI] > 0)
        {
            dailyTemps[dayI] /= hourCount[dayI];
            validDays++;
        }
        else
        {
            // Handle missing days by interpolation
            if (dayI > 0 && dayI < 364)
            {
                // Find nearest valid days
                label prevDay = dayI - 1;
                label nextDay = dayI + 1;
                while (prevDay > 0 && hourCount[prevDay] == 0) prevDay--;
                while (nextDay < 365 && hourCount[nextDay] == 0) nextDay++;
                
                if (hourCount[prevDay] > 0 && hourCount[nextDay] > 0)
                {
                    // Linear interpolation
                    scalar weight = scalar(dayI - prevDay) / scalar(nextDay - prevDay);
                    dailyTemps[dayI] = dailyTemps[prevDay] * (1.0 - weight) + 
                                       dailyTemps[nextDay] * weight;
                }
            }
        }
    }
    
    Info << "EPW file processed: " << validDays << " valid days found" << endl;

    // Calculate running mean (ASHRAE 55 formula)
    scalar alpha = 0.8;
    scalar runningMean = 0.0;
    scalar totalWeight = 0.0;

    label targetDay = label(dayOfYear) - 1;

    for (int i = 1; i <= 30; i++)
    {
        label pastDay = (targetDay - i + 365) % 365;
        scalar weight = alpha * Foam::pow(1.0 - alpha, i - 1);
        runningMean += weight * dailyTemps[pastDay];
        totalWeight += weight;
    }

    scalar T_rm = runningMean / totalWeight;
    
    // Sanity check
    if (T_rm < -50 || T_rm > 50)
    {
        WarningInFunction
            << "Calculated running mean temperature (" << T_rm 
            << "C) seems unrealistic!" << endl;
    }

    return T_rm;
}

void validatePMVInputs(scalar met, scalar clo, scalar wme, scalar rh)
{
    if (met < 0.8 || met > 4.0)
    {
        FatalErrorInFunction
            << "Metabolic rate out of valid range (0.8-4.0 met): " << met
            << exit(FatalError);
    }

    if (clo < 0.0 || clo > 2.0)
    {
        FatalErrorInFunction
            << "Clothing insulation out of valid range (0-2.0 clo): " << clo
            << exit(FatalError);
    }

    if (wme < 0.0 || wme > met)
    {
        FatalErrorInFunction
            << "External work out of valid range (0-" << met << " met): " << wme
            << exit(FatalError);
    }

    if (rh < 0.0 || rh > 100.0)
    {
        FatalErrorInFunction
            << "Relative humidity out of valid range (0-100%): " << rh
            << exit(FatalError);
    }
}

Foam::scalar calculateWaterVapourPressure(scalar temperature, scalar relativeHumidity)
{
    scalar tempCelsius = temperature - 273.0;
    return relativeHumidity * 10.0
         * Foam::exp(16.6536 - (4030.183 / (tempCelsius + 235.0)));
}

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
    using namespace Foam::ASHRAE55Constants;

    scalar airTempC = airTemp - 273.0;
    scalar tcla = airTemp + (35.5 - airTempC) / (3.5 * 6.45 * (icl + 0.1));
    scalar xn = tcla / 100.0;
    scalar xf = xn;

    scalar p1 = icl * fcl;
    scalar p2 = p1 * radiationCoeff;
    scalar p3 = p1 * 100.0;
    scalar p4 = p1 * airTemp;
    scalar p5 = 308.7 - 0.028 * (metabolicRate - externalWork)
              + p2 * Foam::pow((radiationTemp + 273.0) / 100.0, 4);

    label iterationCount = 0;

    do
    {
        iterationCount++;
        xf = (xf + xn) / 2.0;

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
    }
    while (mag(xn - xf) > clothingConvergenceTol);

    return Tuple2<scalar>(100.0 * xn - 273.0, xn);
}

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
    using namespace Foam::ASHRAE55Constants;

    scalar airTempC = airTemp - 273.0;
    scalar icl = 0.155 * clothingInsulation;
    scalar fcl = (icl < clothingThreshold)
        ? (1.0 + clothingFactor1 * icl)
        : (clothingFactor2 + clothingFactor3 * icl);

    scalar pa = calculateWaterVapourPressure(airTemp, relativeHumidity);

    Tuple2<scalar> tclResult = calculateClothingSurfaceTemperature
    (
        airTemp,
        velocity,
        icl,
        fcl,
        radiationTemp,
        metabolicRate,
        externalWork
    );

    scalar tcl = tclResult.first();
    scalar xn = tclResult.second();

    scalar hl1 = skinDiffusionCoeff * (basePressure - (thermalCoeff * (metabolicRate - externalWork)) - pa);
    scalar hl2 = 0.0;

    if ((metabolicRate - externalWork) > baseMetabolicRate)
    {
        hl2 = 0.42 * ((metabolicRate - externalWork) - baseMetabolicRate);
    }

    scalar hl3 = respirationCoeff1 * metabolicRate * (respirationHumidity - pa);
    scalar hl4 = respirationCoeff2 * metabolicRate * (respirationTemp - airTempC);
    scalar hl5 = radiationCoeff * fcl
        * (Foam::pow(xn, 4) - Foam::pow((radiationTemp + 273.0) / 100.0, 4));

    scalar hcf = forcedConvectionCoeff * Foam::sqrt(velocity);
    scalar hcn = naturalConvectionCoeff * Foam::pow(mag(tcl - airTempC), naturalConvectionExp);
    scalar hc = Foam::max(hcf, hcn);
    scalar hl6 = fcl * hc * (tcl - airTempC);

    scalar ts = thermalSensCoeff1 * Foam::exp(-thermalSensCoeff2 * metabolicRate) + thermalSensCoeff3;

    return ts * ((metabolicRate - externalWork) - hl1 - hl2 - hl3 - hl4 - hl5 - hl6);
}

Foam::scalar calculatePPD(scalar pmv)
{
    using namespace Foam::ASHRAE55Constants;

    return 100.0 - ppdBase * Foam::exp(-ppdCoeff1 * Foam::pow(pmv, 4) - ppdCoeff2 * Foam::pow(pmv, 2));
}

Foam::autoPtr<volScalarField> loadRelativeHumidityField
(
    const fvMesh& mesh,
    const Time& runTime,
    bool& humidityAvailable,
    word& humiditySource
)
{
    wordList humidityCandidates(5);
    humidityCandidates[0] = "thermo:relHum";
    humidityCandidates[1] = "thermoRelHum";
    humidityCandidates[2] = "relHum";
    humidityCandidates[3] = "RH";
    humidityCandidates[4] = "relativeHumidity";

    forAll(humidityCandidates, i)
    {
        IOobject humidityHeader
        (
            humidityCandidates[i],
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        );

        if (humidityHeader.typeHeaderOk<volScalarField>(true))
        {
            humidityAvailable = true;
            humiditySource = humidityCandidates[i];
            return autoPtr<volScalarField>(new volScalarField(humidityHeader, mesh));
        }
    }

    humidityAvailable = false;
    humiditySource = "ASHRAE55Dict RH";
    return autoPtr<volScalarField>();
}

// Function to calculate and display comfort statistics
void displayComfortStatistics
(
    const volScalarField& ASHRAELevel80,
    const volScalarField& ASHRAELevel90,
    const volScalarField& TOp,
    const volScalarField& T,
    const volVectorField& U,
    const scalar& T_rm_out
)
{
    // Calculate statistics
    label cells80 = 0;
    label cells90 = 0;
    label totalCells = 0;
    
    scalar minTOp = 1000.0;
    scalar maxTOp = -1000.0;
    scalar avgTOp = 0.0;
    scalar avgAirTemp = 0.0;
    scalar avgVelocity = 0.0;
    
    forAll(ASHRAELevel80, cellI)
    {
        if (T[cellI] > 283.15 && T[cellI] < 306.65) // Valid temperature range
        {
            totalCells++;
            avgTOp += TOp[cellI];
            avgAirTemp += T[cellI];
            avgVelocity += mag(U[cellI]);
            
            if (TOp[cellI] < minTOp) minTOp = TOp[cellI];
            if (TOp[cellI] > maxTOp) maxTOp = TOp[cellI];
            
            if (ASHRAELevel80[cellI] > 0.5) cells80++;
            if (ASHRAELevel90[cellI] > 0.5) cells90++;
        }
    }
    
    if (totalCells > 0)
    {
        avgTOp /= totalCells;
        avgAirTemp /= totalCells;
        avgVelocity /= totalCells;
        
        Info<< nl << "================== ASHRAE 55 Comfort Analysis ==================" << endl;
        Info<< "Running mean outdoor temperature: " << T_rm_out << " C" << endl;
        Info<< "Neutral temperature (t_cmf): " << (0.31 * T_rm_out + 17.8) << " C" << endl;
        Info<< nl << "Domain Statistics:" << endl;
        Info<< "  Total cells evaluated: " << totalCells << endl;
        Info<< "  Average air temperature: " << (avgAirTemp - 273.15) << " C" << endl;
        Info<< "  Average operative temperature: " << (avgTOp - 273.15) << " C" << endl;
        Info<< "  Operative temperature range: " << (minTOp - 273.15) 
            << " to " << (maxTOp - 273.15) << " C" << endl;
        Info<< "  Average air velocity: " << avgVelocity << " m/s" << endl;
        Info<< nl << "Comfort Compliance:" << endl;
        Info<< "  80% acceptability: " << cells80 << " cells (" 
            << (scalar(cells80)/totalCells*100.0) << "%)" << endl;
        Info<< "  90% acceptability: " << cells90 << " cells (" 
            << (scalar(cells90)/totalCells*100.0) << "%)" << endl;
        
        // Comfort zone ranges
        scalar t_cmf = 0.31 * T_rm_out + 17.8;
        Info<< nl << "Comfort Zone Ranges (without cooling effect):" << endl;
        Info<< "  80% acceptability: " << (t_cmf - 3.5) << " to " 
            << (t_cmf + 3.5) << " C" << endl;
        Info<< "  90% acceptability: " << (t_cmf - 2.5) << " to " 
            << (t_cmf + 2.5) << " C" << endl;
        
        // Check if cooling effect is active
        label cellsWithCooling = 0;
        forAll(U, cellI)
        {
            if (mag(U[cellI]) >= 0.6 && TOp[cellI] > 298.15) cellsWithCooling++;
        }
        
        if (cellsWithCooling > 0)
        {
            Info<< nl << "Cooling effect active in " << cellsWithCooling 
                << " cells (" << (scalar(cellsWithCooling)/totalCells*100.0) 
                << "%)" << endl;
        }
        
        Info<< "================================================================" << nl << endl;
    }
    else
    {
        Info<< nl << "Warning: No cells in valid temperature range for comfort analysis!" << endl;
    }
}

void displayPMVStatistics
(
    const volScalarField& ASHRAEAcceptable,
    const volScalarField& PMV,
    const volScalarField& PPD,
    const volScalarField& TOp,
    const volScalarField& T,
    const volVectorField& U,
    scalar averageRH
)
{
    label acceptableCells = 0;
    label totalCells = 0;
    scalar avgPMV = 0.0;
    scalar avgPPD = 0.0;
    scalar avgTOp = 0.0;
    scalar avgAirTemp = 0.0;
    scalar avgVelocity = 0.0;

    forAll(PMV, cellI)
    {
        totalCells++;
        avgPMV += PMV[cellI];
        avgPPD += PPD[cellI];
        avgTOp += TOp[cellI];
        avgAirTemp += T[cellI];
        avgVelocity += mag(U[cellI]);

        if (ASHRAEAcceptable[cellI] > 0.5)
        {
            acceptableCells++;
        }
    }

    if (totalCells > 0)
    {
        avgPMV /= totalCells;
        avgPPD /= totalCells;
        avgTOp /= totalCells;
        avgAirTemp /= totalCells;
        avgVelocity /= totalCells;

        Info<< nl << "================= ASHRAE 55 PMV Analysis =================" << endl;
        Info<< "Acceptability criterion: -0.5 <= PMV <= +0.5" << endl;
        Info<< nl << "Domain Statistics:" << endl;
        Info<< "  Total cells evaluated: " << totalCells << endl;
        Info<< "  Average air temperature: " << (avgAirTemp - 273.15) << " C" << endl;
        Info<< "  Average operative temperature: " << (avgTOp - 273.15) << " C" << endl;
        Info<< "  Average air velocity: " << avgVelocity << " m/s" << endl;
        Info<< "  Relative humidity used: " << averageRH << " %" << endl;
        Info<< nl << "Comfort Results:" << endl;
        Info<< "  Average PMV: " << avgPMV << endl;
        Info<< "  Average PPD: " << avgPPD << " %" << endl;
        Info<< "  Acceptable cells: " << acceptableCells << " ("
            << (scalar(acceptableCells)/totalCells*100.0) << "%)" << endl;
        Info<< "==========================================================" << nl << endl;
    }
}

// * * * * * * * * * * * * * * * * * * Program  * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    Foam::timeSelector::addOptions();

    // Add command line options
    argList::addOption
    (
        "epw",
        "fileName",
        "EPW weather file to calculate running mean outdoor temperature"
    );

    argList::addOption
    (
        "dayOfYear", 
        "scalar",
        "Day of year (1-365) for EPW calculation (default: 180 = June 29)"
    );

    argList::addBoolOption
    (
        "solarData",
        "Use detailed solar calculations from EPW file for accurate MRT (requires -epw option)"
    );

    argList::addOption
    (
        "latitude",
        "scalar", 
        "Site latitude in degrees (overrides EPW file value if provided)"
    );

    argList::addOption
    (
        "longitude", 
        "scalar",
        "Site longitude in degrees (overrides EPW file value if provided)"
    );

    argList::addOption
    (
        "hour",
        "scalar",
        "Hour of day for solar calculations (0-24, default: 12.0)"
    );

    argList::addOption
    (
        "runningMean",
        "scalar", 
        "Directly specify running mean outdoor temperature in C"
    );

    argList::addOption
    (
        "model",
        "word",
        "ASHRAE 55 model: adaptive or pmv (default: adaptive)"
    );

    argList::addOption
    (
        "met",
        "scalar",
        "Metabolic rate in met for PMV mode"
    );

    argList::addOption
    (
        "clo",
        "scalar",
        "Clothing insulation in clo for PMV mode"
    );

    argList::addOption
    (
        "wme",
        "scalar",
        "External work in met for PMV mode"
    );

    argList::addOption
    (
        "rh",
        "scalar",
        "Relative humidity in percent for PMV mode if no humidity field is available"
    );

    #include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"

    //- Get times list
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"

    IOdictionary ASHRAE55Dict
    (
        IOobject
        (
            "ASHRAE55Dict",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    word modelName("adaptive");
    scalar met = 1.2;
    scalar clo = 0.5;
    scalar wme = 0.0;
    scalar RH = 50.0;

    if (ASHRAE55Dict.found("model"))
    {
        modelName = word(ASHRAE55Dict.lookup("model"));
    }
    if (ASHRAE55Dict.found("met"))
    {
        met = readScalar(ASHRAE55Dict.lookup("met"));
    }
    if (ASHRAE55Dict.found("clo"))
    {
        clo = readScalar(ASHRAE55Dict.lookup("clo"));
    }
    if (ASHRAE55Dict.found("wme"))
    {
        wme = readScalar(ASHRAE55Dict.lookup("wme"));
    }
    if (ASHRAE55Dict.found("RH"))
    {
        RH = readScalar(ASHRAE55Dict.lookup("RH"));
    }

    if (args.found("model"))
    {
        modelName = args.get<word>("model");
    }
    if (args.found("met"))
    {
        met = args.get<scalar>("met");
    }
    if (args.found("clo"))
    {
        clo = args.get<scalar>("clo");
    }
    if (args.found("wme"))
    {
        wme = args.get<scalar>("wme");
    }
    if (args.found("rh"))
    {
        RH = args.get<scalar>("rh");
    }

    if (modelName != "adaptive" && modelName != "pmv")
    {
        FatalErrorInFunction
            << "Unknown ASHRAE 55 model '" << modelName
            << "'. Supported values are 'adaptive' and 'pmv'."
            << exit(FatalError);
    }

    if (modelName == "pmv")
    {
        validatePMVInputs(met, clo, wme, RH);
        Info << "Using ASHRAE 55 PMV model" << endl;
        Info << "  met = " << met << ", clo = " << clo
             << ", wme = " << wme << ", RH = " << RH << "%" << endl;
    }
    else
    {
        Info << "Using ASHRAE 55 adaptive model" << endl;
    }

    // Determine running mean outdoor temperature and solar calculation settings
    scalar T_rm_out = 20.0; // Default value
    bool useSolarCalculations = args.found("solarData");
    scalar latitude = 50.0;  // Default: Central Europe
    scalar longitude = 8.0;  // Default: Central Europe
    scalar hourOfDay = 12.0; // Default: noon
    fileName epwFile;
    scalar dayOfYear = 180;
    EPWData epwData;  // Store complete EPW weather data
    bool useEPWCoordinates = true; // Flag to track coordinate source

    if (args.found("hour"))
        hourOfDay = args.get<scalar>("hour");

    if (args.found("runningMean"))
    {
        T_rm_out = args.get<scalar>("runningMean");
        Info << "Using directly specified running mean: " << T_rm_out << " C" << endl;
        
        // If no EPW file, must use command line coordinates
        if (args.found("latitude"))
            latitude = args.get<scalar>("latitude");
        if (args.found("longitude"))
            longitude = args.get<scalar>("longitude");
        useEPWCoordinates = false;
    }
    else if (args.found("epw"))
    {
        epwFile = args.get<fileName>("epw");
        
        if (args.found("dayOfYear"))
            dayOfYear = args.get<scalar>("dayOfYear");
            
        Info << "Reading EPW file: " << epwFile << endl;
        Info << "Day of year: " << dayOfYear << endl;
        
        // Read complete EPW data including radiation and location
        epwData = readEPWFile(epwFile);
        
        // Use EPW coordinates by default, but allow command line override
        latitude = epwData.latitude;
        longitude = epwData.longitude;
        
        // Check for command line overrides
        if (args.found("latitude"))
        {
            latitude = args.get<scalar>("latitude");
            useEPWCoordinates = false;
            Info << "Overriding EPW latitude with command line value: " << latitude << endl;
        }
        if (args.found("longitude"))
        {
            longitude = args.get<scalar>("longitude");
            useEPWCoordinates = false;
            Info << "Overriding EPW longitude with command line value: " << longitude << endl;
        }
        
        T_rm_out = calculateRunningMeanFromEPW(epwFile, dayOfYear);
        Info << "Calculated running mean: " << T_rm_out << " C" << endl;
        
        if (useSolarCalculations)
        {
            Info << nl << "Solar calculations enabled" << endl;
            Info << "Using coordinates from " << (useEPWCoordinates ? "EPW file" : "command line") << ":" << endl;
            Info << "  Latitude: " << latitude << " degrees" << endl;
            Info << "  Longitude: " << longitude << " degrees" << endl;
            if (useEPWCoordinates)
            {
                Info << "  Location: " << epwData.locationName << endl;
            }
            Info << "Note: For best results, use OpenFOAM 2412+ native solar radiation model instead" << endl;
        }
    }
    else if (modelName == "adaptive")
    {
        Info << "Warning: No EPW file or running mean specified." << endl;
        Info << "Using default value: " << T_rm_out << " C" << endl;
        Info << "Usage examples:" << endl;
        Info << "  # With OpenFOAM solar model (G contains solar radiation):" << endl;
        Info << "  ASHRAE55 -epw weather.epw -dayOfYear 180" << endl;
        Info << "  ASHRAE55 -runningMean 22.5" << endl;
        Info << "  # With EPW solar calculations (auto-detects location from EPW):" << endl;
        Info << "  ASHRAE55 -epw weather.epw -solarData -hour 14" << endl;
        Info << "  # Override EPW coordinates if needed:" << endl; 
        Info << "  ASHRAE55 -epw weather.epw -solarData -latitude 52.5 -longitude 13.4" << endl;
        
        if (useSolarCalculations)
        {
            FatalErrorInFunction
                << "Solar calculations require EPW file (-epw option)"
                << exit(FatalError);
        }
    }
    else if (useSolarCalculations)
    {
        FatalErrorInFunction
            << "Solar calculations require EPW file (-epw option)"
            << exit(FatalError);
    }

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        Info<< "Time = " << runTime.timeName() << endl;

        #include "createFields.H"

        volScalarField T(THeader, mesh);
        const fvPatchList& Patches = T.mesh().boundary();

        scalar STemp(20);

        volVectorField U(UHeader, mesh);

        scalar t_cmf(0), ce(0);
        scalar averageRHUsed = 0.0;
        bool humidityAvailable = false;
        word humiditySource("ASHRAE55Dict RH");
        autoPtr<volScalarField> humidityField;

        if (modelName == "pmv")
        {
            humidityField = loadRelativeHumidityField(mesh, runTime, humidityAvailable, humiditySource);
            Info << "PMV mode humidity source: " << humiditySource;
            if (!humidityAvailable)
            {
                Info << " (" << RH << "% fallback)";
            }
            Info << endl;
        }

        //- Check if radiation field G is available
        Info << "G field check: G.headerOk() = " << G.headerOk() 
             << ", useSolarCalculations = " << useSolarCalculations 
             << ", epwFile.size() = " << epwFile.size() << endl;
        
        //- Radiation Model not available? Use area-weighted wall temperature
        if (G.headerOk()!=1 && !useSolarCalculations)
        {
            STemp = radiationTemperature(mesh, Patches);
            Info << "Using area-weighted wall temperature for MRT: " << STemp << " C" << endl;
        }

        forAll (mesh.cells(), cellI)
        {
            ASHRAELevel80[cellI] = 0;
            ASHRAELevel90[cellI] = 0;
            PMV[cellI] = 0;
            PPD[cellI] = 0;
            ASHRAEAcceptable[cellI] = 0;
            
            //- Use radiation model if available
            if (G.headerOk() == 1)
            {
                if (G[cellI] < 0)
                    G[cellI] = 0;

                //- Limit maximum radiation to avoid numerical issues
                if ( G[cellI] > 50000 )
                    G[cellI] = 50000;

                //- For solar radiation: Convert to realistic mean radiant temperature
                //- CRITICAL: Stefan-Boltzmann inversion gives unrealistic values for solar irradiance
                //- Use empirical correlations based on outdoor thermal comfort measurements
                
                scalar T_radiant_equiv = Foam::pow( G[cellI] / ( 4.0 * 5.67e-8), 0.25) - 273.15;
                
                if (T_radiant_equiv > 60.0)  // Unrealistic for outdoor MRT
                {
                    //- Use empirical solar effect correlation
                    //- Based on measured outdoor MRT vs solar irradiance data
                    //- Typical outdoor MRT = air temperature + solar effect (5-20C)
                    scalar solar_effect = 5.0 + Foam::sqrt(G[cellI] / 100.0);  // Empirical: 5-20C rise
                    solar_effect = min(solar_effect, 20.0);  // Max 20C solar heating
                    
                    STemp = (T[cellI] - 273.15) + solar_effect;
                }
                else
                {
                    //- Use Stefan-Boltzmann result for low radiation
                    STemp = T_radiant_equiv;
                }
                
                //- Final safety limits for extreme conditions
                STemp = min(STemp, 65.0);   // Max MRT = 65C (hot asphalt limit)
                STemp = max(STemp, -30.0);  // Min MRT = -30C (extreme cold)
            }
            //- Use EPW solar calculations if radiation model not available and -solarData specified
            else if (G.headerOk() != 1 && useSolarCalculations && epwFile.size() > 0)
            {
                if (cellI == 0)
                {
                    Info << "Entering EPW solar calculation branch" << endl;
                }
                //- Use the specified hour of day with actual EPW radiation data
                bool showDebugInfo = (cellI == 0);  // Only show debug for first cell
                STemp = calculateSolarMRT(
                    epwData,
                    dayOfYear,
                    hourOfDay,  // Use command line specified hour
                    latitude,
                    longitude,
                    T[cellI] - 273.15,  // Air temperature in Celsius
                    showDebugInfo
                );
                
                if (cellI == 0)
                {
                    Info<< "Solar MRT calculation: Air temp = " << (T[cellI] - 273.15) 
                        << " C, Solar MRT = " << STemp << " C" << endl;
                }
            }
            else
            {
                // No radiation field and no solar calculations - use air temperature as MRT
                STemp = T[cellI] - 273.15;
                if (cellI == 0)
                {
                    Info << "No radiation data available, using air temperature as MRT" << endl;
                }
            }

            //- Calculate operative temperature: average of air and mean radiant temperature
            TOp[cellI] = ( T[cellI] + (STemp + 273.15)) / 2.0;

            if (modelName == "adaptive")
            {
                //- Only evaluate comfort if temperature is in reasonable range (10C to 33.5C)
                if ( (T[cellI] > 283.15) && (T[cellI] < 306.65) )
                {
                    ce = 0;

                    //- Calculate cooling effect of elevated air speed when Top > 25C
                    if ( (mag(U[cellI]) >= 0.6) && (TOp[cellI] > 298.15) )
                    {
                        if (mag(U[cellI]) < 0.9)
                        {
                            ce = 1.2;
                        }
                        else if (mag(U[cellI]) < 1.2)
                        {
                            ce = 1.8;
                        }
                        else
                        {
                            ce = 2.2;
                        }
                    }

                    t_cmf = (0.31 * T_rm_out) + 17.8;
                    scalar t_cmf_K = t_cmf + 273.15;

                    if ((TOp[cellI] >= (t_cmf_K - 3.5)) && (TOp[cellI] <= (t_cmf_K + 3.5 + ce)))
                        ASHRAELevel80[cellI] = 1;

                    if ((TOp[cellI] >= (t_cmf_K - 2.5)) && (TOp[cellI] <= (t_cmf_K + 2.5 + ce)))
                        ASHRAELevel90[cellI] = 1;
                }
            }
            else
            {
                scalar cellRH = RH;

                if (humidityAvailable && humidityField.valid())
                {
                    scalar rawHumidity = humidityField()[cellI];
                    cellRH = (rawHumidity <= 1.0) ? rawHumidity * 100.0 : rawHumidity;
                    cellRH = Foam::max(0.0, Foam::min(100.0, cellRH));
                }

                scalar pmv = calculatePMV
                (
                    T[cellI],
                    mag(U[cellI]),
                    cellRH,
                    STemp,
                    met * Foam::ASHRAE55Constants::baseMetabolicRate,
                    clo,
                    wme * Foam::ASHRAE55Constants::baseMetabolicRate
                );

                PMV[cellI] = pmv;
                PPD[cellI] = calculatePPD(pmv);

                if (pmv >= -0.5 && pmv <= 0.5)
                {
                    ASHRAEAcceptable[cellI] = 1;
                }

                averageRHUsed += cellRH;
            }
        }

        if (modelName == "adaptive")
        {
            displayComfortStatistics(ASHRAELevel80, ASHRAELevel90, TOp, T, U, T_rm_out);
            ASHRAELevel80.write();
            ASHRAELevel90.write();
        }
        else
        {
            averageRHUsed =
                (mesh.cells().size() > 0)
              ? averageRHUsed / scalar(mesh.cells().size())
              : RH;
            displayPMVStatistics(ASHRAEAcceptable, PMV, PPD, TOp, T, U, averageRHUsed);
            PMV.write();
            PPD.write();
            ASHRAEAcceptable.write();
        }

        TOp.write();

        Info << "Done" << endl;
    }

    return 0;
}


// ************************************************************************* //
