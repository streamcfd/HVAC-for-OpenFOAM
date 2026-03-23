# ASHRAE55Foam

## Overview

ASHRAE55Foam is an OpenFOAM post-processing utility that calculates thermal comfort metrics according to ANSI/ASHRAE Standard 55-2020. It supports two ASHRAE 55 pathways in one standalone tool:

- `adaptive` for naturally ventilated spaces, based on running mean outdoor temperature
- `pmv` for conditioned indoor spaces, based on PMV/PPD with ASHRAE 55 acceptability output

The utility evaluates thermal comfort from CFD results, computes operative temperature, and can extract location data from EPW weather files when the adaptive model is used.

## Features

- **Global Location Support**
  - Automatically extracts latitude, longitude, and time zone from EPW files
  - Works for any location worldwide without manual coordinate input
  - Supports command-line override of location parameters

- **Thermal Comfort Assessment**
  - Calculates Operative Temperature (TOp) - the average of air temperature and mean radiant temperature
  - Evaluates compliance with ASHRAE 55 adaptive comfort model at 80% and 90% acceptability levels
  - Evaluates PMV-based indoor acceptability for conditioned spaces
  - Accounts for cooling effects from elevated air speeds (>0.6 m/s)

- **Solar Radiation Integration**
  - Optimized for OpenFOAM 2412+ with native solar radiation models
  - Automatically reads radiation field 'G' from OpenFOAM radiation models (solarLoad, P1, fvDOM)
  - Includes fallback EPW-based solar calculations with proper time zone corrections

- **Climate Data Processing**
  - Parses EPW (EnergyPlus Weather) files to calculate running mean outdoor temperature
  - Uses 30-day exponentially weighted running mean per ASHRAE 55 requirements
  - Extracts hourly radiation data for accurate MRT calculations

## Model Selection

ASHRAE55Foam supports two modes:

### 1. `adaptive`

Use this mode for naturally ventilated spaces. The tool evaluates ASHRAE 55 adaptive comfort with:

- `TOp`
- `ASHRAELevel80`
- `ASHRAELevel90`

This mode uses `-runningMean` or `-epw` to obtain the running mean outdoor temperature.

### 2. `pmv`

Use this mode for indoor conditioned spaces when you want an ASHRAE-55-style PMV check in this tool. The tool evaluates:

- `TOp`
- `PMV`
- `PPD`
- `ASHRAEAcceptable`

This mode does not require EPW data unless you also want EPW-based solar/MRT handling.

## Installation

### Prerequisites

- OpenFOAM (v2412 or later, optimized for v2412+)
- Standard OpenFOAM development environment

### Build Instructions

```bash
# Navigate to the ASHRAE55 directory
cd $FOAM_RUN/HVAC-for-OpenFOAM/ASHRAE55

# Clean previous builds (optional)
wclean

# Compile the utility
wmake
```

## Usage

### Basic Usage

```bash
ASHRAE55Foam [OPTIONS]
```

### Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `-epw <fileName>` | EPW weather file for running mean calculation and location data | - |
| `-dayOfYear <scalar>` | Day of year (1-365) for EPW calculation | 180 |
| `-hour <scalar>` | Hour of day for solar calculations (0-24) | 12.0 |
| `-solarData` | Use detailed solar calculations from EPW (only if no OpenFOAM solar model) | false |
| `-latitude <scalar>` | Site latitude in degrees (overrides EPW value) | EPW value |
| `-longitude <scalar>` | Site longitude in degrees (overrides EPW value) | EPW value |
| `-runningMean <scalar>` | Directly specify running mean outdoor temperature in C | - |
| `-model <word>` | ASHRAE 55 model: `adaptive` or `pmv` | `adaptive` |
| `-met <scalar>` | Metabolic rate in met for PMV mode | `1.2` |
| `-clo <scalar>` | Clothing insulation in clo for PMV mode | `0.5` |
| `-wme <scalar>` | External work in met for PMV mode | `0.0` |
| `-rh <scalar>` | Relative humidity in % for PMV mode if no field is available | `50` |

### Examples

1. **Using EPW file with automatic location detection:**
   ```bash
   # Automatically uses location from EPW file
   ASHRAE55Foam -epw weather.epw -dayOfYear 200
   ```

2. **Bangkok example with solar calculations:**
   ```bash
   # Automatically detects Bangkok coordinates from EPW
   ASHRAE55Foam -epw THA_Bangkok.484560_IWEC.epw -dayOfYear 180 -hour 14 -solarData
   ```

3. **Override EPW coordinates if needed:**
   ```bash
   ASHRAE55Foam -epw weather.epw -solarData -latitude 40.7 -longitude -74.0 -hour 15
   ```

4. **Direct specification of running mean temperature:**
   ```bash
   ASHRAE55Foam -runningMean 25.5
   ```

5. **Adaptive mode with explicit model selection:**
   ```bash
   ASHRAE55Foam -model adaptive -runningMean 25.5
   ```

6. **Indoor PMV mode for conditioned spaces:**
   ```bash
   ASHRAE55Foam -model pmv -met 1.2 -clo 0.7 -rh 50
   ```

7. **Indoor PMV mode using defaults from `constant/ASHRAE55Dict`:**
   ```bash
   ASHRAE55Foam -model pmv
   ```

## Required Input Fields

The following fields must be present in your OpenFOAM case:

- **T** (Temperature) - MUST_READ
- **U** (Velocity) - MUST_READ

Optional fields:

- **G** (Radiation) - READ_IF_PRESENT (automatically used if available)
- **thermo:relHum**, **thermoRelHum**, **relHum**, **RH**, or **relativeHumidity** - optional humidity sources for `pmv` mode

## Output Fields

The utility creates these fields depending on the selected model:

1. **TOp** - Operative temperature [K]
2. **ASHRAELevel80** - Binary field (0/1) indicating 80% adaptive acceptability compliance
3. **ASHRAELevel90** - Binary field (0/1) indicating 90% adaptive acceptability compliance
4. **PMV** - Predicted Mean Vote for PMV mode
5. **PPD** - Predicted Percentage of Dissatisfied for PMV mode
6. **ASHRAEAcceptable** - Binary field (0/1) for PMV acceptability (`-0.5 <= PMV <= 0.5`)

## Technical Details

### Adaptive Comfort Model

The ASHRAE 55-2020 adaptive comfort model is implemented as:

- Neutral temperature: `t_cmf = 0.31 × T_rm_out + 17.8` (C)
- 80% acceptability: ±3.5 C from neutral temperature
- 90% acceptability: ±2.5 C from neutral temperature
- Valid for naturally ventilated spaces with running mean outdoor temperature between 10 C and 33.5 C

### PMV Mode

For mechanically conditioned indoor spaces, `ASHRAE55Foam` can run in `pmv` mode. In this mode it:

- Computes PMV and PPD from air temperature, MRT, air speed, humidity, clothing, and metabolic rate
- Reads humidity from common OpenFOAM relative-humidity field names when available
- Falls back to `RH` from `ASHRAE55Dict` or `-rh` if no humidity field exists
- Marks cells as acceptable when `PMV` is between `-0.5` and `+0.5`

Default PMV inputs:

- `met = 1.2`
- `clo = 0.5`
- `wme = 0.0`
- `RH = 50`

These defaults can be set in `constant/ASHRAE55Dict` or overridden on the command line.

### Cooling Effect from Air Speed

When operative temperature > 25 C and air speed ≥ 0.6 m/s:
- 0.6-0.9 m/s: 1.2K cooling effect
- 0.9-1.2 m/s: 1.8K cooling effect
- >1.2 m/s: 2.2K cooling effect

### Mean Radiant Temperature Calculation

- **With OpenFOAM radiation model**: Converts irradiance G to MRT using empirical correlations
- **With EPW solar data**: Calculates solar position based on location and time zone
- **Without radiation data**: Uses area-weighted wall temperature average
- Solar effects are limited to realistic ranges (5-20 C temperature rise)
- Time zone corrections ensure accurate solar calculations globally

### Location Data from EPW

The tool automatically extracts from EPW files:
- Latitude and longitude for solar position calculations
- Time zone (GMT offset) for accurate local solar time
- Location name for reference
- Hourly weather and radiation data

## Applications

- HVAC system design and evaluation
- Building thermal comfort assessment
- Urban microclimate studies
- Natural ventilation design
- Compliance checking with ASHRAE 55-2020 standards

## Notes

- The tool is designed to work globally with any standard EPW weather file
- MRT formulas are universally valid but not specifically optimized for tropical regions
- For most accurate results, use OpenFOAM v2412+ with native solar radiation models
- `constant/ASHRAE55Dict` can be used to set default `model`, `met`, `clo`, `wme`, and `RH`
- `adaptive` and `pmv` are separate ASHRAE 55 evaluation paths in this tool; `pmv` does not depend on `comfortFoam`

## Version Information

- **Version**: 1.3
- **Author**: Thomas Tian
- **License**: GPL-3.0
- **Updates**: Added automatic global location support from EPW files

## References

1. ANSI/ASHRAE Standard 55-2020: Thermal Environmental Conditions for Human Occupancy
2. de Dear, R. J., & Brager, G. S. (1998). Developing an adaptive model of thermal comfort and preference

## Support

For issues, questions, or contributions, please visit the [GitHub repository](https://github.com/your-repo/HVAC-for-OpenFOAM).
