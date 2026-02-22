# ISO7243Foam

## Overview

`ISO7243Foam` is an OpenFOAM post-processing utility for ISO 7243 style heat stress assessment using WBGT (Wet Bulb Globe Temperature).

It is designed for workflows where CFD is run first (for example with `buoyantHumiditySimpleFoam`) and heat stress metrics are computed afterward.

## Key Features

- Calculates WBGT for both ISO 7243 formulas:
  - Indoor/no solar: `WBGT = 0.7*Tnwb + 0.3*Tg`
  - Outdoor/solar: `WBGT = 0.7*Tnwb + 0.2*Tg + 0.1*Ta`
- Reads humidity as post-processing input with priority for `buoyantHumiditySimpleFoam` outputs:
  - `thermo:relHum`, `thermoRelHum`, `relHum`, `RH`, `relativeHumidity`
- Supports direct input fields for `Tnwb` and `Tg`
- Provides practical fallbacks when `Tnwb` or `Tg` are missing
- Supports `-cellSet` and `-cellZone` analysis regions
- Full time-selection support (`-latestTime`, `-time`, etc.)

## Build

```bash
cd $FOAM_RUN/../HVAC-for-OpenFOAM/ISO7243
wmake
```

## Usage

```bash
ISO7243Foam
```

### Formula Selection

- Use indoor/no-solar formula:

```bash
ISO7243Foam -indoorNoSolar
```

- Use outdoor/solar formula:

```bash
ISO7243Foam -outdoorSolar
```

### Region-Restricted Analysis

```bash
ISO7243Foam -cellZone occupiedZone
ISO7243Foam -cellSet workerVolume
```

## Configuration

Create `constant/ISO7243Dict` in your case (template provided in this folder).

Main options:

- `RH`: fallback relative humidity [%] when no humidity field is present
- `useOutdoorSolarFormula`: default formula for `WBGTISO7243_C`
- `tNwbFieldName`, `tgFieldName`, `tmrtFieldName`: custom field names
- `globeTemperatureOffset`: optional correction for estimated `Tg`
- `globeEmissivity`: used if `Tg` is estimated from `qr`
- `lowRiskLimit`, `moderateRiskLimit`, `highRiskLimit`: thresholds for `ISO7243RiskLevel`

## Output Fields

All output values are in degC and stored as dimensionless scalar fields:

- `TnwbISO7243_C`
- `TgISO7243_C`
- `WBGTIndoorISO7243_C`
- `WBGTOutdoorISO7243_C`
- `WBGTISO7243_C`
- `ISO7243RiskLevel`

## Notes

- For strict ISO 7243 compliance, provide physically based `Tnwb` and `Tg` fields.
- If `Tnwb` is missing, the tool estimates it from `Ta` and RH using the Stull approximation.
- If `Tg` is missing, estimation priority is:
  `Tmrt` -> `G` -> `IDefault` -> `qr` -> area-weighted wall temperature.
- `ISO7243RiskLevel` uses user-defined thresholds and is not a direct replacement for full normative workload/acclimatization tables.
