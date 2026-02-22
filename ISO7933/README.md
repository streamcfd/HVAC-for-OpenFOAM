# ISO7933Foam

## Overview

`ISO7933Foam` is an OpenFOAM post-processing utility for **dynamic PHS-oriented** heat strain assessment in ISO 7933 workflows.

It is intended for CFD workflows where environmental fields are solved first (e.g. with `buoyantHumiditySimpleFoam`) and heat-strain metrics are evaluated afterward.

## Scope Note

This implementation is a practical dynamic engineering model. It is **not a certified normative compliance calculator**.

Use it for hotspot detection, ranking, and scenario comparison. For formal compliance work, validate against dedicated reference procedures.

## Build

```bash
cd $FOAM_RUN/../HVAC-for-OpenFOAM/ISO7933
wmake
```

## Usage

```bash
ISO7933Foam
```

Common options:

- `-latestTime`
- `-time <value>`
- `-cellSet <name>`
- `-cellZone <name>`

## Humidity Input (Post-Processing)

Field priority:

1. `thermo:relHum` (from `buoyantHumiditySimpleFoam`)
2. `thermoRelHum`
3. `relHum`
4. `RH`
5. `relativeHumidity`
6. fallback to `RH` from `ISO7933Dict`

Values in `[0,1]` are automatically converted to `%`.

## Radiant Temperature Input

Field priority:

1. `Tr` (or `trFieldName`)
2. `Tmrt` (or `tmrtFieldName`)
3. `Tg` (or `tgFieldName`)
4. `G`
5. `IDefault`
6. `qr`
7. area-weighted wall temperature fallback

## Configuration

Create `constant/ISO7933Dict` in your case (template included in this folder).

Main parameter groups:

- Activity: `met`, `wme`, `metIsMetUnit`
- Clothing: `clo`, `im`
- Dynamics: `exposureDuration`, `timeStep`
- Physiology: `initialCoreTemperature`, `initialSkinTemperature`, `coreTemperatureLimit`
- Sweat model: `maxSweatRate`, `sweatControlGain`, `sweatTimeConstant`
- Dehydration limits: `lossLimit95Fraction`, `lossLimit50Fraction`

## Output Fields

- `TaISO7933_C`
- `TrISO7933_C`
- `PaISO7933_kPa`
- `ISO7933Ereq_Wm2`
- `ISO7933Emax_Wm2`
- `ISO7933Wreq`
- `ISO7933ReqSweatRate_gm2h`
- `ISO7933SweatRate_gm2h`
- `ISO7933Eact_Wm2`
- `ISO7933Tre_C`
- `ISO7933Tsk_C`
- `ISO7933PredWaterLoss_g`
- `ISO7933DlimTre_min`
- `ISO7933DlimLoss95_min`
- `ISO7933DlimLoss50_min`
- `ISO7933HeatStrainIndex`
- `ISO7933RiskLevel`

## Suggested Workflow

```bash
buoyantHumiditySimpleFoam
ISO7933Foam -latestTime
```
