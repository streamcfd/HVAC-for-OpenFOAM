# thermostatSource

An `fvOption` source for OpenFOAM that provides thermostat-controlled heating with proportional regulation and optional supply temperature limiting.

## Description

Applies a thermostat-controlled scalar source (heating power) to energy equations (`h` or `e`). The source is regulated based on the temperature at a sensor location, which can be independent from the heater zone.

Operates as a proportional controller (P-Regler), similar to a thermostatic radiator valve:

- **T_sensor <= TLower**: full heating power Q (100%)
- **T_sensor >= TUpper**: no heating (0%)
- **Between TLower and TUpper**: linear interpolation (proportional band)

An optional `TMax` parameter limits the maximum heater zone temperature (Vorlauftemperatur), preventing unrealistic temperatures e.g. when modelling heat pump systems with low supply temperatures.

## Build

```bash
cd thermostatSource
wmake
```

## Usage

Add to `system/controlDict`:

```
libs (libthermostatSource);
```

Add to `constant/fvOptions`:

```
thermostat1
{
    type            thermostatSource;
    selectionMode   cellZone;
    cellZone        heaterZone;

    fields          (h);

    volumeMode      absolute;       // absolute [W] or specific [W/m3]
    Q               1000;           // heating power [W]

    // Sensor configuration
    sensorSelectionMode   point;    // all / cellZone / cellSet / point
    sensorPoint           (2.5 1.5 1.2);

    // Proportional band
    TLower          286.15;         // 13 C - full power below this
    TUpper          288.15;         // 15 C - off above this

    // Optional
    TMax            318.15;         // 45 C - max heater temperature (Vorlauf)
    TName           T;              // temperature field name (default: T)
}
```

## Parameters

| Parameter | Description | Type | Required | Default |
|---|---|---|---|---|
| type | `thermostatSource` | word | yes | - |
| fields | Operand field (`h` or `e`) | wordList | yes | - |
| volumeMode | `absolute` [W] or `specific` [W/m3] | word | yes | - |
| Q | Heating power | scalar | yes | - |
| sensorSelectionMode | `all` / `cellZone` / `cellSet` / `point` | word | yes | - |
| sensorCellZone | Cell zone name for sensor | word | conditional | - |
| sensorCellSet | Cell set name for sensor | word | conditional | - |
| sensorPoint | Point location for sensor | vector | conditional | - |
| TLower | 100% power below this temperature [K] | scalar | yes | - |
| TUpper | 0% power above this temperature [K] | scalar | yes | - |
| TMax | Max heater temperature (Vorlauf limit) [K] | scalar | no | none |
| TName | Temperature field name | word | no | T |

## Runtime Output

At each time step, the source prints:

```
thermostatSource thermostat1: T_sensor=287.3 K, factor=43.5%, Q_eff=435
```

## Compatible Solvers

Any compressible energy solver, e.g.:
- `buoyantSimpleFoam`
- `buoyantPimpleFoam`
- `buoyantHumiditySimpleFoam`
- `buoyantHumidityPimpleFoam`

## License

GPLv3 - see [LICENSE](../LICENSE)
