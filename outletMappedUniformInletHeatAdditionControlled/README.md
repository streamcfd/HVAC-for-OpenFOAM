# outletMappedUniformInletHeatAdditionControlled

Runtime-selectable inlet temperature boundary condition for OpenFOAM (openfoam.com)
that maps outlet temperature and applies controlled heat addition based on a sensor
point in the volume.

This extends the behavior of `outletMappedUniformInletHeatAddition` with a simple
proportional controller, deadband, and optional ramp limiting.

## Build

Source your OpenFOAM environment first, then:

```bash
cd /home/flowsimpro/OpenFOAM/HVAC-for-OpenFOAM/outletMappedUniformInletHeatAdditionControlled
wmake libso
```

The library is built to:

```
$FOAM_USER_LIBBIN/liboutletMappedUniformInletHeatAdditionControlled.so
```

Add to `system/controlDict`:

```
libs ("liboutletMappedUniformInletHeatAdditionControlled.so");
```

## Usage (example in 0/T)

```foam
inlet
{
    type            outletMappedUniformInletHeatAdditionControlled;
    outletPatch     Outlet;
    phi             phi;

    Qbase           1000;        // base heat addition [W]
    Qmin            0;           // minimum allowed [W]
    Qmax            2000;        // maximum allowed [W]

    sensorPoint     (0.1 0.2 0.3);
    Tset            600;         // target temperature [K]
    band            5;           // deadband [K]
    gain            0.01;        // proportional gain
    rampRate        100;         // max |dQ/dt| [W/s]

    TMax            2000;        // safety cap (inherited)
    debug           1;           // optional debug output

    value           uniform 300;
}
```

## Dictionary keys

Required:
- `outletPatch` : name of outlet patch to map
- `Qbase`       : base heat addition
- `sensorPoint` : sensor location (x y z)
- `Tset`        : target temperature

Optional (defaults in parentheses):
- `phi`      (phi)
- `Qmin`     (0)
- `Qmax`     (Qbase)
- `band`     (0)
- `gain`     (0)
- `rampRate` (0)
- `TMin`     (0)
- `TMax`     (5000)
- `debug`    (0)

## Notes

- The sensor temperature is read from the cell containing `sensorPoint` and
  reduced across MPI with a global max (the owning rank wins).
- If the sensor point is outside the mesh, the BC uses `Qbase` and warns once.
