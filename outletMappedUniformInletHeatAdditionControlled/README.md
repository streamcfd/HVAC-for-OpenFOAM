# outletMappedUniformInletHeatAdditionControlled

Runtime-selectable inlet temperature boundary condition for OpenFOAM (openfoam.com)
that maps outlet temperature and applies controlled heat addition based on a sensor
point in the volume.

This extends the behavior of `outletMappedUniformInletHeatAddition` with a
controller, optional deadband, and optional ramp limiting. A linear control mode
allows cooling (negative heat addition).

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

### Heating (default scale mode)

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
    Tset            293.15;      // 20 C
    band            0.5;         // K
    gain            0.01;        // 1/K (scale mode)
    rampRate        0;           // W/s (0 = no limit)

    TMax            2000;        // safety cap (inherited)
    debug           1;           // optional debug output

    value           uniform 293.15;
}
```

### Cooling (linear mode)

```foam
inlet
{
    type            outletMappedUniformInletHeatAdditionControlled;
    outletPatch     Outlet;
    phi             phi;

    controlMode     linear;      // enables cooling behavior
    Qbase           0;
    Qmin            -60000;      // allow negative heat addition (cooling)
    Qmax            0;

    sensorPoint     (0.1 0.2 0.3);
    Tset            293.15;      // 20 C
    band            0;           // K
    gain            2000;        // W/K (linear mode)
    rampRate        0;           // W/s (0 = no limit)

    TMax            2000;
    value           uniform 293.15;
}
```

## Dictionary keys

Required:
- `outletPatch` : name of outlet patch to map
- `Qbase`       : base heat addition
- `sensorPoint` : sensor location (x y z)
- `Tset`        : target temperature

Optional (defaults in parentheses):
- `phi`         (phi)
- `Qmin`        (0)
- `Qmax`        (Qbase)
- `band`        (0)
- `gain`        (0)
- `controlMode` (scale)
- `rampRate`    (0)
- `TMin`        (0)
- `TMax`        (5000)
- `debug`       (0)

## Control law

If Tsens <= Tset + band:
- Qeff = Qbase

If Tsens > Tset + band:
- scale mode:  Qeff = clamp(Qbase * max(0, 1 - gain*err), Qmin, Qmax)
  gain units: 1/K
- linear mode: Qeff = clamp(Qbase - gain*err, Qmin, Qmax)
  gain units: W/K (allows cooling when Qmin < 0)

## Notes

- The sensor temperature is read from the cell containing `sensorPoint` and
  reduced across MPI with a global max (the owning rank wins).
- If the sensor point is outside the mesh, the BC uses `Qbase` and warns once.
- `rampRate` limits |dQ/dt| in W/s; set to 0 for immediate changes.
