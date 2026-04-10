# heatFluxRadiation `fvDOM` Parallel-Plates Validation

System test for `heatFluxRadiation + fvDOM` with two black diffusive
parallel plates.

## Setup

- Solver: `buoyantSimpleFoam`
- Geometry: box `0.1 x 1 x 1 m`
- `y` and `z` sides: cyclic
- `hot`: `fixedValue 333.15 K`
- `cold`: `heatFluxRadiation`, `Q = 0`, initial `293.15 K`
- Radiation model: `fvDOM`
- Medium: `absorptionEmissionModel none`
- Surfaces: black (`absorptivity = emissivity = 1`)

## Theory

For infinite black parallel plates, the net radiative heat flux onto the cold
plate is:

`qr = sigma * (T_hot^4 - T_cold^4)`

With `T_hot = 333.15 K` and `T_cold = 293.15 K`, this gives a positive heat
flux onto the cold plate. In the first solver step, `Qconv = 0` because the
cold patch and its adjacent cell both start at `293.15 K`. The BC update for
the following step is therefore:

`Delta T = deltaT * A * qr / (mass * Cp)`

with `A = 1 m^2`, `deltaT = 1`, `mass = 10`, `Cp = 1000`.

Numerically:

- `qr_theory = 279.741473 W/m^2`
- `T_cold,2,theory = 293.177974 K`

## Result

In time directory `1`, `fvDOM` gives:

- `qr(cold) = 279.743 W/m^2`

in `1/qr`, which differs from theory by only about `0.0015 W/m^2`.

The `heatFluxRadiation` patch temperature in `1/T` still remains at the
initial value, because the radiative flux is only computed during that step.
In the next step the BC uses this `qr` directly:

- `T_cold,2 = 293.178 K`

in `2/T`, which matches the theoretical Euler step to within about
`2.6e-5 K`.
