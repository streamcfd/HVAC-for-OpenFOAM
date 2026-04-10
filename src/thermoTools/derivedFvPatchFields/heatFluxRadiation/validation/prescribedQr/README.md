# heatFluxRadiation `qr` Validation

Minimal theory case for validating the `qr` term in `heatFluxRadiation`.

## Setup

- Solver: `buoyantSimpleFoam`
- Geometry: box `0.1 x 1 x 1 m`
- Patch `hot`: area `A = 1 m^2`
- Patch `cold`: `fixedValue 293.15 K`
- Initial state: everything at `293.15 K`
- Gravity: `0`
- Flow: laminar, `U = 0`
- Additional source: `Q constant 0`
- Radiative flux: field `qr`, prescribed as `1000 W/m^2` on `hot`

## Theory for the first iteration

At the start, `T_hot = T_cell = 293.15 K`, so `Qconv = 0` in the first
update. With `deltaT = 1`, `mass = 10 kg`, `Cp = 1000 J/kg/K`,
`A = 1 m^2`, the implemented update gives:

`Delta T = deltaT * A * qr / (mass * Cp)`

which becomes:

`Delta T = 1 * 1 * 1000 / (10 * 1000) = 0.1 K`

Expected patch temperature after one step:

`T_hot,1 = 293.15 + 0.1 = 293.25 K`

## Result

After one iteration the case writes:

- `refValue uniform 293.25;`

to `1/T` for patch `hot`.

With the same setup and reversed sign `qr = -1000 W/m^2`, the expected
first-step value is `293.05 K`.
