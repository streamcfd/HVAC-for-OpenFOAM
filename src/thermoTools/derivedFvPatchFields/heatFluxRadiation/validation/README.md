# heatFluxRadiation Validation Cases

This directory contains two focused validation cases for
`heatFluxRadiation`.

## Cases

- `prescribedQr`
  Validates the `qr` contribution of the boundary condition itself against an
  exact first-step energy balance with a prescribed `qr` field.

- `fvDOMParallelPlates`
  Validates the coupled behaviour of `fvDOM` and `heatFluxRadiation` for two
  black parallel plates. The radiative heat flux is checked in time `1`, and
  the corresponding patch-temperature update is checked in time `2`.
