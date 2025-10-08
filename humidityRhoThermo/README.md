# HVAC-for-OpenFOAM
HVACTools for OpenFOAM - Make OpenFOAM for HVACTool Simulation great again.

Humidity Libary:
•	Switch to the following folder :
cd $FOAM_SRC/thermophysicalModels/basic

•	Add the following lines after the liquidThermo.C in under Make folder in files:
humidityRhoThermo/humidityRhoThermo.C
humidityRhoThermo/humidityRhoThermos.C 
derivedFvPatchFields/fixedHumidity/fixedHumidityFvPatchScalarField.C

•	Compile the thermo library and boundary condition 
wmake libso

•	(Optional) Define humidity transport coefficients in
    `constant/thermophysicalProperties`, e.g.

```
humidityTransport
{
    molecularDiffusivity   [0 2 -1 0 0 0 0] 2.5e-05;
    turbulentSchmidtNumber 0.7;
}
```

  If omitted, the solver falls back to the defaults shown above.
