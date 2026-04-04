# HVAC-for-OpenFOAM

[![OpenFOAM v2412+](https://img.shields.io/badge/OpenFOAM-v2412+-blue.svg)](https://www.openfoam.com/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Advanced HVAC simulation toolkit for OpenFOAM, providing specialized solvers and utilities for humidity modeling, thermal comfort analysis, and indoor climate assessment according to international standards (ISO 7730, ASHRAE-55, UTCI).

## Features

### Core Capabilities
- **Humidity Transport Modeling**: Full integration of humidity physics into OpenFOAM's thermophysical framework
- **Thermal Comfort Analysis**: Comprehensive comfort metrics (PMV, PPD, UTCI, ASHRAE-55)
- **Solar Radiation Integration**: Support for EPW weather data and OpenFOAM's native solar models
- **Buoyancy-Driven Flows**: Specialized solvers for natural and mixed convection with humidity effects
- **Building Physics**: Advanced boundary conditions for multi-layer walls with thermal mass
- **Wind-Driven Rain**: Lagrangian particle tracking for rain droplet impingement

### Components

#### CFD Solvers
- **buoyantHumidityPimpleFoam**: Transient solver for buoyant, turbulent flow with humidity transport
- **buoyantHumiditySimpleFoam**: Steady-state solver using SIMPLE algorithm
- **buoyantBoussinesqPimpleDyMFoam**: Dynamic mesh solver for Boussinesq approximation
- **windDrivenRainFoam**: Wind-driven rain simulation with multiple droplet phases

#### Comfort Analysis Tools
- **comfortFoam**: ISO 7730 thermal comfort (PMV, PPD, Draft Rating)
- **UTCIFoam**: Universal Thermal Climate Index calculator
- **ASHRAE55Foam**: ASHRAE Standard 55 compliance analysis
- **ISO7243Foam**: ISO 7243 WBGT heat stress assessment
- **ISO7933Foam**: ISO 7933 heat strain screening
- **AoAFoam**: Age of Air calculation for ventilation assessment

#### fvOptions
- **thermostatSource**: Thermostat-controlled heating source with proportional regulation, sensor placement, and supply temperature limiting (Vorlauf)

#### Boundary Conditions
- **buildingElementBC**: Advanced thermal boundary condition for building walls with multi-layer support

#### Libraries
- **humidityRhoThermo**: Thermophysical model extension for humidity calculations
- **solarCalculator**: Solar position and radiation calculations

## Documentation

### Solver Documentation
- [buoyantHumidityPimpleFoam](buoyantHumidityPimpleFoam/README.md) - Transient humidity solver
- [buoyantHumiditySimpleFoam](buoyantHumiditySimpleFoam/README.md) - Steady-state humidity solver
- [humidityRhoThermo](humidityRhoThermo/README.md) - Humidity thermophysical library
- [windDrivenRainFoam](windDrivenRainFoam/README.md) - Wind-driven rain simulation

### Comfort Analysis Tools
- [comfortFoam](comfortFoam/README.md) - ISO 7730 comfort metrics
- [UTCIFoam](UTCIFoam/README.md) - UTCI calculation with weather data support
- [ASHRAE55Foam](ASHRAE55/README.md) - ASHRAE-55 adaptive comfort model
- [ISO7243Foam](ISO7243/README.md) - ISO 7243 WBGT post-processing
- [ISO7933Foam](ISO7933/README.md) - ISO 7933 heat strain screening

### fvOptions
- [thermostatSource](thermostatSource/README.md) - Thermostat-controlled heating with sensor and Vorlauf limit

### Boundary Conditions
- [buildingElementBC](buildingElementBC/README.md) - Multi-layer wall boundary condition with radiation

## OpenFOAM-v2512 Patch Set

This repository contains a ready-to-copy patch tree for OpenFOAM v2512 at:

- `OpenFOAM-v2512/src/thermophysicalModels/radiation/`

The patch extends solar radiation handling for partially transparent surfaces (`0 < t < 1`) by:

- enabling transmissivity from `transparent` boundary radiation properties
- propagating ray transmissivity through `solarLoad` face shading
- scaling direct and reflected solar load contributions consistently

To apply in an OpenFOAM-v2512 source tree, copy the matching files from this folder into the OpenFOAM source and rebuild the radiation libraries with `wmake`.

## Technical Details

### Humidity Modeling
The humidity transport is integrated into OpenFOAM's thermophysical framework through:
- Extended equation of state with humidity effects
- Humidity-dependent thermophysical properties
- Custom boundary conditions (e.g., `fixedHumidity`)
- Fields: relative humidity, specific humidity, water vapor partial pressure

### Comfort Calculations
All comfort tools follow established standards:
- **ISO 7730**: Fanger's PMV/PPD model with local draft rating
- **ISO 7243**: WBGT-based heat stress assessment
- **ISO 7933**: Heat strain screening based on required evaporation
- **UTCI**: 6th-order polynomial approximation with radiation effects
- **ASHRAE-55**: Adaptive comfort model for naturally ventilated spaces

## Applications

- HVAC system design and optimization
- Indoor air quality assessment
- Natural ventilation design
- Building energy simulation with dynamic thermal mass
- Urban microclimate studies
- Thermal comfort evaluation
- Compliance checking with international standards
- Building facade performance analysis
- Rain penetration and moisture risk assessment

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.
