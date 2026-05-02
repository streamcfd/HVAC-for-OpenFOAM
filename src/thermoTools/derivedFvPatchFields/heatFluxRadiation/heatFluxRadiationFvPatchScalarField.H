/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | www.openfoam.com
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::heatFluxRadiationFvPatchScalarField

Group
    grpThermoBoundaryConditions grpWallBoundaryConditions

Description
    Lumped wall-temperature model with prescribed total power input and/or
    prescribed heat flux, plus optional radiative heat-flux coupling via
    \c qr and optional ambient coupling through \c h/\c Ta.

    The patch uses a single temperature value for the whole wall and updates it
    from the net power balance

        m Cp dT/dt = Q + integral(q dA) + Qrad + Qamb - Qconv

    where \c Q is a prescribed total power [W], \c q is an optional prescribed
    heat flux [W/m^2], \c Qrad is integrated from the optional \c qr field,
    \c Qamb is an optional ambient heat exchange term based on \c h and \c Ta,
    and \c Qconv is the conductive/convective heat loss to the adjacent fluid
    obtained from the local temperature gradient.

Usage
    \table
        Property       | Description                         | Required | Default
        mode           | Input mode: power, flux, coefficient | no      | legacy
        kappaMethod    | thermal conductivity option         | yes      |
        Q              | Total power input [W], required for mode power | conditional | none
        q              | Heat flux input [W/m^2], required for mode flux | conditional | none
        h              | Ambient heat transfer coefficient [W/m^2/K], required for mode coefficient | conditional | none
        Ta             | Ambient temperature [K], required for mode coefficient | conditional | none
        emissivity     | Emissivity for radiative exchange to ambient | no | 0
        thicknessLayers | Optional layer thicknesses [m]     | no       | none
        kappaLayers    | Optional layer conductivities [W/m/K] | no     | none
        Cp             | Heat capacity [J/kg/K]              | yes      |
        massPerArea    | Areal mass [kg/m^2]                 | yes      |
        qr             | Radiative flux field name           | no       | none
        qrRelaxation   | Relaxation for the radiative field  | no       | 1
        relaxation     | Relaxation for the temperature step | no       | 1
        minTemperature | Optional lower patch-temperature bound [K] | no | none
        maxTemperature | Optional upper patch-temperature bound [K] | no | none
    \endtable

    Example of the boundary condition specification:
    \verbatim
    person
    {
        type            heatFluxRadiation;
        mode            power;
        kappaMethod     fluidThermo;
        Q               constant 60;
        massPerArea     50;
        Cp              3500;
        qr              qr;
        qrRelaxation    0.7;
        minTemperature  280;
        maxTemperature  320;
        value           uniform 304.15;
    }
    \endverbatim

    A cooling floor can be specified directly with \c q:
    \verbatim
    floor
    {
        type            heatFluxRadiation;
        mode            flux;
        kappaMethod     fluidThermo;
        q               uniform -20;
        massPerArea     720;
        Cp              900;
        qr              qr;
        value           uniform 295.15;
    }
    \endverbatim

    Ambient coupling without prescribed power or heat flux can be specified with
    \c mode coefficient:
    \verbatim
    facade
    {
        type            heatFluxRadiation;
        mode            coefficient;
        kappaMethod     fluidThermo;
        h               constant 5;
        Ta              constant 293.15;
        massPerArea     720;
        Cp              900;
        value           uniform 293.15;
    }
    \endverbatim

SourceFiles
    heatFluxRadiationFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef heatFluxRadiationFvPatchScalarField_H
#define heatFluxRadiationFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"
#include "Function1.H"
#include "PatchFunction1.H"

namespace Foam
{

class heatFluxRadiationFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase
{
    // Private data

        //- Input mode: legacy, power, flux, coefficient
        word mode_;

        //- Prescribed total power [W]
        autoPtr<Function1<scalar>> Q_;

        //- Prescribed heat flux [W/m^2]
        autoPtr<PatchFunction1<scalar>> q_;

        //- Optional ambient heat transfer coefficient [W/m^2/K]
        autoPtr<PatchFunction1<scalar>> h_;

        //- Optional ambient temperature [K]
        autoPtr<Function1<scalar>> Ta_;

        //- Specific heat capacity [J/kg/K]
        scalar Cp_;

        //- Areal mass [kg/m^2]
        scalar massPerArea_;

        //- Optional surface emissivity for radiative transfer to ambient
        scalar emissivity_;

        //- Relaxation for the patch temperature update
        scalar relaxation_;

        //- Cache qr for relaxation
        scalarField qrPrevious_;

        //- Relaxation for qr
        scalar qrRelaxation_;

        //- Name of the radiative heat flux
        const word qrName_;

        //- Optional lower patch-temperature bound
        scalar minTemperature_;

        //- Optional upper patch-temperature bound
        scalar maxTemperature_;

        //- Thickness of optional external layers
        scalarList thicknessLayers_;

        //- Conductivity of optional external layers
        scalarList kappaLayers_;

        //- Current time index (used for updating)
        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("heatFluxRadiation");


    // Constructors

        //- Construct from patch and internal field
        heatFluxRadiationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        heatFluxRadiationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        heatFluxRadiationFvPatchScalarField
        (
            const heatFluxRadiationFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        heatFluxRadiationFvPatchScalarField
        (
            const heatFluxRadiationFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        heatFluxRadiationFvPatchScalarField
        (
            const heatFluxRadiationFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return fvPatchField<scalar>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap(const fvPatchField<scalar>&, const labelList&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            void write(Ostream&) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
