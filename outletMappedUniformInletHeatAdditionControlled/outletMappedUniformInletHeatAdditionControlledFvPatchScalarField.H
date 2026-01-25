/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 OpenCFD Ltd.
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
    Foam::outletMappedUniformInletHeatAdditionControlledFvPatchScalarField

Group
    grpInletBoundaryConditions

Description
    This temperature boundary condition averages the temperature over the
    outlet patch specified by name "outletPatch" and applies an extra
    heat source. This is set as a uniform temperature value on this patch.
    Additionally TMin/TMax limits can be applied.

    A feedback controller is added to adjust the heat addition based on a
    sensor point inside the domain.
    The control mode can be "scale" (default) or "linear".
    "scale": Qeff = Qbase*factor with gain in 1/K.
    "linear": Qeff = Qbase - gain*err with gain in W/K (allows cooling).

    \heading Patch usage

    \table
        Property    | Description                    | Required | Default value
        outletPatch | name of outlet patch           | yes      |
        Qbase       | base heat addition             | yes      |
        Qmin        | minimum heat addition          | no       | 0
        Qmax        | maximum heat addition          | no       | Qbase
        sensorPoint | sensor point coordinate        | yes      |
        Tset        | target temperature             | yes      |
        band        | deadband                       | no       | 0
        controlMode | control mode (scale|linear)    | no       | scale
        gain        | proportional gain (mode dep.)  | no       | 0
        rampRate    | max |dQ/dt|                    | no       | 0
        phi         | flux field name                | no       | phi
        TMin        | min temperature limit          | no       | 0.0
        TMax        | max temperature limit          | no       | 5000.0
        debug       | print Tsens and Qeff           | no       | 0
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            outletMappedUniformInletHeatAdditionControlled;
        outletPatch     Outlet;
        phi             phi;
        Qbase           10;
        Qmin            0;
        Qmax            10;
        TMin            0;
        TMax            5000;
        sensorPoint     (0 0 0);
        Tset            500;
        band            5;
        gain            0.01;
        rampRate        1e3;
        debug           1;
        value           uniform 0;
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    outletMappedUniformInletHeatAdditionControlledFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef outletMappedUniformInletHeatAdditionControlledFvPatchScalarField_H
#define outletMappedUniformInletHeatAdditionControlledFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "fvPatchFields.H"
#include "Switch.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
  Class outletMappedUniformInletHeatAdditionControlledFvPatchScalarField
  Declaration
\*---------------------------------------------------------------------------*/

class outletMappedUniformInletHeatAdditionControlledFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of the outlet patch to be mapped
        word outletPatchName_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Minimum temperature limit
        scalar TMin_;

        //- Maximum temperature limit
        scalar TMax_;

        //- Sensor point coordinate
        point sensorPoint_;

        //- Target temperature
        scalar Tset_;

        //- Deadband
        scalar band_;

        //- Control mode ("scale" or "linear")
        word controlMode_;

        //- Proportional gain
        scalar gain_;

        //- Base heat addition
        scalar Qbase_;

        //- Minimum heat addition
        scalar Qmin_;

        //- Maximum heat addition
        scalar Qmax_;

        //- Ramp rate limit
        scalar rampRate_;

        //- Current effective heat addition
        scalar Qeff_;

        //- Previous time-step effective heat addition
        scalar Qprev_;

        //- Time index of last update
        label timeIndex_;

        //- Debug output switch
        Switch debug_;

        //- Last time index for debug output
        label debugTimeIndex_;

        //- Warn once if sensor point not found
        bool warnedNoCell_;


public:

    //- Runtime type information
    TypeName("outletMappedUniformInletHeatAdditionControlled");


    // Constructors

        //- Construct from patch and internal field
        outletMappedUniformInletHeatAdditionControlledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        outletMappedUniformInletHeatAdditionControlledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        outletMappedUniformInletHeatAdditionControlledFvPatchScalarField
        (
            const outletMappedUniformInletHeatAdditionControlledFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        outletMappedUniformInletHeatAdditionControlledFvPatchScalarField
        (
            const outletMappedUniformInletHeatAdditionControlledFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        outletMappedUniformInletHeatAdditionControlledFvPatchScalarField
        (
            const outletMappedUniformInletHeatAdditionControlledFvPatchScalarField&,
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

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
