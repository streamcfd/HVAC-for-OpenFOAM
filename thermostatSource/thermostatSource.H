/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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
    Foam::fv::thermostatSource

Group
    grpFvOptionsSources

Description
    Applies a thermostat-controlled scalar source (e.g. heating power) to
    energy equations (h or e). The source is switched on/off based on the
    volume-averaged temperature in a sensor region.

    Operates as a proportional controller (P-Regler), similar to a
    thermostatic radiator valve (Thermostatventil):
    - When T_sensor <= TLower: full heating power Q (100%)
    - When T_sensor >= TUpper: no heating (0%)
    - Between TLower and TUpper: linear interpolation (proportional band)

    \f[
        Q_{eff} = Q \cdot \frac{T_{upper} - T_{sensor}}{T_{upper} - T_{lower}}
    \f]

    The source term is applied as an explicit contribution Su [W] or [W/m3]
    to the selected cells (the "heater" region). The sensor region can be
    a different cellZone/cellSet.

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    thermostat1
    {
        // Mandatory entries
        type            thermostatSource;
        selectionMode   cellZone;
        cellZone        heaterZone;

        // Field to apply source to (h or e)
        fields          (h);

        // Thermostat settings
        volumeMode      absolute;       // absolute [W] or specific [W/m3]
        Q               1000;           // heating power when ON

        // Sensor configuration
        sensorSelectionMode   cellZone;  // all/cellZone/cellSet/points
        sensorCellZone        sensorZone;

        // Proportional band (like a thermostatic valve)
        TLower          290;            // 100% power below this [K]
        TUpper          300;            // 0% power above this [K]

        // Optional entries
        TName           T;              // temperature field name (default: T)
        TMax            318.15;         // max heater temperature (Vorlauf) [K]
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                          | Type   | Reqd | Dflt
      type      | Type name: thermostatSource          | word   | yes  | -
      fields    | Name of operand field (h or e)       | wordList| yes | -
      volumeMode| Volume mode: absolute or specific    | word   | yes  | -
      Q         | Heating power when ON                | scalar | yes  | -
      sensorSelectionMode | Cell selection for sensor  | word   | yes  | -
      sensorCellZone | Cell zone name for sensor       | word   | cndtnl | -
      sensorCellSet  | Cell set name for sensor        | word   | cndtnl | -
      sensorPoint    | Point location for sensor       | vector | cndtnl | -
      TLower    | 100% power below this temperature     | scalar | yes  | -
      TUpper    | 0% power above this temperature      | scalar | yes  | -
      TName     | Name of temperature field            | word   | no   | T
      TMax      | Max heater temperature (Vorlauf)     | scalar | no   | GREAT
    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink
     - \link cellSetOption.H \endlink

SourceFiles
    thermostatSource.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fv_thermostatSource_H
#define Foam_fv_thermostatSource_H

#include "cellSetOption.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class thermostatSource Declaration
\*---------------------------------------------------------------------------*/

class thermostatSource
:
    public fv::cellSetOption
{
public:

    // Public Enumerations

        //- Options for the volume mode type
        enum volumeModeType
        {
            vmAbsolute,
            vmSpecific
        };

        //- Names for volumeModeType
        static const Enum<volumeModeType> volumeModeTypeNames_;

        //- Options for the sensor selection mode
        enum sensorSelectionModeType
        {
            ssmAll,
            ssmCellZone,
            ssmCellSet,
            ssmPoint
        };

        //- Names for sensorSelectionModeType
        static const Enum<sensorSelectionModeType>
            sensorSelectionModeTypeNames_;


private:

    // Private Data

        //- Volume mode
        volumeModeType volumeMode_;

        //- Heating power [W] or [W/m3]
        scalar Q_;

        //- Lower temperature threshold (100% power below this)
        scalar TLower_;

        //- Upper temperature threshold (0% power above this)
        scalar TUpper_;

        //- Name of temperature field
        word TName_;

        //- Max heater temperature (Vorlauftemperatur limit)
        scalar TMax_;

        //- Sensor cell selection mode
        sensorSelectionModeType sensorSelectionMode_;

        //- Sensor zone/set name
        word sensorName_;

        //- Sensor point (for point selection mode)
        point sensorPoint_;

        //- Sensor cell indices
        labelList sensorCells_;


    // Private Member Functions

        //- Set up the sensor cell selection
        void setSensorCells();

        //- Compute volume-averaged temperature in sensor region
        scalar sensorTemperature() const;

        //- Compute max temperature in heater cells
        scalar heaterMaxTemperature() const;

        //- No copy construct
        thermostatSource(const thermostatSource&) = delete;

        //- No copy assignment
        void operator=(const thermostatSource&) = delete;


public:

    //- Runtime type information
    TypeName("thermostatSource");


    // Constructors

        //- Construct from components
        thermostatSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~thermostatSource() = default;


    // Member Functions

    // Evaluation

        //- Add explicit contribution to incompressible equation
        virtual void addSup
        (
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to compressible equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );


    // IO

        //- Read source dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
