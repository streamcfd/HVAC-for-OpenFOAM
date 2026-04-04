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

\*---------------------------------------------------------------------------*/

#include "thermostatSource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "cellSet.H"

// * * * * * * * * * * * * * Static Data Members  * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(thermostatSource, 0);
    addToRunTimeSelectionTable(option, thermostatSource, dictionary);
}
}

const Foam::Enum<Foam::fv::thermostatSource::volumeModeType>
Foam::fv::thermostatSource::volumeModeTypeNames_
({
    { volumeModeType::vmAbsolute, "absolute" },
    { volumeModeType::vmSpecific, "specific" },
});

const Foam::Enum<Foam::fv::thermostatSource::sensorSelectionModeType>
Foam::fv::thermostatSource::sensorSelectionModeTypeNames_
({
    { sensorSelectionModeType::ssmAll, "all" },
    { sensorSelectionModeType::ssmCellZone, "cellZone" },
    { sensorSelectionModeType::ssmCellSet, "cellSet" },
    { sensorSelectionModeType::ssmPoint, "point" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::thermostatSource::setSensorCells()
{
    switch (sensorSelectionMode_)
    {
        case ssmAll:
        {
            sensorCells_ = identity(mesh_.nCells());
            break;
        }
        case ssmCellZone:
        {
            const label zoneID = mesh_.cellZones().findZoneID(sensorName_);

            if (zoneID < 0)
            {
                FatalErrorInFunction
                    << "Sensor cellZone " << sensorName_
                    << " not found in mesh"
                    << exit(FatalError);
            }

            sensorCells_ = mesh_.cellZones()[zoneID];
            break;
        }
        case ssmCellSet:
        {
            cellSet cs(mesh_, sensorName_);
            sensorCells_ = cs.toc();
            break;
        }
        case ssmPoint:
        {
            const label cellI = mesh_.findCell(sensorPoint_);

            if (cellI < 0)
            {
                // Point not on this processor - empty list
                sensorCells_.clear();
            }
            else
            {
                sensorCells_.resize(1);
                sensorCells_[0] = cellI;
            }
            break;
        }
    }

    Info<< "    thermostatSource: sensor region has "
        << returnReduce(sensorCells_.size(), sumOp<label>())
        << " cells" << endl;
}


Foam::scalar Foam::fv::thermostatSource::sensorTemperature() const
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>(TName_);
    const scalarField& V = mesh_.V();

    scalar sumTV = 0;
    scalar sumV = 0;

    for (const label cellI : sensorCells_)
    {
        sumTV += T[cellI] * V[cellI];
        sumV += V[cellI];
    }

    reduce(sumTV, sumOp<scalar>());
    reduce(sumV, sumOp<scalar>());

    if (sumV > VSMALL)
    {
        return sumTV / sumV;
    }

    return 0;
}


Foam::scalar Foam::fv::thermostatSource::heaterMaxTemperature() const
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>(TName_);

    scalar maxT = -GREAT;

    for (const label cellI : cells_)
    {
        maxT = max(maxT, T[cellI]);
    }

    reduce(maxT, maxOp<scalar>());

    return maxT;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::thermostatSource::thermostatSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh),
    volumeMode_(vmAbsolute),
    Q_(0),
    TLower_(0),
    TUpper_(0),
    TName_("T"),
    TMax_(GREAT),
    sensorSelectionMode_(ssmAll),
    sensorName_(),
    sensorPoint_(Zero),
    sensorCells_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::thermostatSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    return this->addSup(volScalarField::null(), eqn, fieldi);
}


void Foam::fv::thermostatSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    // Read current sensor temperature
    const scalar Tsensor = sensorTemperature();

    // Proportional controller (like a thermostatic valve)
    //   T <= TLower: factor = 1 (100% power)
    //   T >= TUpper: factor = 0 (off)
    //   between:     linear interpolation
    scalar factor = 1;
    if (Tsensor >= TUpper_)
    {
        factor = 0;
    }
    else if (Tsensor > TLower_)
    {
        factor = (TUpper_ - Tsensor) / (TUpper_ - TLower_);
    }

    // Apply Vorlauf temperature limit on heater cells
    if (TMax_ < GREAT)
    {
        const scalar THeater = heaterMaxTemperature();
        const scalar band = 2.0;  // 2K proportional band

        if (THeater >= TMax_)
        {
            factor = 0;
        }
        else if (THeater > TMax_ - band)
        {
            factor *= (TMax_ - THeater) / band;
        }
    }

    const scalar Qeff = Q_ * factor;

    Info<< "thermostatSource " << name_
        << ": T_sensor=" << Tsensor << " K"
        << ", factor=" << factor*100 << "%"
        << ", Q_eff=" << Qeff
        << endl;

    if (factor < SMALL)
    {
        return;
    }

    // Apply source
    const dimensionSet SuDims(eqn.dimensions()/dimVolume);

    scalar VDash = 1;
    if (volumeMode_ == vmAbsolute)
    {
        VDash = V_;
    }

    const dimensioned<scalar> SuValue
    (
        "Su",
        SuDims,
        Qeff/VDash
    );

    if (this->useSubMesh())
    {
        auto tsu = DimensionedField<scalar, volMesh>::New
        (
            name_ + fieldNames_[fieldi] + "Su",
            IOobject::NO_REGISTER,
            mesh_,
            dimensioned<scalar>(SuDims, Zero)
        );
        UIndirectList<scalar>(tsu.ref(), cells_) = SuValue.value();
        eqn += tsu;
    }
    else
    {
        eqn += SuValue;
    }
}


bool Foam::fv::thermostatSource::read(const dictionary& dict)
{
    if (fv::cellSetOption::read(dict))
    {
        coeffs_.readEntry("fields", fieldNames_);
        resetApplied();

        volumeMode_ = volumeModeTypeNames_.get("volumeMode", coeffs_);
        coeffs_.readEntry("Q", Q_);
        coeffs_.readEntry("TLower", TLower_);
        coeffs_.readEntry("TUpper", TUpper_);
        coeffs_.readIfPresent("TName", TName_);
        TMax_ = coeffs_.getOrDefault<scalar>("TMax", GREAT);

        if (TLower_ > TUpper_)
        {
            FatalIOErrorInFunction(coeffs_)
                << "TLower (" << TLower_
                << ") must be <= TUpper (" << TUpper_ << ")"
                << exit(FatalIOError);
        }

        // Read sensor configuration
        sensorSelectionMode_ =
            sensorSelectionModeTypeNames_.get
            (
                "sensorSelectionMode",
                coeffs_
            );

        switch (sensorSelectionMode_)
        {
            case ssmCellZone:
            {
                coeffs_.readEntry("sensorCellZone", sensorName_);
                break;
            }
            case ssmCellSet:
            {
                coeffs_.readEntry("sensorCellSet", sensorName_);
                break;
            }
            case ssmPoint:
            {
                coeffs_.readEntry("sensorPoint", sensorPoint_);
                break;
            }
            case ssmAll:
            {
                break;
            }
        }

        setSensorCells();

        Info<< "    thermostatSource settings:" << nl
            << "        Q       = " << Q_ << nl
            << "        TLower  = " << TLower_ << nl
            << "        TUpper  = " << TUpper_ << nl
            << "        TName   = " << TName_ << nl
            << "        TMax    = "
            << (TMax_ < GREAT ? Foam::name(TMax_) : Foam::string("none"))
            << nl
            << "        sensor  = " << sensorSelectionModeTypeNames_
                                        [sensorSelectionMode_]
            << endl;

        return true;
    }

    return false;
}


// ************************************************************************* //
