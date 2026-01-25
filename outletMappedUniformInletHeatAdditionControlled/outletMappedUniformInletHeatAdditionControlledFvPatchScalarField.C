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

\*---------------------------------------------------------------------------*/

#include "outletMappedUniformInletHeatAdditionControlledFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "basicThermo.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::outletMappedUniformInletHeatAdditionControlledFvPatchScalarField::
outletMappedUniformInletHeatAdditionControlledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    outletPatchName_(),
    phiName_("phi"),
    TMin_(0),
    TMax_(5000),
    sensorPoint_(point(0, 0, 0)),
    Tset_(0),
    band_(0),
    controlMode_("scale"),
    gain_(0),
    Qbase_(0),
    Qmin_(0),
    Qmax_(0),
    rampRate_(0),
    Qeff_(0),
    Qprev_(0),
    timeIndex_(db().time().timeIndex()),
    debug_(false),
    debugTimeIndex_(-1),
    warnedNoCell_(false)
{}


Foam::outletMappedUniformInletHeatAdditionControlledFvPatchScalarField::
outletMappedUniformInletHeatAdditionControlledFvPatchScalarField
(
    const outletMappedUniformInletHeatAdditionControlledFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    outletPatchName_(ptf.outletPatchName_),
    phiName_(ptf.phiName_),
    TMin_(ptf.TMin_),
    TMax_(ptf.TMax_),
    sensorPoint_(ptf.sensorPoint_),
    Tset_(ptf.Tset_),
    band_(ptf.band_),
    controlMode_(ptf.controlMode_),
    gain_(ptf.gain_),
    Qbase_(ptf.Qbase_),
    Qmin_(ptf.Qmin_),
    Qmax_(ptf.Qmax_),
    rampRate_(ptf.rampRate_),
    Qeff_(ptf.Qeff_),
    Qprev_(ptf.Qprev_),
    timeIndex_(ptf.timeIndex_),
    debug_(ptf.debug_),
    debugTimeIndex_(ptf.debugTimeIndex_),
    warnedNoCell_(ptf.warnedNoCell_)
{}


Foam::outletMappedUniformInletHeatAdditionControlledFvPatchScalarField::
outletMappedUniformInletHeatAdditionControlledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    outletPatchName_(dict.get<word>("outletPatch")),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    TMin_(dict.getOrDefault<scalar>("TMin", 0)),
    TMax_(dict.getOrDefault<scalar>("TMax", 5000)),
    sensorPoint_(dict.get<point>("sensorPoint")),
    Tset_(dict.get<scalar>("Tset")),
    band_(dict.getOrDefault<scalar>("band", 0)),
    controlMode_(dict.getOrDefault<word>("controlMode", "scale")),
    gain_(dict.getOrDefault<scalar>("gain", 0)),
    Qbase_(dict.get<scalar>("Qbase")),
    Qmin_(dict.getOrDefault<scalar>("Qmin", 0)),
    Qmax_(dict.getOrDefault<scalar>("Qmax", Qbase_)),
    rampRate_(dict.getOrDefault<scalar>("rampRate", 0)),
    Qeff_(dict.getOrDefault<scalar>("Qeff", Qbase_)),
    Qprev_(dict.getOrDefault<scalar>("Qprev", Qeff_)),
    timeIndex_(db().time().timeIndex()),
    debug_(dict.getOrDefault<Switch>("debug", false)),
    debugTimeIndex_(-1),
    warnedNoCell_(false)
{}


Foam::outletMappedUniformInletHeatAdditionControlledFvPatchScalarField::
outletMappedUniformInletHeatAdditionControlledFvPatchScalarField
(
    const outletMappedUniformInletHeatAdditionControlledFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    outletPatchName_(ptf.outletPatchName_),
    phiName_(ptf.phiName_),
    TMin_(ptf.TMin_),
    TMax_(ptf.TMax_),
    sensorPoint_(ptf.sensorPoint_),
    Tset_(ptf.Tset_),
    band_(ptf.band_),
    controlMode_(ptf.controlMode_),
    gain_(ptf.gain_),
    Qbase_(ptf.Qbase_),
    Qmin_(ptf.Qmin_),
    Qmax_(ptf.Qmax_),
    rampRate_(ptf.rampRate_),
    Qeff_(ptf.Qeff_),
    Qprev_(ptf.Qprev_),
    timeIndex_(ptf.timeIndex_),
    debug_(ptf.debug_),
    debugTimeIndex_(ptf.debugTimeIndex_),
    warnedNoCell_(ptf.warnedNoCell_)
{}


Foam::outletMappedUniformInletHeatAdditionControlledFvPatchScalarField::
outletMappedUniformInletHeatAdditionControlledFvPatchScalarField
(
    const outletMappedUniformInletHeatAdditionControlledFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    outletPatchName_(ptf.outletPatchName_),
    phiName_(ptf.phiName_),
    TMin_(ptf.TMin_),
    TMax_(ptf.TMax_),
    sensorPoint_(ptf.sensorPoint_),
    Tset_(ptf.Tset_),
    band_(ptf.band_),
    controlMode_(ptf.controlMode_),
    gain_(ptf.gain_),
    Qbase_(ptf.Qbase_),
    Qmin_(ptf.Qmin_),
    Qmax_(ptf.Qmax_),
    rampRate_(ptf.rampRate_),
    Qeff_(ptf.Qeff_),
    Qprev_(ptf.Qprev_),
    timeIndex_(ptf.timeIndex_),
    debug_(ptf.debug_),
    debugTimeIndex_(ptf.debugTimeIndex_),
    warnedNoCell_(ptf.warnedNoCell_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::outletMappedUniformInletHeatAdditionControlledFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const Time& time = db().time();

    if (timeIndex_ != time.timeIndex())
    {
        timeIndex_ = time.timeIndex();
        Qprev_ = Qeff_;
    }

    const volScalarField& T =
    (
        dynamic_cast<const volScalarField&>(this->internalField())
    );

    // Sensor evaluation
    const label cellId = mesh.findCell(sensorPoint_);
    scalar TsensLocal = -GREAT;
    if (cellId != -1)
    {
        TsensLocal = T[cellId];
    }

    const scalar Tsens = returnReduce(TsensLocal, maxOp<scalar>());
    const bool sensorFound = (Tsens > -0.5*GREAT);

    // Control law for effective heat addition
    scalar Qeff = Qbase_;

    if (sensorFound)
    {
        const scalar Tband = Tset_ + band_;
        const bool linearMode = (controlMode_ == "linear");

        if (Tsens > Tband)
        {
            const scalar err = Tsens - Tband;
            if (linearMode)
            {
                Qeff = clamp(Qbase_ - gain_*err, Qmin_, Qmax_);
            }
            else
            {
                const scalar factor = max(scalar(0), 1 - gain_*err);
                Qeff = clamp(Qbase_*factor, Qmin_, Qmax_);
            }
        }
        else
        {
            Qeff = Qbase_;
        }

        if (rampRate_ > 0)
        {
            const scalar dt = time.deltaTValue();
            if (dt > SMALL)
            {
                const scalar dQ = rampRate_*dt;
                Qeff = clamp(Qeff, Qprev_ - dQ, Qprev_ + dQ);
            }
        }
    }
    else
    {
        Qeff = Qbase_;
        if (!warnedNoCell_)
        {
            WarningInFunction
                << "sensorPoint " << sensorPoint_
                << " not found in mesh. Using Qbase." << endl;
            warnedNoCell_ = true;
        }
    }

    Qeff_ = Qeff;

    const fvPatch& fvp = patch();

    label outletPatchID =
        fvp.patch().boundaryMesh().findPatchID(outletPatchName_);

    if (outletPatchID < 0)
    {
        FatalErrorInFunction
            << "Unable to find outlet patch " << outletPatchName_
            << abort(FatalError);
    }

    const fvPatch& outletPatch = fvp.boundaryMesh()[outletPatchID];

    const fvPatchField<scalar>& outletPatchField =
        T.boundaryField()[outletPatchID];

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const scalarField& outletPatchPhi = phi.boundaryField()[outletPatchID];
    const scalar sumOutletPatchPhi = gSum(outletPatchPhi);

    if (sumOutletPatchPhi > SMALL)
    {
        scalar averageOutletField =
            gWeightedSum(outletPatchPhi, outletPatchField)/sumOutletPatchPhi;

        const basicThermo& thermo =
            db().lookupObject<basicThermo>(basicThermo::dictName);

        const scalarField& pp = thermo.p().boundaryField()[outletPatchID];
        const scalarField& pT = thermo.T().boundaryField()[outletPatchID];

        const scalarField Cpf(thermo.Cp(pp, pT, outletPatchID));

        scalar totalPhiCp = sumOutletPatchPhi*gAverage(Cpf);

        operator==(clamp(averageOutletField + Qeff_/totalPhiCp, TMin_, TMax_));
    }
    else
    {
        scalar averageOutletField =
            gWeightedAverage(outletPatch.magSf(), outletPatchField);

        operator==(averageOutletField);
    }

    if (debug_)
    {
        const label ti = time.timeIndex();
        if (debugTimeIndex_ != ti)
        {
            debugTimeIndex_ = ti;
            Info<< typeName << " " << patch().name()
                << " Tsens=" << Tsens
                << " Qeff=" << Qeff_
                << " Qbase=" << Qbase_
                << " Tset=" << Tset_
                << " band=" << band_
                << " gain=" << gain_
                << " mode=" << controlMode_
                << " rampRate=" << rampRate_
                << endl;
        }
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::outletMappedUniformInletHeatAdditionControlledFvPatchScalarField::
write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);

    os.writeEntry("outletPatch", outletPatchName_);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);

    os.writeEntry("Qbase", Qbase_);
    os.writeEntry("Qmin", Qmin_);
    os.writeEntry("Qmax", Qmax_);

    os.writeEntry("sensorPoint", sensorPoint_);
    os.writeEntry("Tset", Tset_);
    os.writeEntry("band", band_);
    if (controlMode_ != "scale")
    {
        os.writeEntry("controlMode", controlMode_);
    }
    os.writeEntry("gain", gain_);
    os.writeEntry("rampRate", rampRate_);

    os.writeEntry("TMin", TMin_);
    os.writeEntry("TMax", TMax_);

    if (debug_)
    {
        os.writeEntry("debug", debug_);
    }

    os.writeEntry("Qeff", Qeff_);
    os.writeEntry("Qprev", Qprev_);

    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        outletMappedUniformInletHeatAdditionControlledFvPatchScalarField
    );
}


// ************************************************************************* //
