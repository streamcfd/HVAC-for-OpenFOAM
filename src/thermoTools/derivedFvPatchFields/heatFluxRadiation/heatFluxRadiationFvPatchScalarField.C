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

\*---------------------------------------------------------------------------*/

#include "heatFluxRadiationFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

namespace Foam
{

heatFluxRadiationFvPatchScalarField::heatFluxRadiationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch()),
    Q_(nullptr),
    Cp_(0),
    mass_(0),
    relaxation_(1),
    qrPrevious_(),
    qrRelaxation_(1),
    qrName_("none"),
    minTemperature_(-VGREAT),
    maxTemperature_(VGREAT),
    curTimeIndex_(-1)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


heatFluxRadiationFvPatchScalarField::heatFluxRadiationFvPatchScalarField
(
    const heatFluxRadiationFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    temperatureCoupledBase(patch(), ptf),
    Q_(ptf.Q_.clone()),
    Cp_(ptf.Cp_),
    mass_(ptf.mass_),
    relaxation_(ptf.relaxation_),
    qrPrevious_(),
    qrRelaxation_(ptf.qrRelaxation_),
    qrName_(ptf.qrName_),
    minTemperature_(ptf.minTemperature_),
    maxTemperature_(ptf.maxTemperature_),
    curTimeIndex_(-1)
{
    if (qrName_ != "none")
    {
        qrPrevious_.resize(mapper.size());
        qrPrevious_.map(ptf.qrPrevious_, mapper);
    }
}


heatFluxRadiationFvPatchScalarField::heatFluxRadiationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    Q_(Function1<scalar>::New("Q", dict, &db())),
    Cp_(dict.get<scalar>("Cp")),
    mass_(dict.get<scalar>("mass")),
    relaxation_(dict.getOrDefault<scalar>("relaxation", 1)),
    qrPrevious_(),
    qrRelaxation_(dict.getOrDefault<scalar>("qrRelaxation", 1)),
    qrName_(dict.getOrDefault<word>("qr", "none")),
    minTemperature_(dict.getOrDefault<scalar>("minTemperature", -VGREAT)),
    maxTemperature_(dict.getOrDefault<scalar>("maxTemperature", VGREAT)),
    curTimeIndex_(-1)
{
    if (Cp_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "Cp must be > 0 for patch " << p.name()
            << exit(FatalIOError);
    }

    if (mass_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "mass must be > 0 for patch " << p.name()
            << exit(FatalIOError);
    }

    if (minTemperature_ > maxTemperature_)
    {
        FatalIOErrorInFunction(dict)
            << "minTemperature (" << minTemperature_ << ") exceeds "
            << "maxTemperature (" << maxTemperature_ << ')'
            << exit(FatalIOError);
    }

    if (qrName_ != "none")
    {
        if (dict.found("qrPrevious"))
        {
            qrPrevious_ = scalarField("qrPrevious", dict, p.size());
        }
        else
        {
            qrPrevious_.resize(p.size(), Zero);
        }
    }

    fvPatchFieldBase::readDict(dict);
    this->readValueEntry(dict, IOobjectOption::MUST_READ);
    refValue() = *this;
    refGrad() = Zero;
    valueFraction() = 1.0;
}


heatFluxRadiationFvPatchScalarField::heatFluxRadiationFvPatchScalarField
(
    const heatFluxRadiationFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    temperatureCoupledBase(tppsf),
    Q_(tppsf.Q_.clone()),
    Cp_(tppsf.Cp_),
    mass_(tppsf.mass_),
    relaxation_(tppsf.relaxation_),
    qrPrevious_(tppsf.qrPrevious_),
    qrRelaxation_(tppsf.qrRelaxation_),
    qrName_(tppsf.qrName_),
    minTemperature_(tppsf.minTemperature_),
    maxTemperature_(tppsf.maxTemperature_),
    curTimeIndex_(-1)
{}


heatFluxRadiationFvPatchScalarField::heatFluxRadiationFvPatchScalarField
(
    const heatFluxRadiationFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    temperatureCoupledBase(patch(), tppsf),
    Q_(tppsf.Q_.clone()),
    Cp_(tppsf.Cp_),
    mass_(tppsf.mass_),
    relaxation_(tppsf.relaxation_),
    qrPrevious_(tppsf.qrPrevious_),
    qrRelaxation_(tppsf.qrRelaxation_),
    qrName_(tppsf.qrName_),
    minTemperature_(tppsf.minTemperature_),
    maxTemperature_(tppsf.maxTemperature_),
    curTimeIndex_(-1)
{}


void heatFluxRadiationFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::autoMap(mapper);
    temperatureCoupledBase::autoMap(mapper);

    if (qrName_ != "none")
    {
        qrPrevious_.autoMap(mapper);
    }
}


void heatFluxRadiationFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const heatFluxRadiationFvPatchScalarField& rhs =
        refCast<const heatFluxRadiationFvPatchScalarField>(ptf);

    temperatureCoupledBase::rmap(rhs, addr);

    if (qrName_ != "none")
    {
        qrPrevious_.rmap(rhs.qrPrevious_, addr);
    }
}


void heatFluxRadiationFvPatchScalarField::updateCoeffs()
{
    if (updated() || (curTimeIndex_ == this->db().time().timeIndex()))
    {
        return;
    }

    scalarField& Tp(*this);
    const scalarField Tp0(Tp);
    const scalarField& magSf = patch().magSf();
    const scalar deltaT(db().time().deltaTValue());

    tmp<scalarField> tkappa(kappa(Tp));
    const scalarField qConv(tkappa.ref()*snGrad());

    scalarField qr(Tp.size(), Zero);
    if (qrName_ != "none")
    {
        qr = lerp
        (
            qrPrevious_,
            patch().lookupPatchField<volScalarField>(qrName_),
            qrRelaxation_
        );
        qrPrevious_ = qr;
    }

    const scalar QConv = gWeightedSum(magSf, qConv);
    const scalar QRad = gWeightedSum(magSf, qr);
    const scalar QPrescribed = Q_->value(this->db().time().timeOutputValue());

    const scalar dTp =
        relaxation_*deltaT*(QPrescribed + QRad - QConv)/(mass_*Cp_);

    Tp += dTp;

    if (minTemperature_ > -VGREAT || maxTemperature_ < VGREAT)
    {
        Tp.clamp_range(minTemperature_, maxTemperature_);
    }

    refGrad() = 0.0;
    refValue() = Tp;
    valueFraction() = 1.0;

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        auto limits = gMinMax(*this);
        const auto avg = gAverage(*this);

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " :"
            << " dT:" << gAverage(Tp - Tp0)
            << " Qprescribed [W]:" << QPrescribed
            << " Qrad [W]:" << QRad
            << " Qconv [W]:" << QConv
            << " wall temperature min:" << limits.min()
            << " max:" << limits.max()
            << " avg:" << avg
            << endl;
    }

    curTimeIndex_ = this->db().time().timeIndex();
}


void heatFluxRadiationFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchField<scalar>::write(os);
    temperatureCoupledBase::write(os);

    Q_->writeData(os);
    os.writeEntry("Cp", Cp_);
    os.writeEntry("mass", mass_);

    if (relaxation_ < 1)
    {
        os.writeEntry("relaxation", relaxation_);
    }

    os.writeEntry("qr", qrName_);
    if (qrName_ != "none")
    {
        os.writeEntry("qrRelaxation", qrRelaxation_);
        qrPrevious_.writeEntry("qrPrevious", os);
    }

    if (minTemperature_ > -VGREAT)
    {
        os.writeEntry("minTemperature", minTemperature_);
    }

    if (maxTemperature_ < VGREAT)
    {
        os.writeEntry("maxTemperature", maxTemperature_);
    }
}


makePatchTypeField
(
    fvPatchScalarField,
    heatFluxRadiationFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
