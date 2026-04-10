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
#include "physicoChemicalConstants.H"

using Foam::constant::physicoChemical::sigma;

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
    q_(nullptr),
    h_(nullptr),
    Ta_(nullptr),
    Cp_(0),
    mass_(0),
    emissivity_(0),
    relaxation_(1),
    qrPrevious_(),
    qrRelaxation_(1),
    qrName_("none"),
    minTemperature_(-VGREAT),
    maxTemperature_(VGREAT),
    thicknessLayers_(),
    kappaLayers_(),
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
    Q_(ptf.Q_.valid() ? ptf.Q_.clone() : nullptr),
    q_(ptf.q_.valid() ? ptf.q_.clone(patch().patch()) : nullptr),
    h_(ptf.h_.valid() ? ptf.h_.clone(patch().patch()) : nullptr),
    Ta_(ptf.Ta_.valid() ? ptf.Ta_.clone() : nullptr),
    Cp_(ptf.Cp_),
    mass_(ptf.mass_),
    emissivity_(ptf.emissivity_),
    relaxation_(ptf.relaxation_),
    qrPrevious_(),
    qrRelaxation_(ptf.qrRelaxation_),
    qrName_(ptf.qrName_),
    minTemperature_(ptf.minTemperature_),
    maxTemperature_(ptf.maxTemperature_),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_),
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
    Q_(dict.found("Q") ? Function1<scalar>::New("Q", dict, &db()) : nullptr),
    q_
    (
        dict.found("q")
      ? PatchFunction1<scalar>::New(patch().patch(), "q", dict)
      : nullptr
    ),
    h_
    (
        dict.found("h")
      ? PatchFunction1<scalar>::New(patch().patch(), "h", dict)
      : nullptr
    ),
    Ta_(dict.found("Ta") ? Function1<scalar>::New("Ta", dict, &db()) : nullptr),
    Cp_(dict.get<scalar>("Cp")),
    mass_(dict.get<scalar>("mass")),
    emissivity_(dict.getOrDefault<scalar>("emissivity", 0)),
    relaxation_(dict.getOrDefault<scalar>("relaxation", 1)),
    qrPrevious_(),
    qrRelaxation_(dict.getOrDefault<scalar>("qrRelaxation", 1)),
    qrName_(dict.getOrDefault<word>("qr", "none")),
    minTemperature_(dict.getOrDefault<scalar>("minTemperature", -VGREAT)),
    maxTemperature_(dict.getOrDefault<scalar>("maxTemperature", VGREAT)),
    thicknessLayers_(),
    kappaLayers_(),
    curTimeIndex_(-1)
{
    if (!Q_.valid() && !q_.valid())
    {
        FatalIOErrorInFunction(dict)
            << "Either Q [W] or q [W/m^2] must be provided for patch "
            << p.name()
            << exit(FatalIOError);
    }

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

    if (h_.valid() != Ta_.valid())
    {
        FatalIOErrorInFunction(dict)
            << "Both h [W/m^2/K] and Ta [K] must be provided together for "
            << "patch " << p.name()
            << exit(FatalIOError);
    }

    if (dict.readIfPresent("thicknessLayers", thicknessLayers_))
    {
        if (!h_.valid())
        {
            FatalIOErrorInFunction(dict)
                << "thicknessLayers/kappaLayers require h and Ta for patch "
                << p.name()
                << exit(FatalIOError);
        }

        dict.readEntry("kappaLayers", kappaLayers_);

        if (thicknessLayers_.size() != kappaLayers_.size())
        {
            FatalIOErrorInFunction(dict)
                << "number of layers for thicknessLayers and kappaLayers "
                << "must be the same for patch " << p.name()
                << exit(FatalIOError);
        }
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
    Q_(tppsf.Q_.valid() ? tppsf.Q_.clone() : nullptr),
    q_(tppsf.q_.valid() ? tppsf.q_.clone(patch().patch()) : nullptr),
    h_(tppsf.h_.valid() ? tppsf.h_.clone(patch().patch()) : nullptr),
    Ta_(tppsf.Ta_.valid() ? tppsf.Ta_.clone() : nullptr),
    Cp_(tppsf.Cp_),
    mass_(tppsf.mass_),
    emissivity_(tppsf.emissivity_),
    relaxation_(tppsf.relaxation_),
    qrPrevious_(tppsf.qrPrevious_),
    qrRelaxation_(tppsf.qrRelaxation_),
    qrName_(tppsf.qrName_),
    minTemperature_(tppsf.minTemperature_),
    maxTemperature_(tppsf.maxTemperature_),
    thicknessLayers_(tppsf.thicknessLayers_),
    kappaLayers_(tppsf.kappaLayers_),
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
    Q_(tppsf.Q_.valid() ? tppsf.Q_.clone() : nullptr),
    q_(tppsf.q_.valid() ? tppsf.q_.clone(patch().patch()) : nullptr),
    h_(tppsf.h_.valid() ? tppsf.h_.clone(patch().patch()) : nullptr),
    Ta_(tppsf.Ta_.valid() ? tppsf.Ta_.clone() : nullptr),
    Cp_(tppsf.Cp_),
    mass_(tppsf.mass_),
    emissivity_(tppsf.emissivity_),
    relaxation_(tppsf.relaxation_),
    qrPrevious_(tppsf.qrPrevious_),
    qrRelaxation_(tppsf.qrRelaxation_),
    qrName_(tppsf.qrName_),
    minTemperature_(tppsf.minTemperature_),
    maxTemperature_(tppsf.maxTemperature_),
    thicknessLayers_(tppsf.thicknessLayers_),
    kappaLayers_(tppsf.kappaLayers_),
    curTimeIndex_(-1)
{}


void heatFluxRadiationFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::autoMap(mapper);
    temperatureCoupledBase::autoMap(mapper);

    if (q_.valid())
    {
        q_->autoMap(mapper);
    }

    if (h_.valid())
    {
        h_->autoMap(mapper);
    }

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

    if (q_.valid())
    {
        q_->rmap(rhs.q_(), addr);
    }

    if (h_.valid())
    {
        h_->rmap(rhs.h_(), addr);
    }

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
    const scalar t = this->db().time().timeOutputValue();

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
    scalar QAmbient = Zero;
    scalar QPrescribed = Zero;

    if (Q_.valid())
    {
        QPrescribed += Q_->value(t);
    }

    if (q_.valid())
    {
        tmp<scalarField> tqFlux(q_->value(t));
        const scalarField& qFlux = tqFlux();
        QPrescribed += gWeightedSum(magSf, qFlux);
    }

    if (h_.valid())
    {
        tmp<scalarField> thtcCoeff(h_->value(t));
        const scalarField& htcCoeff = thtcCoeff();
        const scalar Ta = Ta_->value(t);

        scalar totalSolidRes = 0;
        if (thicknessLayers_.size())
        {
            forAll(thicknessLayers_, iLayer)
            {
                const scalar l = thicknessLayers_[iLayer];
                if (kappaLayers_[iLayer] > 0)
                {
                    totalSolidRes += l/kappaLayers_[iLayer];
                }
            }
        }

        scalarField hrad(Tp.size(), Zero);
        if (emissivity_ > 0)
        {
            const scalar eSig(emissivity_*sigma.value());

            if (totalSolidRes > 0)
            {
                scalarField TpLambda(htcCoeff/(htcCoeff + 1/totalSolidRes));
                scalarField Ts(TpLambda*Ta + (1 - TpLambda)*Tp);

                hrad = eSig*((pow3(Ta) + pow3(Ts)) + Ta*Ts*(Ta + Ts));

                forAll(hrad, i)
                {
                    scalar hradTmp0 = hrad[i];
                    scalar TaLambda =
                        (htcCoeff[i] + hradTmp0)
                       /(htcCoeff[i] + hradTmp0 + 1/totalSolidRes);

                    scalar TsiNew = TaLambda*Ta + (1 - TaLambda)*Tp[i];
                    scalar Tsi = Ts[i];

                    while (mag(Tsi - TsiNew)/max(Tsi, SMALL) > 0.01)
                    {
                        Tsi = TsiNew;
                        scalar hradNew
                        (
                            eSig
                           *(
                                (pow3(Ta) + pow3(Tsi))
                              + Ta*Tsi*(Ta + Tsi)
                            )
                        );

                        TaLambda =
                            (htcCoeff[i] + hradNew)
                           /(htcCoeff[i] + hradNew + 1/totalSolidRes);

                        TsiNew = TaLambda*Ta + (1 - TaLambda)*Tp[i];
                    }

                    hrad[i] =
                        eSig*((pow3(Ta) + pow3(Tsi)) + Ta*Tsi*(Ta + Tsi));
                }
            }
            else
            {
                hrad = eSig*((pow3(Ta) + pow3(Tp)) + Ta*Tp*(Ta + Tp));
            }
        }

        const scalarField hp(1/(1/(htcCoeff + hrad) + totalSolidRes));
        const scalarField qAmbient(hp*(Ta - Tp));

        QAmbient = gWeightedSum(magSf, qAmbient);
    }

    const scalar dTp =
        relaxation_*deltaT*(QPrescribed + QRad + QAmbient - QConv)
       /(mass_*Cp_);

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
            << " Qambient [W]:" << QAmbient
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

    if (Q_.valid())
    {
        Q_->writeData(os);
    }
    if (q_.valid())
    {
        q_->writeData(os);
    }
    if (h_.valid())
    {
        h_->writeData(os);
        Ta_->writeData(os);

        if (emissivity_ > 0)
        {
            os.writeEntry("emissivity", emissivity_);
        }

        if (thicknessLayers_.size())
        {
            thicknessLayers_.writeEntry("thicknessLayers", os);
            kappaLayers_.writeEntry("kappaLayers", os);
        }
    }
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
