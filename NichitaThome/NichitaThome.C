/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "NichitaThome.H"
#include "fvc.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(NichitaThome, 0);
    addToRunTimeSelectionTable(phaseChangeTwoPhaseMixture, NichitaThome, components);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixtures::NichitaThome::NichitaThome
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    phaseChangeTwoPhaseMixture(typeName, U, phi),

    cond_(phaseChangeTwoPhaseMixtureCoeffs_.subDict(type() + "Coeffs").lookup("condensation")),
    evap_(phaseChangeTwoPhaseMixtureCoeffs_.subDict(type() + "Coeffs").lookup("evaporation"))
{
	Info<< "NichitaThome model settings:  " << endl;
	Info<< "Condensation is " << cond_	<< endl;
	Info<< "Evaporation is "  << evap_  << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
// tmp by trzeba dodac
Foam::volVectorField Foam::phaseChangeTwoPhaseMixtures::NichitaThome::calcGradAlphal() const
{
	volScalarField limitedAlpha1 = min(max(alpha1(), scalar(0)), scalar(1));
	return fvc::grad(limitedAlpha1);
}

Foam::volVectorField Foam::phaseChangeTwoPhaseMixtures::NichitaThome::calcGradT() const
{
	return fvc::grad(T_);
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::NichitaThome::mDotAlphal() 
{
    volScalarField limitedAlpha1 = min(max(alpha1(), scalar(0)), scalar(1));
    const dimensionedScalar T1("1K", dimTemperature, 1.0);
	const volScalarField kEff = this->k();
	const volScalarField gradAlphaGradT = calcGradAlphal() & calcGradT();

	mCondNoAlphal_ = -neg(T_ - TSat_)*kEff*gradAlphaGradT/hEvap_;
	mEvapNoAlphal_ =  pos(T_ - TSat_)*kEff*gradAlphaGradT/hEvap_;


	// In NichitaThome model there is no alpha term
	// probably it should be divided here by alphal and (1-alphal) but it
	// could produce errrors. To avoid this and follow the algorithm in alphaEqn.H
	// the modified NichitaThome model is implemented with additional multiplication
	// by (1-alphal) for condensation and alphal for evaporation.
	// Thus, the terms in pEqn and TEqn have to be also multiplied by these terms.
	mCondAlphal_   = mCondNoAlphal_*(1-limitedAlpha1);
	mEvapAlphal_   = mEvapNoAlphal_*limitedAlpha1;
	//mCondAlphal_ = mCondAlphal_;//*pos(scalar(1)-limitedAlpha1)/max(scalar(1)-limitedAlpha1, 1e-6);
	//mEvapAlphal_ = mEvapAlphal_;//*pos(limitedAlpha1)/max(limitedAlpha1, 1e-6);

	mCondNoTmTSat_ = -mCondAlphal_/T1;
	mEvapNoTmTSat_ =  mEvapAlphal_/T1;

	if (cond_ && evap_)
	{
		return Pair<tmp<volScalarField> >
		(
	        mCondNoAlphal_*scalar(1),
		    mEvapNoAlphal_*scalar(1)
		);
	}
	else if (cond_)
	{
		return Pair<tmp<volScalarField> >
		(
	        mCondNoAlphal_*scalar(1),
		    mEvapNoAlphal_*scalar(0)
		);
	}
	else if (evap_)
	{
		return Pair<tmp<volScalarField> >
		(
	        mCondNoAlphal_*scalar(0),
		    mEvapNoAlphal_*scalar(1)
		);
	}
	else 
	{
		return Pair<tmp<volScalarField> >
		(
	        mCondNoAlphal_*scalar(0),
		    mEvapNoAlphal_*scalar(0)
		);
	}
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::NichitaThome::mDotP() const
{
	if (cond_ && evap_)
	{
		return Pair<tmp<volScalarField> >
		(
	        mCondAlphal_*pos(p_-pSat_)/max(p_-pSat_,1E-6*pSat_),
		    mEvapAlphal_*neg(p_-pSat_)/max(pSat_-p_,1E-6*pSat_)
		);
	}
	else if (cond_)
	{
		return Pair<tmp<volScalarField> >
		(
	        mCondAlphal_*pos(p_-pSat_)/max(p_-pSat_,1E-6*pSat_),
		    mEvapAlphal_*scalar(0)
						*neg(p_-pSat_)/max(pSat_-p_,1E-6*pSat_)
		);
	}
	else if (evap_)
	{
		return Pair<tmp<volScalarField> >
		(
	        mCondAlphal_*pos(p_-pSat_)/max(p_-pSat_,1E-6*pSat_)*scalar(0),
		    mEvapAlphal_*neg(p_-pSat_)/max(pSat_-p_,1E-6*pSat_)
		);
	}
	else 
	{
		return Pair<tmp<volScalarField> >
		(
	        mCondAlphal_*scalar(0)*pos(p_-pSat_)/max(p_-pSat_,1E-6*pSat_),
		    mEvapAlphal_*scalar(0)*neg(p_-pSat_)/max(pSat_-p_,1E-6*pSat_)
		);
	}
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::NichitaThome::mDotT() const
{
    const dimensionedScalar T1("1K", dimTemperature, 1.0);
	if (cond_ && evap_)
	{
		return Pair<tmp<volScalarField> >
		(
		    mCondNoTmTSat_/max(TSat_ - T_, 1E-6*TSat_)*T1,
		    mEvapNoTmTSat_/max(T_ - TSat_, 1E-6*TSat_)*T1
		);
	}
	else if (cond_)
	{
		return Pair<tmp<volScalarField> >
		(
		    mCondNoTmTSat_/max(TSat_ - T_, 1E-6*TSat_)*T1,
		    mEvapNoTmTSat_/max(T_ - TSat_, 1E-6*TSat_)*T1*scalar(0)
		);
	}
	else if (evap_)
	{
		return Pair<tmp<volScalarField> >
		(
		    mCondNoTmTSat_/max(TSat_ - T_, 1E-6*TSat_)*T1*scalar(0),
		    mEvapNoTmTSat_/max(T_ - TSat_, 1E-6*TSat_)*T1
		);
	}
	else
	{
		return Pair<tmp<volScalarField> >
		(
		    mCondNoTmTSat_/max(TSat_ - T_, 1E-6*TSat_)*T1*scalar(0),
		    mEvapNoTmTSat_/max(T_ - TSat_, 1E-6*TSat_)*T1*scalar(0)
		);
	}
}

bool Foam::phaseChangeTwoPhaseMixtures::NichitaThome::read()
{
    if (phaseChangeTwoPhaseMixture::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("condensation") >> cond_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("evaporation") >> evap_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
