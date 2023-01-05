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
    phaseChangeTwoPhaseMixture(typeName, U, phi)

{
	Info<< "NichitaThome mass transfer model is used.  " << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
// tmp by trzeba dodac
Foam::volVectorField Foam::phaseChangeTwoPhaseMixtures::NichitaThome::calcGradAlphal() 
{
	calcLimitedAlphal();
	return fvc::grad(limitedAlphal_);
}

Foam::volVectorField Foam::phaseChangeTwoPhaseMixtures::NichitaThome::calcGradT() const
{
	return fvc::grad(T_);
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::NichitaThome::mDotAlphal() const
{
	return Pair<tmp<volScalarField>>
	(
		tmp<volScalarField>(mCondNoAlphal_),
		tmp<volScalarField>(mEvapNoAlphal_)
    );
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::NichitaThome::mDotP() const
{
	return Pair<tmp<volScalarField> >
	(
        mCondAlphal_*pos(p_-pSat_)/max(p_-pSat_,1E-6*pSat_),
	    mEvapAlphal_*neg(p_-pSat_)/max(pSat_-p_,1E-6*pSat_)
	);
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::NichitaThome::mDotT() const
{
    const dimensionedScalar T1("1K", dimTemperature, 1.0);
	return Pair<tmp<volScalarField> >
	(
	    mCondNoTmTSat_/max(TSat_ - T_, 1E-6*TSat_)*T1,
	    mEvapNoTmTSat_/max(T_ - TSat_, 1E-6*TSat_)*T1
	);
}

void Foam::phaseChangeTwoPhaseMixtures::NichitaThome::correct()
{
	const volScalarField gradAlphaGradT = calcGradAlphal() & calcGradT();
    const dimensionedScalar T1("1K", dimTemperature, 1.0);
	const volScalarField kEff = this->k();

	// In NichitaThome model there is no alpha term
	// probably it should be divided here by alphal and (1-alphal) but it
	// could produce errrors. To avoid this and follow the algorithm in alphaEqn.H
	// the modified NichitaThome model is implemented with additional multiplication
	// by (1-alphal) for condensation and alphal for evaporation.
	// Thus, the terms in pEqn and TEqn have to be also multiplied by these terms.
	if (cond_)
	{
		mCondNoAlphal_ = -neg(T_ - TSat_)*kEff*gradAlphaGradT/hEvap_;
		//mCondAlphal_ = mCondAlphal_;//*pos(scalar(1)-limitedAlpha1)/max(scalar(1)-limitedAlpha1, 1e-6);
		mCondAlphal_   = mCondNoAlphal_*(1-limitedAlphal_);
		mCondNoTmTSat_ = -mCondAlphal_/T1;
	}

	if (evap_)
	{
		mEvapNoAlphal_ =  pos(T_ - TSat_)*kEff*gradAlphaGradT/hEvap_;
		//mEvapAlphal_ = mEvapAlphal_;//*pos(limitedAlpha1)/max(limitedAlpha1, 1e-6);
		mEvapAlphal_   = mEvapNoAlphal_*limitedAlphal_;
		mEvapNoTmTSat_ =  mEvapAlphal_/T1;
	}
}

bool Foam::phaseChangeTwoPhaseMixtures::NichitaThome::read()
{
    if (phaseChangeTwoPhaseMixture::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");
        //phaseChangeTwoPhaseMixtureCoeffs_.lookup("condensation") >> cond_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
