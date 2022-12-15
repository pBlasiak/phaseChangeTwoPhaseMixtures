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

#include "Tanasawa.H"
#include "fvc.H"
#include "fvm.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

#include "volFields.H"
//#include "FaceCellWave.H"
//#include "smoothData.H"
//#include "sweepData.H"
//#include "fvMatrices.H"
//#include "fvcVolumeIntegrate.H"
//#include "fvmLaplacian.H"
//#include "fvmSup.H"
#include "zeroGradientFvPatchField.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(Tanasawa, 0);
    addToRunTimeSelectionTable(phaseChangeTwoPhaseMixture, Tanasawa, components);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixtures::Tanasawa::Tanasawa
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    phaseChangeTwoPhaseMixture(typeName, U, phi),

    gamma_("gamma", phaseChangeTwoPhaseMixtureCoeffs_.subDict(type() + "Coeffs")),
   	mCoeff_(2.0*gamma_/(2.0 - gamma_)/sqrt(2.0*M_PI*R_)*hEvap_*rho2())
{
	Info<< "Tanasawa model settings:  " << endl;
	Info<< "gamma = "		  << gamma_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Tanasawa::mDotAlphal() const
{
	return Pair<tmp<volScalarField>>
	(
		tmp<volScalarField>(mCondNoAlphal_),
		tmp<volScalarField>(mEvapNoAlphal_)
    );
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Tanasawa::mDotP() const
{
	return Pair<tmp<volScalarField> >
	(
	 // New
        //-mCondAlphal_*neg(T_-TSat_)/max(p_-pSat_,1E-6*pSat_),
	    //mEvapAlphal_*pos(T_-TSat_)/max(pSat_-p_,1E-6*pSat_)
	 // implicit implementation results in zero velocity at the interface
        //mCondAlphal_*pos(p_-pSat_)/max(p_-pSat_,1E-6*pSat_),
	    //mEvapAlphal_*neg(p_-pSat_)/max(pSat_-p_,1E-6*pSat_)
	 // explicit implementation results in large errors
        -mCondAlphal_,
	    mEvapAlphal_
	);
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Tanasawa::mDotT() const
{
	return Pair<tmp<volScalarField> >
	(
	    tmp<volScalarField>(mCondNoTmTSat_),
	    tmp<volScalarField>(mEvapNoTmTSat_)
	);
}

void Foam::phaseChangeTwoPhaseMixtures::Tanasawa::correct()
{
	phaseChangeTwoPhaseMixture::correct();

	calcMagGradLimitedAlphal();

	// In Tanasawa model there is no alpha term
	// probably it should be divided here by alphal and (1-alphal) but it
	// could produce errrors. To avoid this and follow the algorithm in alphaEqn.H
	// the modified Tanasawa model is implemented with additional multiplication
	// by (1-alphal) for condensation and alphal for evaporation.
	// Thus, the terms in pEqn and TEqn have to be also multiplied by these terms.
	// One can try to implement invAlpha which is zero for alphal = 0 otherwise
	// 1/alphal and multiply it by mCondAlphal_ and mEvapAlphal_ .
	// Minus sign "-" to provide mc > 0  and mv < 0
	if (cond_)
	{
		mCondNoAlphal_ = -mCoeff_*neg(T_ - TSat_)*(T_ - TSat_)*magGradLimitedAlphal_
			/sqrt(pow(TSat_,3.0));
		mCondAlphal_   =  mCondNoAlphal_*(1-limitedAlphal_);
		mCondNoTmTSat_ = -mCoeff_*neg(T_ - TSat_)*magGradLimitedAlphal_
			/sqrt(pow(TSat_,3.0))*(1-limitedAlphal_);
	}

	if (evap_)
	{
		mEvapNoAlphal_ = -mCoeff_*pos(T_ - TSat_)*(T_ - TSat_)*magGradLimitedAlphal_
			/sqrt(pow(TSat_,3.0));
		mEvapAlphal_   =  mEvapNoAlphal_*limitedAlphal_;
		mEvapNoTmTSat_ =  mCoeff_*pos(T_ - TSat_)*magGradLimitedAlphal_
			/sqrt(pow(TSat_,3.0))*limitedAlphal_;
	}
}

bool Foam::phaseChangeTwoPhaseMixtures::Tanasawa::read()
{
    if (phaseChangeTwoPhaseMixture::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("gamma") >> gamma_;

		mCoeff_ = 2.0*gamma_/(2.0 - gamma_)/sqrt(2.0*M_PI*R_)*hEvap_*rho2();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
