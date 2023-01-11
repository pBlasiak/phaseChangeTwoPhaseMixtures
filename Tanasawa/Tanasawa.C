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
    R_("R", dimGasConstant, phaseChangeTwoPhaseMixtureCoeffs_),
   	RintCoeff_{(2.0 - gamma_)*sqrt(2.0*M_PI*R_)/(2.0*gamma_*pow(satProps_->hEvap(),2)*rho2())},
	Rint_{RintCoeff_*pow(TSat(), 3./2)}
{
	Info<< "Tanasawa model settings:  " << endl;
	Info<< "gamma = "		  << gamma_ << endl;
	Info<< "R = "             << R_     << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseChangeTwoPhaseMixtures::Tanasawa::j()
{
	// Minus sign "-" to provide mc > 0  and mv < 0
	if (cond_)
	{
		jc_ = -neg(T() - TSat())*(T() - TSat())/Rint_/hEvap();
	}

	//if (evap_)
	//{
	//	mEvapNoAlphal_ = -mCoeff_*pos(T_ - TSat_)*(T_ - TSat_)*magGradLimitedAlphal_
	//		/sqrt(pow(TSat_,3.0));
	//	mEvapAlphal_   =  mEvapNoAlphal_*limitedAlphal_;
	//	mEvapNoTmTSat_ =  mCoeff_*pos(T_ - TSat_)*magGradLimitedAlphal_
	//		/sqrt(pow(TSat_,3.0))*limitedAlphal_;
	//}
}

Foam::tmp<Foam::volScalarField>
Foam::phaseChangeTwoPhaseMixtures::Tanasawa::alphaSourceSp() 
{
	alphaSourceSp_ = (je_ - jc_)*HW_->magGradLimitedAlphal();
	return alphaSourceSp_;
	//return Pair<tmp<volScalarField>>
	//(
	//	tmp<volScalarField>(-mCondAlphal_),
	//	tmp<volScalarField>(mEvapNoAlphal_)
    //);
}

Foam::tmp<Foam::volScalarField>
Foam::phaseChangeTwoPhaseMixtures::Tanasawa::alphaSourceSu() 
{
	alphaSourceSu_ =  jc_*HW_->magGradLimitedAlphal();
	return alphaSourceSu_;
	//return Pair<tmp<volScalarField>>
	//(
	//	tmp<volScalarField>(-mCondAlphal_),
	//	tmp<volScalarField>(mEvapNoAlphal_)
    //);
}

Foam::tmp<Foam::volScalarField>
Foam::phaseChangeTwoPhaseMixtures::Tanasawa::pSourceSp()
{
	return pSourceSp_;
	//return Pair<tmp<volScalarField> >
	//(
	// // New
    //    //-mCondAlphal_*neg(T_-TSat_)/max(p_-pSat_,1E-6*pSat_),
	//    //mEvapAlphal_*pos(T_-TSat_)/max(pSat_-p_,1E-6*pSat_)
	// // implicit implementation results in zero velocity at the interface
    //    //mCondAlphal_*pos(p_-pSat_)/max(p_-pSat_,1E-6*pSat_),
	//    //mEvapAlphal_*neg(p_-pSat_)/max(pSat_-p_,1E-6*pSat_)
	// // explicit implementation results in large errors
    //    -mCondAlphal_,
	//    mEvapAlphal_
	//);
}

Foam::tmp<Foam::volScalarField>
Foam::phaseChangeTwoPhaseMixtures::Tanasawa::pSourceSu()
{
	if (isHardtWondra())
	{
		HW_->spread(jc_, je_);
		return HW_->rhoSourcel(); 
	}
	else
	{
		pSourceSu_ = (jc_ - je_)*HW_->magGradLimitedAlphal();
		return pSourceSu_;
	}
}

Foam::tmp<Foam::volScalarField>
Foam::phaseChangeTwoPhaseMixtures::Tanasawa::TSourceSp() 
{
	return TSourceSp_;
}

Foam::tmp<Foam::volScalarField>
Foam::phaseChangeTwoPhaseMixtures::Tanasawa::TSourceSu() 
{
	if (isHardtWondra())
	{
		return HW_->hSourcel(); 
	}
	else
	{
		TSourceSu_ = (jc_ - je_)*HW_->magGradLimitedAlphal()*hEvap()*TSat();
		return TSourceSu_;
	}
}

void Foam::phaseChangeTwoPhaseMixtures::Tanasawa::correct()
{
	phaseChangeTwoPhaseMixture::correct();

	calcRint(); 
	j();
}

bool Foam::phaseChangeTwoPhaseMixtures::Tanasawa::read()
{
    if (phaseChangeTwoPhaseMixture::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("gamma") >> gamma_;

		RintCoeff_ = (2.0 - gamma_)*sqrt(2.0*M_PI*R_)/(2.0*gamma_*pow(satProps_->hEvap(),2)*rho2());

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
