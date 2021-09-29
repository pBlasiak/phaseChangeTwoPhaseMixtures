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
//#include "fvc.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

//#include "fvCFD.H"
#include "fvcSmooth.H"
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

    cond_(phaseChangeTwoPhaseMixtureCoeffs_.subDict(type() + "Coeffs").lookup("condensation")),
    evap_(phaseChangeTwoPhaseMixtureCoeffs_.subDict(type() + "Coeffs").lookup("evaporation")),
    gamma_(phaseChangeTwoPhaseMixtureCoeffs_.subDict(type() + "Coeffs").lookup("gamma")),
   	mCoeff_(2.0*gamma_/(2.0 - gamma_)/sqrt(2.0*M_PI*R_)*hEvap_*rho2())
{
	Info<< "Tanasawa model settings:  " << endl;
	Info<< "Condensation is " << cond_	<< endl;
	Info<< "Evaporation is "  << evap_  << endl;
	Info<< "gamma = "		  << gamma_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::volScalarField Foam::phaseChangeTwoPhaseMixtures::Tanasawa::calcGradAlphal() const
{
	volScalarField limitedAlpha1 = min(max(alpha1(), scalar(0)), scalar(1));
	return mag(fvc::grad(limitedAlpha1));
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Tanasawa::mDotAlphal()
{
    //const dimensionedScalar T0("0", dimTemperature, 0.0);
    volScalarField limitedAlpha1 = min(max(alpha1(), scalar(0)), scalar(1));
	volScalarField gradAlphal = mag(fvc::grad(limitedAlpha1));

//    volScalarField psi
//    (
//        IOobject
//        (
//            "psi",
//            alpha1().time().timeName(),
//            alpha1().db(),
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//		alpha1().mesh(),
//        dimensionedScalar("psi0", dimensionSet(1,-3,-1,0,0,0,0), 0 ),
//		"zeroGradient"
//    );
//	
	//dimensionedScalar DPsi = phaseChangeTwoPhaseMixtureCoeffs_.lookup("DPsi");
	scalar spread_ = 3;
    dimensionedScalar DPsi
    (
        "D",
        dimArea,
        spread_/sqr(gAverage(alpha1().mesh().nonOrthDeltaCoeffs()))
    );
////- alpha cutoff value for source term distribution 
////  (no source terms in cells with cutoff < alpha1 < 1-cutoff)
scalar cutoff = 1e-3;
//
//////- Calculation of source terms according to Hardt/Wondra (JCP, 2008)
////volScalarField psi0Tild = mag(fvc::grad(alpha1));
////
////dimensionedScalar intPsi0Tild = fvc::domainIntegrate(psi0Tild);
////dimensionedScalar intAlphaPsi0Tild = fvc::domainIntegrate(alpha1*psi0Tild);
////
////dimensionedScalar N ("N", dimensionSet(0,0,0,0,0,0,0), 2.0);
////if (intAlphaPsi0Tild.value() > 1e-99)
////{
////	N = intPsi0Tild/intAlphaPsi0Tild;
////}
////
////volScalarField jEvap = (T-Tsat)/(Rph*hEvap);
////
////psi0 = N*jEvap*alpha1*psi0Tild;
//
//
//
//
	// minus sign "-" to provide mc > 0  and mv < 0
	mCondNoAlphal_ = -mCoeff_*neg(T_ - TSat_)*(T_ - TSat_)*gradAlphal/sqrt(pow(TSat_,3.0));
	mEvapNoAlphal_ = -mCoeff_*pos(T_ - TSat_)*(T_ - TSat_)*gradAlphal/sqrt(pow(TSat_,3.0));

	fvc::spreadSource
	(
		mEvapNoAlphal_,
		mEvapNoAlphal_,
		alpha1(),
		alpha2(),
		DPsi,
		cutoff
	);
//
//dimensionedScalar intPsi0 = fvc::domainIntegrate(mEvapNoAlphal_);
////- Smearing of source term field
//fvScalarMatrix psiEqn
//(
//	fvm::Sp(scalar(1),psi) - fvm::laplacian(DPsi,psi) == mEvapNoAlphal_
//);
//
//Info<< "SOLVING" << endl;
//
//psiEqn.solve();
//
////- Cut cells with cutoff < alpha1 < 1-cutoff and rescale remaining source term field
//dimensionedScalar intPsiLiquid ("intPsiLiquid", dimensionSet(1,0,-1,0,0,0,0), 0.0);
//dimensionedScalar intPsiVapor ("intPsiLiquid", dimensionSet(1,0,-1,0,0,0,0), 0.0);
//
//forAll(alpha1().mesh().C(),iCell)
//{
//	if (limitedAlpha1[iCell] < cutoff)
//	{
//		intPsiVapor.value() += (1.0-limitedAlpha1[iCell])*psi[iCell]*alpha1().mesh().V()[iCell];
//	}
//	else if (limitedAlpha1[iCell] > 1.0-cutoff)
//	{
//		intPsiLiquid.value() += limitedAlpha1[iCell]*psi[iCell]*alpha1().mesh().V()[iCell];
//	}
//}
//
////- Calculate Nl and Nv
//dimensionedScalar Nl ("Nl", dimensionSet(0,0,0,0,0,0,0), 2.0);
//dimensionedScalar Nv ("Nv", dimensionSet(0,0,0,0,0,0,0), 2.0);
//
//reduce(intPsiLiquid.value(),sumOp<scalar>());
//reduce(intPsiVapor.value(),sumOp<scalar>());
//
//if (intPsiLiquid.value() > 1e-99)
//{
//    Nl = intPsi0/intPsiLiquid;
//}
//if (intPsiVapor.value() > 1e-99)
//{
//    Nv = intPsi0/intPsiVapor;
//}
//
//        
////- Set source terms in cells with alpha1 < cutoff or alpha1 > 1-cutoff
//forAll(alpha1().mesh().C(),iCell)
//{
//	if (limitedAlpha1[iCell] < cutoff)
//	{
//		mEvapNoAlphal_[iCell] = Nv.value()*(1.0-limitedAlpha1[iCell])*psi[iCell];
//	}
//	else if (limitedAlpha1[iCell] > 1.0-cutoff)
//	{
//		mEvapNoAlphal_[iCell] = -Nl.value()*limitedAlpha1[iCell]*psi[iCell];
//	}
//	else
//	{
//		mEvapNoAlphal_[iCell] = 0.0;
//	}
//}
	// In Tanasawa model there is no alpha term
	// probably it should be divided here by alphal and (1-alphal) but it
	// could produce errrors. To avoid this and follow the algorithm in alphaEqn.H
	// the modified Tanasawa model is implemented with additional multiplication
	// by (1-alphal) for condensation and alphal for evaporation.
	// Thus, the terms in pEqn and TEqn have to be also multiplied by these terms.
	// One can try to implement invAlpha which is zero for alphal = 0 otherwise
	// 1/alphal and multiply it by mCondAlphal_ and mEvapAlphal_
	mCondAlphal_   = mCondNoAlphal_*(1-limitedAlpha1);
	mEvapAlphal_   = mEvapNoAlphal_*limitedAlpha1;

	// minus sign to provide mc > 0  and mv > 0
	mCondNoTmTSat_ = -mCoeff_*neg(T_ - TSat_)*gradAlphal/sqrt(pow(TSat_,3.0))*(1-limitedAlpha1);
	mEvapNoTmTSat_ =  mCoeff_*pos(T_ - TSat_)*gradAlphal/sqrt(pow(TSat_,3.0))*limitedAlpha1;

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
Foam::phaseChangeTwoPhaseMixtures::Tanasawa::mDotP() const
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
Foam::phaseChangeTwoPhaseMixtures::Tanasawa::mDotT() const
{
	if (cond_ && evap_)
	{
		return Pair<tmp<volScalarField> >
		(
		    mCondNoTmTSat_*scalar(1),
		    mEvapNoTmTSat_*scalar(1)
		);
	}
	else if (cond_)
	{
		return Pair<tmp<volScalarField> >
		(
		    mCondNoTmTSat_*scalar(1),
		    mEvapNoTmTSat_*scalar(0)
		);
	}
	else if (evap_)
	{
		return Pair<tmp<volScalarField> >
		(
		    mCondNoTmTSat_*scalar(0),
		    mEvapNoTmTSat_*scalar(1)
		);
	}
	else
	{
		return Pair<tmp<volScalarField> >
		(
		    mCondNoTmTSat_*scalar(0),
		    mEvapNoTmTSat_*scalar(0)
		);
	}
}

bool Foam::phaseChangeTwoPhaseMixtures::Tanasawa::read()
{
    if (phaseChangeTwoPhaseMixture::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("condensation") >> cond_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("evaporation") >> evap_;
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
