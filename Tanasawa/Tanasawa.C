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

    cond_("condensation", phaseChangeTwoPhaseMixtureCoeffs_.subDict(type() + "Coeffs")),
    evap_("evaporation", phaseChangeTwoPhaseMixtureCoeffs_.subDict(type() + "Coeffs")),
    gamma_("gamma", phaseChangeTwoPhaseMixtureCoeffs_.subDict(type() + "Coeffs")),
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
        mCondAlphal_*pos(p_-pSat_)/max(p_-pSat_,1E-6*pSat_),
	    mEvapAlphal_*neg(p_-pSat_)/max(pSat_-p_,1E-6*pSat_)
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

    //const dimensionedScalar T0("0", dimTemperature, 0.0);
    volScalarField limitedAlpha1 = min(max(alpha1(), scalar(0)), scalar(1));
//	volScalarField gradAlphal = mag(fvc::grad(limitedAlpha1));
volScalarField psi0Tild = mag(fvc::grad(limitedAlpha1));

scalar cutoff = 1e-3;

dimensionedScalar intPsi0Tild = fvc::domainIntegrate(psi0Tild);
dimensionedScalar intAlphaPsi0Tild = fvc::domainIntegrate(limitedAlpha1*psi0Tild);

dimensionedScalar N ("N", dimensionSet(0,0,0,0,0,0,0), 2.0);
if (intAlphaPsi0Tild.value() > 1e-99)
{
	N = intPsi0Tild/intAlphaPsi0Tild;
}


//psi0 = N*jEvap*alpha1*psi0Tild;

	// minus sign "-" to provide mc > 0  and mv < 0
	if (cond_)
	{
		mCondNoAlphal_ = -N*mCoeff_*neg(T_ - TSat_)*(T_ - TSat_)*psi0Tild/sqrt(pow(TSat_,3.0));
		// In Tanasawa model there is no alpha term
		// probably it should be divided here by alphal and (1-alphal) but it
		// could produce errrors. To avoid this and follow the algorithm in alphaEqn.H
		// the modified Tanasawa model is implemented with additional multiplication
		// by (1-alphal) for condensation and alphal for evaporation.
		// Thus, the terms in pEqn and TEqn have to be also multiplied by these terms.
		// One can try to implement invAlpha which is zero for alphal = 0 otherwise
		// 1/alphal and multiply it by mCondAlphal_ and mEvapAlphal_
		mCondAlphal_   =  mCondNoAlphal_*(1-limitedAlpha1);
		mCondNoTmTSat_ = -N*mCoeff_*neg(T_ - TSat_)*psi0Tild/sqrt(pow(TSat_,3.0))*(1-limitedAlpha1);
	}

	if (evap_)
	{
		mEvapNoAlphal_ = -N*mCoeff_*pos(T_ - TSat_)*(T_ - TSat_)*psi0Tild/sqrt(pow(TSat_,3.0));
		mEvapAlphal_   =  mEvapNoAlphal_*limitedAlpha1;
		mEvapNoTmTSat_ =  N*mCoeff_*pos(T_ - TSat_)*psi0Tild/sqrt(pow(TSat_,3.0))*limitedAlpha1;
	}
dimensionedScalar intPsi0l = fvc::domainIntegrate(mCondAlphal_);
dimensionedScalar intPsi0v = fvc::domainIntegrate(mEvapAlphal_);

//volScalarField psil = mCondAlphal_;
//volScalarField psiv = mEvapAlphal_;
    //volScalarField psil
    //(
    //    IOobject
    //    (
    //        "psil",
    //        runTime.timeName(),
    //        mesh,
    //        IOobject::NO_READ,
    //        IOobject::NO_WRITE
    //    ),
	//	mCondAlphal_
    //);
const fvMesh& mesh = T_.mesh();
    volScalarField psiv
    (
        IOobject
        (
            "psiv",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(mEvapAlphal_.dimensions(), Zero),
        zeroGradientFvPatchField<scalar>::typeName
    );
//
//- Smearing of source term field
    dimensionedScalar DPsi
    (
        "DPsi",
        dimArea,
        cutoff/sqr(gAverage(T_.mesh().nonOrthDeltaCoeffs()))
    );
//fvScalarMatrix psilEqn
//(
//	fvm::Sp(scalar(1),psil) - fvm::laplacian(DPsi,psil) == mCondAlphal_
//);
//
//psilEqn.solve();

fvScalarMatrix psivEqn
(
	fvm::Sp(scalar(1),psiv) - fvm::laplacian(DPsi,psiv) == mEvapAlphal_
);

psivEqn.solve();

//- Cut cells with cutoff < alpha1 < 1-cutoff and rescale remaining source term field
dimensionedScalar intPsiLiquidEvaporation ("intPsiLiquidEvaporation", dimensionSet(1,0,-1,0,0,0,0), 0.0);
dimensionedScalar intPsiVaporEvaporation ("intPsiLiquidEvaporation", dimensionSet(1,0,-1,0,0,0,0), 0.0);


forAll(mesh.C(),iCell)
{
	if (limitedAlpha1[iCell] < cutoff)
	{
		intPsiVaporEvaporation.value() += (1.0-limitedAlpha1[iCell])*psiv[iCell]*mesh.V()[iCell];
	}
	else if (limitedAlpha1[iCell] > 1.0-cutoff)
	{
		intPsiLiquidEvaporation.value() += limitedAlpha1[iCell]*psiv[iCell]*mesh.V()[iCell];
	}
}

//- Calculate Nl and Nv
dimensionedScalar Nl ("Nl", dimensionSet(0,0,0,0,0,0,0), 2.0);
dimensionedScalar Nv ("Nv", dimensionSet(0,0,0,0,0,0,0), 2.0);

reduce(intPsiLiquidEvaporation.value(),sumOp<scalar>());
reduce(intPsiVaporEvaporation.value(),sumOp<scalar>());

if (intPsiLiquidEvaporation.value() > 1e-99)
{
    Nl = intPsi0v/intPsiLiquidEvaporation;
}
if (intPsiVaporEvaporation.value() > 1e-99)
{
    Nv = intPsi0v/intPsiVaporEvaporation;
}

        
//- Set source terms in cells with alpha1 < cutoff or alpha1 > 1-cutoff
forAll(mesh.C(),iCell)
{
	if (limitedAlpha1[iCell] < cutoff)
	{
		mEvapAlphal_[iCell] = Nv.value()*(1.0-limitedAlpha1[iCell])*psiv[iCell];
	}
	else if (limitedAlpha1[iCell] > 1.0-cutoff)
	{
		mEvapAlphal_[iCell] = -Nl.value()*limitedAlpha1[iCell]*psiv[iCell];
	}
	else
	{
		mEvapAlphal_[iCell] = 0.0;
	}
}

//- Evaporation source term in energy equation
//hESource = -N*alpha1*psi0Tild/Rph;
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
