/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "phaseChangeTwoPhaseMixture.H"
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
    defineTypeNameAndDebug(phaseChangeTwoPhaseMixture, 0);
    defineRunTimeSelectionTable(phaseChangeTwoPhaseMixture, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixture::phaseChangeTwoPhaseMixture
(
    const word& type,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    thermalIncompressibleTwoPhaseMixture(U, phi),
    //phaseChangeTwoPhaseMixtureCoeffs_(subDict(type + "Coeffs")),
    phaseChangeTwoPhaseMixtureCoeffs_
	(
	    IOdictionary
	    (
	        IOobject
            (
                "phaseChangeProperties", // dictionary name
                U.time().constant(),     // dict is found in "constant"
                U.db(),                  // registry for the dict
                IOobject::MUST_READ,     // must exist, otherwise failure
                IOobject::NO_WRITE       // dict is only read by the solver
            )
	    )
	),
	HW_(phaseChangeTwoPhaseMixtureCoeffs_.lookupOrDefault("HardtWondra", true)),
	cutoff_(phaseChangeTwoPhaseMixtureCoeffs_.lookupOrDefault("cutoff", 1e-3)),
	limitedAlphalCalculated_(false),
	magGradLimitedAlphalCalculated_(false),
    correctionTerm_
    (
        IOobject
        (
            "correctionTerm",
            U.time().timeName(),
            U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U.mesh(),
        //dimensionedScalar("correctionTerm", dimensionSet(1, 0, -1, 0, 0, 0, 0), 0.0)
        dimensionedScalar("correctionTerm", dimensionSet(1,-1,-3,-1,0,0,0), 0.0)
    ),
    //cond_("condensation", phaseChangeTwoPhaseMixtureCoeffs_.subDict(type() + "Coeffs")),
    //evap_("evaporation", phaseChangeTwoPhaseMixtureCoeffs_.subDict(type() + "Coeffs")),
    cond_("condensation", phaseChangeTwoPhaseMixtureCoeffs_),
    evap_("evaporation", phaseChangeTwoPhaseMixtureCoeffs_),
	limitedAlphal_(min(max(alpha1(), scalar(0)), scalar(1))),
	magGradLimitedAlphal_(mag(fvc::grad(limitedAlphal_))),
    mCond_
    (
        IOobject
        (
            "mCond",
            U.time().timeName(),
            U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("mCond", dimensionSet(1, 0, -1, 0, 0, 0, 0), 0.0)
    ),
    mEvap_
    (
        IOobject
        (
            "mEvap",
            U.time().timeName(),
            U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("mEvap", dimensionSet(1, 0, -1, 0, 0, 0, 0), 0.0)
    ),
	mCondAlphal_
	(
	    IOobject
	    (
	        "mCondAlphal",
	        U.time().timeName(),
	        U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
	    ),
	    U.mesh(),
	    dimensionedScalar("mCondAlphal", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
	),
	mEvapAlphal_
	(
	    IOobject
	    (
	        "mEvapAlphal",
	        U.time().timeName(),
	        U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
	    ),
	    U.mesh(),
	    dimensionedScalar("mEvapAlphal", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
	),
	mCondNoAlphal_
	(
	    IOobject
	    (
	        "mCondNoAlphal",
	        U.time().timeName(),
	        U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
	    ),
	    U.mesh(),
	    dimensionedScalar("mCondNoAlphal", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
	),
	mEvapNoAlphal_
	(
	    IOobject
	    (
	        "mEvapNoAlphal",
	        U.time().timeName(),
	        U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
	    ),
	    U.mesh(),
	    dimensionedScalar("mEvapNoAlphal", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
	),
	mCondNoTmTSat_
	(
	    IOobject
	    (
	        "mCondNoTmTSat",
	        U.time().timeName(),
	        U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
	    ),
	    U.mesh(),
	    dimensionedScalar("mCondNoTmTSat", dimensionSet(1, -3, -1, -1, 0, 0, 0), 0.0)
	),
	mEvapNoTmTSat_
	(
	    IOobject
	    (
	        "mEvapNoTmTSat",
	        U.time().timeName(),
	        U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
	    ),
	    U.mesh(),
	    dimensionedScalar("mEvapNoTmTSat", dimensionSet(1, -3, -1, -1, 0, 0, 0), 0.0)
	),
	mCondP_
	(
	    IOobject
	    (
	        "mCondP",
	        U.time().timeName(),
	        U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
	    ),
	    U.mesh(),
	    dimensionedScalar("mCondP", dimensionSet(0, -2, 1, 0, 0, 0, 0), 0.0)
	),
	mEvapP_
	(
	    IOobject
	    (
	        "mEvapP",
	        U.time().timeName(),
	        U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
	    ),
	    U.mesh(),
	    dimensionedScalar("mEvapP", dimensionSet(0, -2, 1, 0, 0, 0, 0), 0.0)
	),
	mCondT_
	(
	    IOobject
	    (
	        "mCondT",
	        U.time().timeName(),
	        U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
	    ),
	    U.mesh(),
	    dimensionedScalar("mCondT", dimensionSet(1, -3, -1, -1, 0, 0, 0), 0.0)
	),
	mEvapT_
	(
	    IOobject
	    (
	        "mEvapT",
	        U.time().timeName(),
	        U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
	    ),
	    U.mesh(),
	    dimensionedScalar("mEvapT", dimensionSet(1, -3, -1, -1, 0, 0, 0), 0.0)
	),
	p_(U.db().lookupObject<volScalarField>("p")),
	T_(U.db().lookupObject<volScalarField>("T")),
    TSatG_("TSatGlobal", dimTemperature, phaseChangeTwoPhaseMixtureCoeffs_),
    TSat_
    (
        IOobject
        (
            "TSat",
            U.time().timeName(),
            U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U.mesh(),
		TSatG_
    ),
    pSat_("pSat", dimPressure, phaseChangeTwoPhaseMixtureCoeffs_),
    hEvap_("hEvap", dimEnergy/dimMass, phaseChangeTwoPhaseMixtureCoeffs_),
    R_("R", dimGasConstant, phaseChangeTwoPhaseMixtureCoeffs_),
    TSatLocalPressure_(readBool(phaseChangeTwoPhaseMixtureCoeffs_.lookup("TSatLocalPressure"))),
	printPhaseChange_(readBool(phaseChangeTwoPhaseMixtureCoeffs_.lookup("printPhaseChange")))
{
	Info<< "TSatGlobal = "				 << TSatG_ << endl;
	Info<< "pSat = "		  			 << pSat_ << endl;
	Info<< "hEvap = "		  			 << hEvap_ << endl;
	Info<< "R = "			  			 << R_ << endl;
	Info<< "TSatLocalPressure = "        << TSatLocalPressure_ << endl;
	Info<< "printPhaseChange = "         << printPhaseChange_ << endl;

	Info<< "Condensation is " << cond_	 << endl;
	Info<< "Evaporation is "  << evap_   << endl;
	Info<< "Hardt-Wondra algorithm is: " << HW_ << endl;
	Info<< "Cutoff is set as: "          << cutoff_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::phaseChangeTwoPhaseMixture::calcLimitedAlphal()
{
	if(!limitedAlphalCalculated_)
	{
		limitedAlphal_ = min(max(alpha1(), scalar(0)), scalar(1));
		limitedAlphalCalculated_ = true;
	}
}

void Foam::phaseChangeTwoPhaseMixture::calcMagGradLimitedAlphal()
{
	if (!magGradLimitedAlphalCalculated_)
	{
		calcLimitedAlphal();
		magGradLimitedAlphal_ = mag(fvc::grad(limitedAlphal_));
		magGradLimitedAlphalCalculated_ = true;
	}
}

void Foam::phaseChangeTwoPhaseMixture::HardtWondra()
{
	if(!HW_)
	{
		magGradLimitedAlphalCalculated_ = false;
		return;
	}
	else
	{
		calcMagGradLimitedAlphal();

		//- Smearing of source term field
		dimensionedScalar DPsi
		(
		    "DPsi",
		    dimArea,
		    cutoff_/sqr(gAverage(T_.mesh().nonOrthDeltaCoeffs()))
		);

		dimensionedScalar intPsi0Tild = fvc::domainIntegrate(magGradLimitedAlphal_);
		dimensionedScalar intAlphaPsi0Tild = fvc::domainIntegrate((limitedAlphal_)*magGradLimitedAlphal_);
			
		dimensionedScalar N("N", dimensionSet(0,0,0,0,0,0,0), 2.0);
		if (intAlphaPsi0Tild.value() > 1e-99)
		{
			N = intPsi0Tild/intAlphaPsi0Tild;
		}

		if (cond_)
		{
			mCondNoAlphal_ *= N;
			mCondNoTmTSat_ *= N;

			dimensionedScalar intPsi0l = fvc::domainIntegrate(mCondAlphal_);

			const fvMesh& mesh = T_.mesh();
			volScalarField psil
			(
			    IOobject
			    (
			        "psil",
			        mesh.time().timeName(),
			        mesh,
			        IOobject::NO_READ,
			        IOobject::NO_WRITE
			    ),
			    mesh,
			    dimensionedScalar(mCondAlphal_.dimensions(), Zero),
			    zeroGradientFvPatchField<scalar>::typeName
			);
			
			fvScalarMatrix psilEqn
			(
				fvm::Sp(scalar(1),psil) - fvm::laplacian(DPsi,psil) == mCondAlphal_
			);
			
			psilEqn.solve();
			
			//- Cut cells with cutoff < alpha1 < 1-cutoff and rescale remaining source term field
			dimensionedScalar intPsiLiquidCondensation ("intPsiLiquidCondensation", dimensionSet(1,0,-1,0,0,0,0), 0.0);
			dimensionedScalar intPsiVaporCondensation ("intPsiVaporCondensation", dimensionSet(1,0,-1,0,0,0,0), 0.0);


			forAll(mesh.C(),iCell)
			{
				if ((limitedAlphal_[iCell]) < cutoff_)
				{
					intPsiVaporCondensation.value() += (1.0-limitedAlphal_[iCell])*psil[iCell]*mesh.V()[iCell];
				}
				else if ((limitedAlphal_[iCell]) > 1.0-cutoff_)
				{
					intPsiLiquidCondensation.value() += (limitedAlphal_[iCell])*psil[iCell]*mesh.V()[iCell];
				}
			}
			
			//- Calculate Nl and Nv
			dimensionedScalar Nl ("Nl", dimensionSet(0,0,0,0,0,0,0), 2.0);
			dimensionedScalar Nv ("Nv", dimensionSet(0,0,0,0,0,0,0), 2.0);
			
			reduce(intPsiLiquidCondensation.value(),sumOp<scalar>());
			reduce(intPsiVaporCondensation.value(),sumOp<scalar>());
			
			if (intPsiLiquidCondensation.value() > 1e-99)
			{
			    Nl = intPsi0l/intPsiLiquidCondensation;
			}
			if (intPsiVaporCondensation.value() > 1e-99)
			{
			    Nv = intPsi0l/intPsiVaporCondensation;
			}
			
			        
			//- Set source terms in cells with alpha1 < cutoff or alpha1 > 1-cutoff
			forAll(mesh.C(),iCell)
			{
				if (limitedAlphal_[iCell] < cutoff_)
				{
					mCondAlphal_[iCell] = -Nv.value()*(1.0-limitedAlphal_[iCell])*psil[iCell];
				}
				else if (limitedAlphal_[iCell] > 1.0-cutoff_)
				{
					mCondAlphal_[iCell] = Nl.value()*(limitedAlphal_[iCell])*psil[iCell];
				}
				else
				{
					mCondAlphal_[iCell] = 0.0;
				}
			}
			
			//- Evaporation source term in energy equation
			//hESource = -N*alpha1*psi0Tild/Rph;
		}

		if (evap_)
		{
			mEvapNoAlphal_ *= N;
			//mEvapAlphal_   *= N;
			mEvapNoTmTSat_ *= N;

			dimensionedScalar intPsi0v = fvc::domainIntegrate(mEvapAlphal_);

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
			
			fvScalarMatrix psivEqn
			(
				fvm::Sp(scalar(1),psiv) - fvm::laplacian(DPsi,psiv) == mEvapAlphal_
			);
			
			psivEqn.solve();
			
			//- Cut cells with cutoff < alpha1 < 1-cutoff and rescale remaining source term field
			dimensionedScalar intPsiLiquidEvaporation ("intPsiLiquidEvaporation", dimensionSet(1,0,-1,0,0,0,0), 0.0);
			dimensionedScalar intPsiVaporEvaporation ("intPsiVaporEvaporation", dimensionSet(1,0,-1,0,0,0,0), 0.0);


			forAll(mesh.C(),iCell)
			{
				if (limitedAlphal_[iCell] < cutoff_)
				{
					intPsiVaporEvaporation.value() += (1.0-limitedAlphal_[iCell])*psiv[iCell]*mesh.V()[iCell];
				}
				else if (limitedAlphal_[iCell] > 1.0-cutoff_)
				{
					intPsiLiquidEvaporation.value() += limitedAlphal_[iCell]*psiv[iCell]*mesh.V()[iCell];
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
				if (limitedAlphal_[iCell] < cutoff_)
				{
					mEvapAlphal_[iCell] = -Nv.value()*(1.0-limitedAlphal_[iCell])*psiv[iCell];
				}
				else if (limitedAlphal_[iCell] > 1.0-cutoff_)
				{
					mEvapAlphal_[iCell] = Nl.value()*limitedAlphal_[iCell]*psiv[iCell];
				}
				else
				{
					mEvapAlphal_[iCell] = 0.0;
				}
			}
			
			//- Evaporation source term in energy equation
			//hESource = -N*alpha1*psi0Tild/Rph;
			//correctionTerm_ = psiv*(-Nv*(1.0-limitedAlphal_)*cp2_+Nl*limitedAlphal_*cp1_);

			///////////////////////////////////////////////////////////////////////////////

			//dimensionedScalar intPsiNoAlphal0v = fvc::domainIntegrate(mEvapNoAlphal_);

			//volScalarField psivNoAlphal
			//(
			//    IOobject
			//    (
			//        "psivNoAlphal",
			//        mesh.time().timeName(),
			//        mesh,
			//        IOobject::NO_READ,
			//        IOobject::NO_WRITE
			//    ),
			//    mesh,
			//    dimensionedScalar(mEvapAlphal_.dimensions(), Zero),
			//    zeroGradientFvPatchField<scalar>::typeName
			//);
			//
			//fvScalarMatrix psivNoAlphalEqn
			//(
			//	fvm::Sp(scalar(1),psivNoAlphal) - fvm::laplacian(DPsi,psivNoAlphal) == mEvapNoAlphal_
			//);
			//
			//psivNoAlphalEqn.solve();
			//
			////- Cut cells with cutoff < alpha1 < 1-cutoff and rescale remaining source term field
			//dimensionedScalar intPsiNoAlphalLiquidEvaporation ("intPsiNoAlphalLiquidEvaporation", dimensionSet(1,0,-1,0,0,0,0), 0.0);
			//dimensionedScalar intPsiNoAlphalVaporEvaporation ("intPsiNoAlphalVaporEvaporation", dimensionSet(1,0,-1,0,0,0,0), 0.0);


			//forAll(mesh.C(),iCell)
			//{
			//	if (limitedAlphal_[iCell] < cutoff_)
			//	{
			//		intPsiNoAlphalVaporEvaporation.value() += (1.0-limitedAlphal_[iCell])*psivNoAlphal[iCell]*mesh.V()[iCell];
			//	}
			//	else if (limitedAlphal_[iCell] > 1.0-cutoff_)
			//	{
			//		intPsiNoAlphalLiquidEvaporation.value() += limitedAlphal_[iCell]*psivNoAlphal[iCell]*mesh.V()[iCell];
			//	}
			//}
			//
			////- Calculate Nl and Nv
			//dimensionedScalar NlNoAlphal ("NlNoAlphal", dimensionSet(0,0,0,0,0,0,0), 2.0);
			//dimensionedScalar NvNoAlphal ("NvNoAlphal", dimensionSet(0,0,0,0,0,0,0), 2.0);
			//
			//reduce(intPsiNoAlphalLiquidEvaporation.value(),sumOp<scalar>());
			//reduce(intPsiNoAlphalVaporEvaporation.value(),sumOp<scalar>());
			//
			//if (intPsiNoAlphalLiquidEvaporation.value() > 1e-99)
			//{
			//    NlNoAlphal = intPsiNoAlphal0v/intPsiNoAlphalLiquidEvaporation;
			//}
			//if (intPsiNoAlphalVaporEvaporation.value() > 1e-99)
			//{
			//    NvNoAlphal = intPsiNoAlphal0v/intPsiNoAlphalVaporEvaporation;
			//}
			//
			//        
			////- Set source terms in cells with alpha1 < cutoff or alpha1 > 1-cutoff
			//forAll(mesh.C(),iCell)
			//{
			//	if (limitedAlphal_[iCell] < cutoff_)
			//	{
			//		mEvapNoAlphal_[iCell] = -NvNoAlphal.value()*(1.0-limitedAlphal_[iCell])*psivNoAlphal[iCell];
			//	}
			//	else if (limitedAlphal_[iCell] > 1.0-cutoff_)
			//	{
			//		mEvapNoAlphal_[iCell] = NlNoAlphal.value()*limitedAlphal_[iCell]*psivNoAlphal[iCell];
			//	}
			//	else
			//	{
			//		mEvapNoAlphal_[iCell] = 0.0;
			//	}
			//}
		}

		magGradLimitedAlphalCalculated_ = false;
		limitedAlphalCalculated_ = false;
	}
}

void Foam::phaseChangeTwoPhaseMixture::calcTSatLocal() 
{
	if (TSatLocalPressure_)
	{
		//Info <<"TSat is calculated based on local pressure field." << endl;
	    TSat_ = 1.0/(1.0/TSatG_ - R_/hEvap_*log(max(p_/pSat_,1E-08)));
	}
	//else
	//{
	//	Info <<"TSat is constant, TSat = " << TSatG_ << endl;
	//}
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixture::vDotAlphal() const
{
    volScalarField alphalCoeff(1.0/rho1() - alpha1()*(1.0/rho1() - 1.0/rho2()));
    Pair<tmp<volScalarField> > mDotAlphal = this->mDotAlphal();

    return Pair<tmp<volScalarField> >
    (
        alphalCoeff*mDotAlphal[0],
        alphalCoeff*mDotAlphal[1]
    );
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixture::vDotP() const
{
    dimensionedScalar pCoeff(1.0/rho1() - 1.0/rho2());
    Pair<tmp<volScalarField> > mDotP = this->mDotP();

    return Pair<tmp<volScalarField> >
	(
	    pCoeff*mDotP[0], 
		pCoeff*mDotP[1]
	);
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixture::vDotT() const
{
	Pair<tmp<volScalarField> > mDotT = this->mDotT();

	return Pair<tmp<volScalarField> >
	(
			hEvap_*mDotT[0],
			hEvap_*mDotT[1]
	);
}

void Foam::phaseChangeTwoPhaseMixture::correct()
{
	thermalIncompressibleTwoPhaseMixture::correct();
	calcTSatLocal();


    const fvMesh& mesh = alpha1().mesh();

	if (printPhaseChange_)
	{
    	Info<< "****Condensation rate: "
    	    << gSum((mCondAlphal_*mesh.V())())*hEvap_.value() << " W" << endl;
    	Info<< "****Evaporation rate: "
    	    << gSum((mEvapAlphal_*mesh.V())())*hEvap_.value() << " W" << endl;
	}
}

bool Foam::phaseChangeTwoPhaseMixture::read()
{
    if (incompressibleTwoPhaseMixture::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");
        lookup("TSatGlobal") >> TSatG_;
        lookup("TSatLocalPressure") >> TSatLocalPressure_;
        lookup("pSat") >> pSat_;
        lookup("hEvap") >> hEvap_;
        lookup("R") >> R_;
        lookup("printPhaseChange") >> printPhaseChange_;
        lookup("HardtWondra") >> HW_;
        lookup("cutoff") >> cutoff_;
        lookup("condensation") >> cond_;
        lookup("evaporation") >> evap_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
