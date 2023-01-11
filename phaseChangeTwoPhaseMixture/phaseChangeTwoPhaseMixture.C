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
    cond_{phaseChangeTwoPhaseMixtureCoeffs_.get<Switch>("condensation")},
    evap_{phaseChangeTwoPhaseMixtureCoeffs_.get<Switch>("evaporation")},
	isHW_(phaseChangeTwoPhaseMixtureCoeffs_.getOrDefault<Switch>("HardtWondra", true)),
    //satProps_(autoPtr<SaturationProperties>()),
    satProps_
    (
        SaturationProperties::New
        (
            U,
            phi
        )
    ),
    HW_(autoPtr<HardtWondra>()),
    jc_
    (
        IOobject
        (
            "jc",
            U.time().timeName(),
            U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("jc", dimensionSet(1, -2, -1, 0, 0, 0, 0), 0.0)
    ),
    je_
    (
        IOobject
        (
            "je",
            U.time().timeName(),
            U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("je", dimensionSet(1, -2, -1, 0, 0, 0, 0), 0.0)
    ),
	//TODO: nie wiem czy jednostki OK
	alphaSourceSp_
	(
	    IOobject
	    (
	        "alphaSourceSp",
	        U.time().timeName(),
	        U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
	    ),
	    U.mesh(),
	    dimensionedScalar("alphaSourceSp", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
	),
	alphaSourceSu_
	(
	    IOobject
	    (
	        "alphaSourceSu",
	        U.time().timeName(),
	        U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
	    ),
	    U.mesh(),
	    dimensionedScalar("alphaSourceSu", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
	),
	pSourceSp_
	(
	    IOobject
	    (
	        "pSourceSp",
	        U.time().timeName(),
	        U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
	    ),
	    U.mesh(),
	    dimensionedScalar("pSourceSp", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
	),
	pSourceSu_
	(
	    IOobject
	    (
	        "pSourceSu",
	        U.time().timeName(),
	        U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
	    ),
	    U.mesh(),
	    dimensionedScalar("pSourceSu", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
	),
	TSourceSp_
	(
	    IOobject
	    (
	        "TSourceSp",
	        U.time().timeName(),
	        U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
	    ),
	    U.mesh(),
	    dimensionedScalar("TSourceSp", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
	),
	TSourceSu_
	(
	    IOobject
	    (
	        "TSourceSu",
	        U.time().timeName(),
	        U.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
	    ),
	    U.mesh(),
	    dimensionedScalar("TSourceSu", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
	),
	//mCondNoTmTSat_
	//(
	//    IOobject
	//    (
	//        "mCondNoTmTSat",
	//        U.time().timeName(),
	//        U.db(),
	//		IOobject::NO_READ,
	//		IOobject::NO_WRITE
	//    ),
	//    U.mesh(),
	//    dimensionedScalar("mCondNoTmTSat", dimensionSet(1, -3, -1, -1, 0, 0, 0), 0.0)
	//),
	//mEvapNoTmTSat_
	//(
	//    IOobject
	//    (
	//        "mEvapNoTmTSat",
	//        U.time().timeName(),
	//        U.db(),
	//		IOobject::NO_READ,
	//		IOobject::NO_WRITE
	//    ),
	//    U.mesh(),
	//    dimensionedScalar("mEvapNoTmTSat", dimensionSet(1, -3, -1, -1, 0, 0, 0), 0.0)
	//),
	//mCondP_
	//(
	//    IOobject
	//    (
	//        "mCondP",
	//        U.time().timeName(),
	//        U.db(),
	//		IOobject::NO_READ,
	//		IOobject::NO_WRITE
	//    ),
	//    U.mesh(),
	//    dimensionedScalar("mCondP", dimensionSet(0, -2, 1, 0, 0, 0, 0), 0.0)
	//),
	//mEvapP_
	//(
	//    IOobject
	//    (
	//        "mEvapP",
	//        U.time().timeName(),
	//        U.db(),
	//		IOobject::NO_READ,
	//		IOobject::NO_WRITE
	//    ),
	//    U.mesh(),
	//    dimensionedScalar("mEvapP", dimensionSet(0, -2, 1, 0, 0, 0, 0), 0.0)
	//),
	//mCondT_
	//(
	//    IOobject
	//    (
	//        "mCondT",
	//        U.time().timeName(),
	//        U.db(),
	//		IOobject::NO_READ,
	//		IOobject::NO_WRITE
	//    ),
	//    U.mesh(),
	//    dimensionedScalar("mCondT", dimensionSet(1, -3, -1, -1, 0, 0, 0), 0.0)
	//),
	//mEvapT_
	//(
	//    IOobject
	//    (
	//        "mEvapT",
	//        U.time().timeName(),
	//        U.db(),
	//		IOobject::NO_READ,
	//		IOobject::NO_WRITE
	//    ),
	//    U.mesh(),
	//    dimensionedScalar("mEvapT", dimensionSet(1, -3, -1, -1, 0, 0, 0), 0.0)
	//),
	printPhaseChange_(readBool(phaseChangeTwoPhaseMixtureCoeffs_.lookup("printPhaseChange")))
{
    //satProps_.reset
    //(
    //    SaturationProperties::New
    //    (
    //        U,
    //        phi
    //    )
    //),
	HW_.reset
	(
		new HardtWondra(alpha1(), satProps_.ref(), *this)
	);
	Info<< "Condensation is   "   << cond_   << endl;
	Info<< "Evaporation is    "   << evap_   << endl;
	Info<< "Hardt-Wondra algorithm is: " << isHW_ << endl;
	Info<< "printPhaseChange = "         << printPhaseChange_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixture::alphaSource() 
{
    volScalarField alphalCoeff(1.0/rho1() - alpha1()*(1.0/rho1() - 1.0/rho2()));
    tmp<volScalarField> Sp = this->alphaSourceSp();
    tmp<volScalarField> Su = this->alphaSourceSu();

    return Pair<tmp<volScalarField>>
    (
        alphalCoeff*Sp,
        alphalCoeff*Su
    );
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixture::pSource() 
{
    dimensionedScalar pCoeff(1.0/rho1() - 1.0/rho2());
    tmp<volScalarField> Sp = this->pSourceSp();
    tmp<volScalarField> Su = this->pSourceSu();

    return Pair<tmp<volScalarField>>
	(
	    pCoeff*Sp, 
		pCoeff*Su
	);
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixture::TSource()
{
	tmp<volScalarField> Sp = this->TSourceSp();
	tmp<volScalarField> Su = this->TSourceSu();

	return Pair<tmp<volScalarField>>
	(
			satProps_->hEvap()*Sp,
			satProps_->hEvap()*Su
	);
}

void Foam::phaseChangeTwoPhaseMixture::correct()
{
	thermalIncompressibleTwoPhaseMixture::correct();

	satProps_->calcTSat();

	HW_->correct(); 

    const fvMesh& mesh = alpha1().mesh();

	if (printPhaseChange_)
	{
		//TODO: zastanowic sie jak to liczyc i czy jest potrzebne
    	//Info<< "****Condensation rate: "
    	//    << gSum((mCondAlphal_*mesh.V())())*satProps_->hEvap().value() << " W" << endl;
    	//Info<< "****Evaporation rate: "
    	//    << gSum((mEvapAlphal_*mesh.V())())*satProps_->hEvap().value() << " W" << endl;
	}
}

bool Foam::phaseChangeTwoPhaseMixture::read()
{
    if (incompressibleTwoPhaseMixture::read())
    {
		//TODO: cos tu jest raczej zle
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");
        lookup("HardtWondra") >> isHW_;
        lookup("printPhaseChange") >> printPhaseChange_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
