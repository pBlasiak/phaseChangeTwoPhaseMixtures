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
	isHW_(phaseChangeTwoPhaseMixtureCoeffs_.getOrDefault<Switch>("HardtWondra", true)),
    satProps_
    (
        SaturationProperties::New
        (
            U,
            phi, 
			phaseChangeTwoPhaseMixtureCoeffs_.get<word>("satPropModel")
        )
    ),
    HW_(new HardtWondra(alpha1(), satProps_.ref())),
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
	printPhaseChange_(readBool(phaseChangeTwoPhaseMixtureCoeffs_.lookup("printPhaseChange")))
{
	Info<< "printPhaseChange = "         << printPhaseChange_ << endl;
	Info<< "Hardt-Wondra algorithm is: " << isHW_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
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
        lookup("spread") >> spread_;
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
