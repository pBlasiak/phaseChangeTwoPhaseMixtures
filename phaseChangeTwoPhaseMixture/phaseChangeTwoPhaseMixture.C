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
    TSatG_("TSatGlobal", dimTemperature, phaseChangeTwoPhaseMixtureCoeffs_.lookup("TSatGlobal")),
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
    pSat_("pSat", dimPressure, phaseChangeTwoPhaseMixtureCoeffs_.lookup("pSat")),
    hEvap_("hEvap", dimEnergy/dimMass, phaseChangeTwoPhaseMixtureCoeffs_.lookup("hEvap")),
    R_("R", dimGasConstant, phaseChangeTwoPhaseMixtureCoeffs_.lookup("R")),
    TSatLocalPressure_(readBool(phaseChangeTwoPhaseMixtureCoeffs_.lookup("TSatLocalPressure"))),
	printPhaseChange_(readBool(phaseChangeTwoPhaseMixtureCoeffs_.lookup("printPhaseChange")))
{
	Info<< "TSatGlobal = "				<< TSatG_ << endl;
	Info<< "pSat = "		  			<< pSat_ << endl;
	Info<< "hEvap = "		  			<< hEvap_ << endl;
	Info<< "R = "			  			<< R_ << endl;
	Info<< "TSatLocalPressure = "       << TSatLocalPressure_ << endl;
	Info<< "printPhaseChange = "        << printPhaseChange_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

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
Foam::phaseChangeTwoPhaseMixture::vDotAlphal() 
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

    volScalarField limitedAlpha1 = min(max(alpha1(), scalar(0)), scalar(1));
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

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
