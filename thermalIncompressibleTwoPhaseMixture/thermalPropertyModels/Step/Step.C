/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "Step.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPropertyModels
{
    defineTypeNameAndDebug(Step, 0);
    addToRunTimeSelectionTable(thermalProperty, Step, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalPropertyModels::Step::Step
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    thermalProperty(typeName, U, phi)

    //Cc_(thermalPropertyCoeffs_.subDict(type() + "Coeffs").lookup("Cc")),
    //Cv_(thermalPropertyCoeffs_.subDict(type() + "Coeffs").lookup("Cv")),

    //mcCoeff_(Cc_*rho2()),
    //mvCoeff_(Cv_*rho1())
{
	//Info<< "Phase change relaxation time factors for the Step model:\n" 
	//	<< "Cc = " << Cc_ << endl
	//	<< "Cv = " << Cv_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> 
Foam::thermalPropertyModels::Step::calcThermProp
(
	const thermalIncompressibleTwoPhaseMixture* titpm,
	th T1,
	th T2 
) const
{
	const volScalarField limitedAlpha1
	(
		min(max(titpm->alpha1(), scalar(0)), scalar(1))
	);

	const dimensionedScalar thpr1 = (*titpm.*T1)();
	const dimensionedScalar thpr2 = (*titpm.*T2)();

	volScalarField thPr
	(
        IOobject
        (
            "thPr",
            U_.time().timeName(),
            U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U_.mesh(),
		thpr1
	);

	forAll(limitedAlpha1, celli)
	{
		if (limitedAlpha1[celli] < 0.5)
		{
			thPr[celli] = thpr2.value(); 
		}
	}

    return tmp<volScalarField>
    (
		new volScalarField(thPr)
	);
}

Foam::tmp<Foam::volScalarField> 
Foam::thermalPropertyModels::Step::calcThermProp
(
	const thermalIncompressibleTwoPhaseMixture* titpm,
	const volScalarField& T1,
	const volScalarField& T2
) const
{
	const volScalarField limitedAlpha1
	(
		min(max(titpm->alpha1(), scalar(0)), scalar(1))
	);

	volScalarField thPr
	(
        IOobject
        (
            "thPr",
            U_.time().timeName(),
            U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
		T1
	);

	forAll(limitedAlpha1, celli)
	{
		if (limitedAlpha1[celli] < 0.5)
		{
			thPr[celli] = T2[celli];
		}
	}

    return tmp<volScalarField>
    (
		new volScalarField(thPr)
	);
}

bool Foam::thermalPropertyModels::Step::read()
{
    //if (thermalProperty::read())
    //{
        //thermalPropertyCoeffs_ = subDict(type() + "Coeffs");

        //thermalPropertyCoeffs_.lookup("Cc") >> Cc_;
        //thermalPropertyCoeffs_.lookup("Cv") >> Cv_;

        //mcCoeff_ = Cc_*rho2();
        //mvCoeff_ = Cv_*rho1();

        return true;
    //}
    //else
    //{
    //    return false;
    //}
}


// ************************************************************************* //
