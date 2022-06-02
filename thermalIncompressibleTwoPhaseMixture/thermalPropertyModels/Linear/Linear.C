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

#include "Linear.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPropertyModels
{
    defineTypeNameAndDebug(Linear, 0);
    addToRunTimeSelectionTable(thermalProperty, Linear, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalPropertyModels::Linear::Linear
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
	//Info<< "Phase change relaxation time factors for the Linear model:\n" 
	//	<< "Cc = " << Cc_ << endl
	//	<< "Cv = " << Cv_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> 
Foam::thermalPropertyModels::Linear::calcThermProp
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

    return tmp<volScalarField>
    (
		new volScalarField
        (
            "linearThermProp",
            limitedAlpha1*thpr1
          + (scalar(1) - limitedAlpha1)*thpr2
        )
	);
}

Foam::tmp<Foam::volScalarField> 
Foam::thermalPropertyModels::Linear::calcThermProp
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

    return tmp<volScalarField>
    (
		new volScalarField
        (
            "linearThermProp",
            limitedAlpha1*T1
          + (scalar(1) - limitedAlpha1)*T2
        )
	);
}

bool Foam::thermalPropertyModels::Linear::read()
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
