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

#include "HarmonicNu.H"
#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kinematicViscosityModels
{
    defineTypeNameAndDebug(HarmonicNu, 0);
    addToRunTimeSelectionTable(kinematicViscosity, HarmonicNu, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kinematicViscosityModels::HarmonicNu::HarmonicNu
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    kinematicViscosity(typeName, U, phi)

{ }


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::kinematicViscosityModels::HarmonicNu::correct
(
	Foam::thermalIncompressibleTwoPhaseMixture* titpm
)
{
    viscosityModel& nuModel1 = const_cast<viscosityModel&>(titpm->nuModel1());
    viscosityModel& nuModel2 = const_cast<viscosityModel&>(titpm->nuModel2());
    nuModel1.correct();
    nuModel2.correct();

    const volScalarField limitedAlpha1
    (
        "limitedAlpha1",
        min(max(titpm->alpha1(), scalar(0)), scalar(1))
    );

	const dimensionedScalar rho1 = titpm->rho1();
	const dimensionedScalar rho2 = titpm->rho2();
	const volScalarField nu1 = nuModel1.nu();
	const volScalarField nu2 = nuModel2.nu();

    titpm->setNu( 
		scalar(1.0)/
		(
			(1.0/(rho1*nu1) - 1.0/(rho2*nu2))*limitedAlpha1
		  + 1.0/(rho2*nu2)
		)/(limitedAlpha1*rho1 + (scalar(1) - limitedAlpha1)*rho2)
		);
}

bool Foam::kinematicViscosityModels::HarmonicNu::read()
{
    //if (kinematicViscosity::read())
    //{
        //kinematicViscosityCoeffs_ = subDict(type() + "Coeffs");

        //kinematicViscosityCoeffs_.lookup("Cc") >> Cc_;
        //kinematicViscosityCoeffs_.lookup("Cv") >> Cv_;

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
