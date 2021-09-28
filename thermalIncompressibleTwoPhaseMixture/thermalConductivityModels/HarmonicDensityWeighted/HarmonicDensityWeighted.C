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

#include "HarmonicDensityWeighted.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalConductivityModels
{
    defineTypeNameAndDebug(HarmonicDensityWeighted, 0);
    addToRunTimeSelectionTable(thermalConductivity, HarmonicDensityWeighted, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalConductivityModels::HarmonicDensityWeighted::HarmonicDensityWeighted
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    thermalConductivity(typeName, U, phi)

    //Cc_(thermalConductivityCoeffs_.subDict(type() + "Coeffs").lookup("Cc")),
    //Cv_(thermalConductivityCoeffs_.subDict(type() + "Coeffs").lookup("Cv")),

    //mcCoeff_(Cc_*rho2()),
    //mvCoeff_(Cv_*rho1())
{
	//Info<< "Phase change relaxation time factors for the HarmonicDensityWeighted model:\n" 
	//	<< "Cc = " << Cc_ << endl
	//	<< "Cv = " << Cv_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> 
Foam::thermalConductivityModels::HarmonicDensityWeighted::k(const Foam::thermalIncompressibleTwoPhaseMixture* titpm) const
{
	const volScalarField limitedAlpha1
	(
		min(max(titpm->alpha1(), scalar(0)), scalar(1))
	);

    return tmp<volScalarField>
    (
		new volScalarField
        (
            "kHarmonicDensityWeighted",
			(
				titpm->k1()*limitedAlpha1/titpm->rho1() 
			  - (scalar(1.0) - limitedAlpha1)*titpm->k2()/titpm->rho2() 
			)/(limitedAlpha1/titpm->rho1() - (scalar(1.0) - limitedAlpha1)/titpm->rho2())
        )
	);
}

bool Foam::thermalConductivityModels::HarmonicDensityWeighted::read()
{
    //if (thermalConductivity::read())
    //{
        //thermalConductivityCoeffs_ = subDict(type() + "Coeffs");

        //thermalConductivityCoeffs_.lookup("Cc") >> Cc_;
        //thermalConductivityCoeffs_.lookup("Cv") >> Cv_;

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
