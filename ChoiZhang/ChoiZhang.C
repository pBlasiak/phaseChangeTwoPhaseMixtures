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

#include "ChoiZhang.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(ChoiZhang, 0);
    addToRunTimeSelectionTable(phaseChangeTwoPhaseMixture, ChoiZhang, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixtures::ChoiZhang::ChoiZhang
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    Lee(U, phi),

    Ccvar_(Cc_)
    //Cv_(phaseChangeTwoPhaseMixtureCoeffs_.subDict(type() + "Coeffs").lookup("Cv")),

{
	Info<< "Phase change relaxation time factors for the " << type() << " model:\n" 
		<< "Cc = " << Ccvar_ << endl
		<< "Cv = " << Cv_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::phaseChangeTwoPhaseMixtures::ChoiZhang::correct()
{
	phaseChangeTwoPhaseMixture::correct();

	//Ccvar_ = Cc_; // jak tego nie ma to Cc najpierw maleje do 1e-150 apotem rosnie do 1e+300
	scalar totvmCond(0);
	scalar totvmEvap(0);
	totvmCond = gSum(mCondAlphal_);
	totvmEvap = gSum(mEvapAlphal_);
	if (!(totvmCond == 0 && mag(totvmEvap) == 0))
		Ccvar_.value() -= Ccvar_.value()*((totvmCond - mag(totvmEvap))
	    	/max(totvmCond, mag(totvmEvap)));
	
    mcCoeff_ = Ccvar_*rho2();
	//Info<< "Updated mcCoeff = " << mcCoeff_.value() << endl;

}


// ************************************************************************* //
