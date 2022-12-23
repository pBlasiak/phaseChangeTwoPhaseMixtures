/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 AUTHOR,AFFILIATION
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

#include "HardtWondra.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//const dataType Foam::HardtWondra::staticData();


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HardtWondra::HardtWondra
(
	const volScalarField& alpha1
)
:
	HWdict_
	(
	    IOdictionary
	    (
	        IOobject
            (
                "phaseChangeProperties", // dictionary name
                alpha1.time().constant(),     // dict is found in "constant"
                alpha1.db(),                  // registry for the dict
                IOobject::MUST_READ,     // must exist, otherwise failure
                IOobject::NO_WRITE       // dict is only read by the solver
            )
	    )
	),
	cutoff_(HWdict_.getOrDefault<scalar>("cutoff", 1e-3)),
	spread_(HWdict_.getOrDefault<scalar>("spread", 3))
{
	Info<< "Spread is set as: "          << spread_ << endl;
	Info<< "Cutoff is set as: "          << cutoff_ << endl;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

//Foam::autoPtr<Foam::HardtWondra>
//Foam::HardtWondra::New()
//{
//    return autoPtr<HardtWondra>(new HardtWondra);
//}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::HardtWondra::~HardtWondra()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
