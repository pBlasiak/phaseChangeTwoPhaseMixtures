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

#include "thermalConductivity.H"
//#include "thermalIncompressibleTwoPhaseMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalConductivity, 0);
    defineRunTimeSelectionTable(thermalConductivity, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalConductivity::thermalConductivity
(
    const word& type,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    thermalConductivityDict_
	(
	    IOdictionary
	    (
	        IOobject
            (
                "transportProperties",   // dictionary name
                U.time().constant(),     // dict is found in "constant"
                U.db(),                  // registry for the dict
                IOobject::MUST_READ,     // must exist, otherwise failure
                IOobject::NO_WRITE       // dict is only read by the solver
            )
	    )
	),
    U_(U),
    phi_(phi)
{ }


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


bool Foam::thermalConductivity::read()
{
    //if (incompressibleTwoPhaseMixture::read())
    //{
    //    thermalConductivityCoeffs_ = subDict(type() + "Coeffs");
    //    lookup("TSatGlobal") >> TSatG_;
    //    lookup("TSatLocalPressure") >> TSatLocalPressure_;
    //    lookup("pSat") >> pSat_;
    //    lookup("hEvap") >> hEvap_;
    //    lookup("R") >> R_;
    //    lookup("printPhaseChange") >> printPhaseChange_;

        return true;
    //}
    //else
    //{
    //    return false;
    //}
}


// ************************************************************************* //
