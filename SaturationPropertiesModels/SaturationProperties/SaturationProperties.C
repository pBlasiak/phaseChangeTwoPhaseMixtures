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

#include "SaturationProperties.H"
//#include "thermalIncompressibleTwoPhaseMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SaturationProperties, 0);
    defineRunTimeSelectionTable(SaturationProperties, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SaturationProperties::SaturationProperties
(
    const word& type,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    SaturationPropertiesDict_
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
    U_(U),
    phi_(phi),
	p_(U.db().lookupObject<volScalarField>("p")),
	T_(U.db().lookupObject<volScalarField>("T")),
    TSatG_("TSatGlobal", dimTemperature, SaturationPropertiesDict_),
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
    pSat_("pSat", dimPressure, SaturationPropertiesDict_),
    hEvap_("hEvap", dimEnergy/dimMass, SaturationPropertiesDict_)
{
	Info<< "TSatGlobal = "				 << TSatG_ << endl;
	Info<< "pSat = "		  			 << pSat_ << endl;
	Info<< "hEvap = "		  			 << hEvap_ << endl; 
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //




// ************************************************************************* //
