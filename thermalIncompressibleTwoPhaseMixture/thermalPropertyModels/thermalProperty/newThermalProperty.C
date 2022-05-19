/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "thermalProperty.H"
//#include "incompressibleTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::thermalProperty>
Foam::thermalProperty::New
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word& thermProp
)
{
    IOdictionary transportPropertiesDict
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    word thermalPropertyTypeName
    (
        transportPropertiesDict.lookup(thermProp + "Model")
    );

    Info<< "Selecting "
        << thermalPropertyTypeName << endl;

    componentsConstructorTable::iterator cstrIter =
        componentsConstructorTablePtr_
            ->find(thermalPropertyTypeName);

    if (cstrIter == componentsConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown thermalProperty type "
            << thermalPropertyTypeName << endl << endl
            << "Valid  thermalPropertys are : " << endl
            << componentsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<thermalProperty>(cstrIter()(U, phi));
}


// ************************************************************************* //
