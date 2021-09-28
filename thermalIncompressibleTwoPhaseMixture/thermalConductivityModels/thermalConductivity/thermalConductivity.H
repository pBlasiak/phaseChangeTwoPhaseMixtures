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

Class
    Foam::thermalConductivity

Description

SourceFiles
    thermalConductivity.C
    newPhaseChangeModel.C

\*---------------------------------------------------------------------------*/

#ifndef thermalConductivity_H
#define thermalConductivity_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class thermalIncompressibleTwoPhaseMixture;

/*---------------------------------------------------------------------------*\
                Class thermalConductivity Declaration
\*---------------------------------------------------------------------------*/

class thermalConductivity
{

protected:

    // Protected data

        dictionary thermalConductivityDict_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;

    // Private Member Functions

        //- Disallow copy construct
        thermalConductivity(const thermalConductivity&);

        //- Disallow default bitwise assignment
        void operator=(const thermalConductivity&);

public:

    //- Runtime type information
    TypeName("thermalConductivity");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            thermalConductivity,
            components,
            (
                const volVectorField& U,
                const surfaceScalarField& phi
            ),
            (U, phi)
        );


    // Selectors

        //- Return a reference to the selected thermal conductivity model
        static autoPtr<thermalConductivity> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Constructors

        //- Construct from components
        thermalConductivity
        (
            const word& type,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~thermalConductivity()
    {}


    // Member Functions


        //- Return the effective laminar thermal conductivity
        virtual tmp<volScalarField> k(const thermalIncompressibleTwoPhaseMixture* titpm) const = 0;

        //- Read the transportProperties dictionary and update
        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
