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
    Foam::thermalProperty

Description

SourceFiles
    thermalProperty.C
    newThermalProperty.C

\*---------------------------------------------------------------------------*/

#ifndef thermalProperty_H
#define thermalProperty_H

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
                Class thermalProperty Declaration
\*---------------------------------------------------------------------------*/

class thermalProperty
{

protected:

    // Protected data

        dictionary thermalPropertyDict_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;

    // Private Member Functions

        //- Disallow copy construct
        thermalProperty(const thermalProperty&);

        //- Disallow default bitwise assignment
        void operator=(const thermalProperty&);

public:

    //- Runtime type information
    TypeName("thermalProperty");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            thermalProperty,
            components,
            (
                const volVectorField& U,
                const surfaceScalarField& phi
            ),
            (U, phi)
        );


    // Selectors

        //- Return a reference to the selected thermal conductivity model
        static autoPtr<thermalProperty> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
			const word& thermProp
        );


    // Constructors

        //- Construct from components
        thermalProperty
        (
            const word& type,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~thermalProperty()
    {}


    // Member Functions

	typedef  const dimensionedScalar& (Foam::thermalIncompressibleTwoPhaseMixture::*th)() const;

        //- Return the effective laminar thermal conductivity
        virtual tmp<volScalarField> calcThermProp
		(
			const thermalIncompressibleTwoPhaseMixture* titpm,
			th T1 = nullptr, 
			th T2 = nullptr
			//const dimensionedScalar& T1,
			//const dimensionedScalar& T2
		) const = 0;
		
        //- Return the effective laminar thermal conductivity
        virtual tmp<volScalarField> calcThermProp
		(
			const thermalIncompressibleTwoPhaseMixture* titpm,
			const volScalarField& T1,
			const volScalarField& T2
		) const = 0;

        //- Read the transportProperties dictionary and update
        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
