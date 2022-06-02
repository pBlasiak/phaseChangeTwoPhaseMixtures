/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::thermalPropertyModels::Step

Description
    Step thermal property model for a mixture. 
	thermProp = if alpha1 >= 0.5 return thermProp1 
	            else return thermProp2

SourceFiles
    Step.C

\*--------------------------------------------------------------------*/

#ifndef Step_H
#define Step_H

#include "thermalProperty.H"
#include "thermalIncompressibleTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPropertyModels
{

/*--------------------------------------------------------------------*\
                              Class Step
\*--------------------------------------------------------------------*/

class Step
:
    public thermalProperty
{
    // Private data

        //dimensionedScalar Cc_;
        //dimensionedScalar Cv_;

public:

    //- Runtime type information
    TypeName("step");


    // Constructors

        //- construct from components
        Step
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~Step()
    {}


    // Member Functions
	
        //- Return the effective laminar thermal property
        virtual tmp<volScalarField> calcThermProp
		(
			const thermalIncompressibleTwoPhaseMixture* titpm,
			th T1 = nullptr, 
			th T2 = nullptr
		) const override;

        virtual tmp<volScalarField> calcThermProp
		(
			const thermalIncompressibleTwoPhaseMixture* titpm,
			const volScalarField& T1,
			const volScalarField& T2
		) const override;

        //- Read the transportProperties dictionary and update
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalPropertyModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
