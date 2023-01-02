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
    Foam::SaturationPropertiesModels::Constant

Description
    Constant thermal property model for a mixture. 
	thermProp = alpha1*thermProp1 + alpha2*thermProp2

SourceFiles
    Constant.C

\*--------------------------------------------------------------------*/

#ifndef Constant_H
#define Constant_H

#include "SaturationProperties.H"
#include "thermalIncompressibleTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace SaturationPropertiesModels
{

/*--------------------------------------------------------------------*\
                              Class Constant
\*--------------------------------------------------------------------*/

class Constant
:
    public SaturationProperties
{
    // Private data

		

public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- construct from components
        Constant
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~Constant()
    {}


    // Member Functions
	
        //- Return the effective laminar thermal property
        virtual void calcTSat () override; 

        //- Read the transportProperties dictionary and update
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SaturationPropertiesModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
