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
    Foam::kinematicViscosityModels::LinearNu

Description
TODO poprawic
    LinearNu thermal conductivity model for a mixture. 
	kEff = alpha1*k1 + alpha2*k2

SourceFiles
    LinearNu.C

\*--------------------------------------------------------------------*/

#ifndef LinearNu_H
#define LinearNu_H

#include "kinematicViscosity.H"
#include "thermalIncompressibleTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kinematicViscosityModels
{

/*--------------------------------------------------------------------*\
                              Class LinearNu
\*--------------------------------------------------------------------*/

class LinearNu
:
    public kinematicViscosity
{
    // Private data

        //dimensionedScalar Cc_;

public:

    //- Runtime type information
    TypeName("linearNu");


    // Constructors

        //- construct from components
        LinearNu
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~LinearNu()
    {}


    // Member Functions
	
        //- Correct the kinematic viscosity
        virtual void correct(thermalIncompressibleTwoPhaseMixture* titpm) override;

        //- Read the transportProperties dictionary and update
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace kinematicViscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
