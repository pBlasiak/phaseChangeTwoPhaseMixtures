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
    Foam::thermalPropertyModels::HarmonicDensityWeighted

Description
    HarmonicDensityWeighted thermal conductivity model for a mixture. 

    Reference:
    \verbatim
	Cosden, Ian A., "A Hybrid Atomistic-Continuum Model for Liquid-Vapor Phase Change" (2013). 
	Publicly Accessible Penn Dissertations. 846.  
    https://repository.upenn.edu/edissertations/846
    \endverbatim

SourceFiles
    HarmonicDensityWeighted.C

\*--------------------------------------------------------------------*/

#ifndef HarmonicDensityWeighted_H
#define HarmonicDensityWeighted_H

#include "thermalProperty.H"
#include "thermalIncompressibleTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPropertyModels
{

/*--------------------------------------------------------------------*\
                              Class HarmonicDensityWeighted
\*--------------------------------------------------------------------*/

class HarmonicDensityWeighted
:
    public thermalProperty
{
    // Private data

        //dimensionedScalar Cc_;
        //dimensionedScalar Cv_;

        //dimensionedScalar mcCoeff_;
        //dimensionedScalar mvCoeff_;


public:

    //- Runtime type information
    TypeName("harmonicDensityWeighted");


    // Constructors

        //- construct from components
        HarmonicDensityWeighted
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~HarmonicDensityWeighted()
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
