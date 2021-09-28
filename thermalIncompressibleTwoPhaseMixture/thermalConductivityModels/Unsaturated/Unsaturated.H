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
    Foam::thermalConductivityModels::Unsaturated

Description
    Unsaturated thermal conductivity model for a mixture. 
	kEff = 2*k_unsaturated

Reference:
@article{SUN20121101,
title = {Development of a vapor–liquid phase change model for volume-of-fluid method in FLUENT},
journal = {International Communications in Heat and Mass Transfer},
volume = {39},
number = {8},
pages = {1101-1106},
year = {2012},
issn = {0735-1933},
doi = {https://doi.org/10.1016/j.icheatmasstransfer.2012.07.020},
url = {https://www.sciencedirect.com/science/article/pii/S0735193312001807},
author = {Dong-Liang Sun and Jin-Liang Xu and Li Wang},
keywords = {Phase change model, Volume-of-fluid method, FLUENT},
abstract = {In this paper, a vapor–liquid phase change model is proposed for volume-of-fluid (VOF) method in FLUENT. This model is suitable for the case in which both unsaturated phase and saturated phase are present. In this model: (1) the unsaturated-phase thermal conductivity λuns and specific heat Cp,uns are actual physical parameters; (2) the saturated-phase thermal conductivity λs and specific heat Cp,s are assumed as zero and Cp,uns, respectively; (3) the interfacial mass-transfer rate m˙s=−m˙uns=2λuns(∇αuns⋅∇T)/L. Finally, the accuracy of the vapor–liquid phase change model is verified by one-dimensional Stefan problem and two-dimensional film boiling problem.}
}

SourceFiles
    Unsaturated.C

\*--------------------------------------------------------------------*/

#ifndef Unsaturated_H
#define Unsaturated_H

#include "thermalConductivity.H"
#include "thermalIncompressibleTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalConductivityModels
{

/*--------------------------------------------------------------------*\
                              Class Unsaturated
\*--------------------------------------------------------------------*/

class Unsaturated
:
    public thermalConductivity
{
    // Private data

        //dimensionedScalar Cc_;
        //dimensionedScalar Cv_;

        //dimensionedScalar mcCoeff_;
        //dimensionedScalar mvCoeff_;


public:

    //- Runtime type information
    TypeName("unsaturated");


    // Constructors

        //- construct from components
        Unsaturated
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~Unsaturated()
    {}


    // Member Functions
	
        //- Return the effective laminar thermal conductivity
        virtual tmp<volScalarField> k(const thermalIncompressibleTwoPhaseMixture* titpm) const override;

        //- Read the transportProperties dictionary and update
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalConductivityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
