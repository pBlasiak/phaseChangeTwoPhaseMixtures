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
    Foam::SaturationPropertiesModels::ClausiusClapeyron

Description
    ClausiusClapeyron thermal property model for a mixture. 
	thermProp = alpha1*thermProp1 + alpha2*thermProp2

SourceFiles
    ClausiusClapeyron.C

\*--------------------------------------------------------------------*/

#ifndef ClausiusClapeyron_H
#define ClausiusClapeyron_H

#include "../SaturationProperties/SaturationProperties.H"
#include "thermalIncompressibleTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace SaturationPropertiesModels
{

/*--------------------------------------------------------------------*\
                              Class ClausiusClapeyron
\*--------------------------------------------------------------------*/

class ClausiusClapeyron
:
    public SaturationProperties
{
    // Private data

        //- Specific gas constant
        dimensionedScalar R_;
		

public:

    //- Runtime type information
    TypeName("ClausiusClapeyron");


    // Constructors

        //- construct from components
        ClausiusClapeyron
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~ClausiusClapeyron()
    {}


    // Member Functions
	
        //- Return const-access to the specific gas constant
        const dimensionedScalar& R() const
        {
            return R_;
        }

        //- Return the effective laminar thermal property
        virtual void calcTSat () override; 

        //- Read the transportProperties dictionary and update
        virtual bool read() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SaturationPropertiesModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
