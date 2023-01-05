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
    Foam::SaturationProperties

Description

SourceFiles
    SaturationProperties.C
    newThermalProperty.C

\*---------------------------------------------------------------------------*/

#ifndef SaturationProperties_H
#define SaturationProperties_H

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
                Class SaturationProperties Declaration
\*---------------------------------------------------------------------------*/

class SaturationProperties
: 
	IOdictionary
{

protected:

    // Protected data

        dictionary SaturationPropertiesDict_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;

        //- Reference to pressure field
        const volScalarField& p_;

        //- Reference to temperature field
        const volScalarField& T_;
		
        //- Global saturation temperature
        dimensionedScalar TSatG_;

        //- Local saturation temperature
        volScalarField TSat_;

        //- Saturation vapour pressure
        dimensionedScalar pSat_;

        //- Latent heat of vaporization
        dimensionedScalar hEvap_;


    // Private Member Functions

        //- Disallow copy construct
        SaturationProperties(const SaturationProperties&);

        //- Disallow default bitwise assignment
        void operator=(const SaturationProperties&);

public:

    //- Runtime type information
    TypeName("SaturationProperties");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            SaturationProperties,
            components,
            (
                const volVectorField& U,
                const surfaceScalarField& phi
            ),
            (U, phi)
        );


    // Selectors

        //- Return a reference to the selected thermal conductivity model
        static autoPtr<SaturationProperties> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
			const word& satProp
        );


    // Constructors

        //- Construct from components
        SaturationProperties
        (
            const word& type,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~SaturationProperties()
    {}


    // Member Functions
	
        //- Return const-access to the pressure field
        const volScalarField& p() const
        {
            return p_;
        }
		
        //- Return const-access to the temperature field
        const volScalarField& T() const
        {
            return T_;
        }

        //- Return const-access to the global saturation temperature
        const dimensionedScalar& TSatG() const
        {
            return TSatG_;
        }

        //- Return const-access to the local saturation temperature
        const volScalarField& TSat() const
        {
            return TSat_;
        }
		
        //- Return const-access to the global saturation vapour pressure
        const dimensionedScalar& pSat() const
        {
            return pSat_;
        }

        //- Return const-access to the latent heat of vaporization
        const dimensionedScalar& hEvap() const
        {
            return hEvap_;
        }

        //- Calculates saturation temperature
        virtual void calcTSat() = 0;


        //- Read the saturation properties and update
        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
