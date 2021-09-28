/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::thermalIncompressibleTwoPhaseMixture

Description
    A two-phase incompressible transportModel

SourceFiles
    thermalIncompressibleTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef thermalIncompressibleTwoPhaseMixture_H
#define thermalIncompressibleTwoPhaseMixture_H

//#include "incompressibleTwoPhaseMixture.H"
#include "immiscibleIncompressibleSmoothedTwoPhaseMixture.H"
#include "thermalConductivity.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class thermalIncompressibleTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class thermalIncompressibleTwoPhaseMixture
:
    public immiscibleIncompressibleSmoothedTwoPhaseMixture
{
protected:

    // Protected data


        //- Thermal variables
        dimensionedScalar k1_;
        dimensionedScalar k2_;

        dimensionedScalar cp1_;
        dimensionedScalar cp2_;

		autoPtr<thermalConductivity> thermalConductivityModel_;


public:

    TypeName("thermalIncompressibleTwoPhaseMixture");


    // Constructors

        //- Construct from U and phi
        thermalIncompressibleTwoPhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~thermalIncompressibleTwoPhaseMixture()
    {}


    // Access function

        //- Return const-access to phase1 kappa
        const dimensionedScalar& k1() const
        {
            return k1_;
        };

        //- Return const-access to phase2 kappa
        const dimensionedScalar& k2() const
        {
            return k2_;
        };

        //- Return const-access to phase1 Cp
        const dimensionedScalar& cp1() const
        {
            return cp1_;
        };

        //- Return const-access to phase1 Cp
        const dimensionedScalar& cp2() const
        {
            return cp2_;
        };

        //- Return the effective laminar thermal conductivity
		tmp<volScalarField> k() const
		{
             return thermalConductivityModel_->k(this);
		}

        //- Return the effective face-interpolated laminar thermal conductivity
		tmp<surfaceScalarField> kf() const
		{
	        return tmp<surfaceScalarField>
	        (
	        	new surfaceScalarField
	        	(
	        		"kf",
	        		(
	        		    fvc::interpolate(k())
	        		)
	        	)
	        );
		}

		//- Return density of the mixture
		tmp<volScalarField> rho() const;

		//- Return specific heat of the mixture
		tmp<volScalarField> cp() const;

        //- Correct the transport and interface properties
        virtual void correct()
        {
            immiscibleIncompressibleSmoothedTwoPhaseMixture::correct();
        }

        //- Read base transportProperties dictionary
        virtual bool read();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
