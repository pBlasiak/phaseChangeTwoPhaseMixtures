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

#include "incompressibleTwoPhaseMixture.H"
//#include "immiscibleIncompressibleSmoothedTwoPhaseMixture.H"
#include "thermalProperty.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class thermalIncompressibleTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class thermalIncompressibleTwoPhaseMixture
:
    public incompressibleTwoPhaseMixture
{
protected:

    // Protected data


        //- Thermal variables
        dimensionedScalar k1_;
        dimensionedScalar k2_;

        dimensionedScalar cp1_;
        dimensionedScalar cp2_;

		autoPtr<thermalProperty> nuModel_;
		autoPtr<thermalProperty> rhoModel_;
		autoPtr<thermalProperty> kModel_;
		autoPtr<thermalProperty> cpModel_;
		//autoPtr<thermalProperty> rhoCpModel_; //????

		volScalarField& rho_;
		volScalarField& k_;
		volScalarField& cp_;

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

		//- Calculates kinematic viscosity of the mixture
		void calcNu() 
		{
			tmp<volScalarField> tnu1 = nuModel1_->nu();
			const volScalarField& nu1 = tnu1();
			tmp<volScalarField> tnu2 = nuModel2_->nu();
			const volScalarField& nu2 = tnu2();

			nu_ = nuModel_->calcThermProp
			(
			   this,
			   nu1,
			   nu2
			);
		} 

		//- Calculates density of the mixture
		void calcRho() 
		{
            rho_ = rhoModel_->calcThermProp
			(
			   this,
			   &thermalIncompressibleTwoPhaseMixture::rho1,
			   &thermalIncompressibleTwoPhaseMixture::rho2
			);
		} 

        //- Calculates the effective laminar thermal conductivity
		void calcK() 
		{
             k_ = kModel_->calcThermProp
			 (
				this,
				&thermalIncompressibleTwoPhaseMixture::k1,
				&thermalIncompressibleTwoPhaseMixture::k2
			 );
		}

        ////- Return the effective face-interpolated laminar thermal conductivity
		//tmp<surfaceScalarField> kf() const
		//{
	    //    return tmp<surfaceScalarField>
	    //    (
	    //    	new surfaceScalarField
	    //    	(
	    //    		"kf",
	    //    		(
	    //    		    fvc::interpolate(k())
	    //    		)
	    //    	)
	    //    );
		//}

		//- Calculates specific heat of the mixture
		void calcCp() 
		{
             cp_ = cpModel_->calcThermProp
			 (
				this,
				&thermalIncompressibleTwoPhaseMixture::cp1,
				&thermalIncompressibleTwoPhaseMixture::cp2
			 );
		}

		//- Returns const-access to mixture density
		tmp<volScalarField> rho() const
		{
			 return rho_;
		}

		//- Returns const-access to mixture thermal conductivity
		tmp<volScalarField> k() const
		{
			 return k_;
		}

		//- Returns const-access to mixture specific heat capacity at constant pressure
		tmp<volScalarField> cp() const
		{
			 return cp_;
		}

        //- Correct the transport properties
        virtual void correct() override;

        //- Read base transportProperties dictionary
        virtual bool read();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
