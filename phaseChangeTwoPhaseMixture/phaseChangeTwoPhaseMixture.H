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
    Foam::phaseChangeTwoPhaseMixture

Description

SourceFiles
    phaseChangeTwoPhaseMixture.C
    newPhaseChangeModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseChangeTwoPhaseMixture_H
#define phaseChangeTwoPhaseMixture_H

#include "thermalIncompressibleTwoPhaseMixture.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Pair.H"
#include "SaturationProperties.H"
#include "HardtWondra.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class phaseChangeTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class phaseChangeTwoPhaseMixture
:
    public thermalIncompressibleTwoPhaseMixture
{
protected:

    // Protected data

        dictionary phaseChangeTwoPhaseMixtureCoeffs_;

		//- Flag if condensation enabled
    	Switch cond_;

		//- Flag if evaporation enabled
    	Switch evap_;
	
	    //- Hardt-Wondra algorithm activation
		//  default is true
		Switch isHW_;

		//- Pointer to saturation properties
		autoPtr<SaturationProperties> satProps_;

		//- Pointer to Hardt-Wondra algorithm
		autoPtr<HardtWondra> HW_;

private:

    // Private data

    // Private Member Functions

        //- Disallow copy construct
        phaseChangeTwoPhaseMixture(const phaseChangeTwoPhaseMixture&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const phaseChangeTwoPhaseMixture&) = delete;


protected:

    // Protected data

		//- Condensation mass flux, kg/(m2 s) (+)
		volScalarField jc_;
		
		//- Evaporation mass flux, kg/(m2 s) (+)
		volScalarField je_;

		//- Implicit part of the source term for alphaEqn
		volScalarField alphaSourceSp_;

		//- Explicit part of the source term for alphaEqn
		volScalarField alphaSourceSu_;

		//- Implicit part of the source term for pEqn
		volScalarField pSourceSp_;

		//- Explicit part of the source term for pEqn
		volScalarField pSourceSu_;

		//- Implicit part of the source term for TEqn
		volScalarField TSourceSp_;

		//- Explicit part of the source term for TEqn
		volScalarField TSourceSu_;

		Switch printPhaseChange_;
		
public:

    //- Runtime type information
    TypeName("phaseChangeTwoPhaseMixture");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseChangeTwoPhaseMixture,
            components,
            (
                const volVectorField& U,
                const surfaceScalarField& phi
            ),
            (U, phi)
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<phaseChangeTwoPhaseMixture> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Constructors

        //- Construct from components
        phaseChangeTwoPhaseMixture
        (
            const word& type,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~phaseChangeTwoPhaseMixture()
    {}


    // Member Functions
	
		Switch isCond() const
		{
			return cond_;
		}

		Switch isEvap() const
		{
			return evap_;
		}

		Switch isHardtWondra() const
		{
			return isHW_;
		}

        //- Return const-access to the pressure field
		const volScalarField& p() const
		{
			return satProps_->p();
		}

        //- Return const-access to the temperature field
		const volScalarField& T() const
		{
			return satProps_->T();
		}

        //- Return const-access to the global saturation temperature
        const dimensionedScalar& TSatG() const
        {
            return satProps_->TSatG();
        }

        //- Return const-access to the saturation temperature
		const volScalarField& TSat() const
		{
			return satProps_->TSat();
		}

        //- Return const-access to the saturation pressure
		const dimensionedScalar& pSat() const
		{
			return satProps_->pSat();
		}

        //- Return const-access to the enthalpy of evaporation
		const dimensionedScalar& hEvap() const
		{
			return satProps_->hEvap();
		}

		//- Calculates phase change mass fluxes
		virtual void j() = 0;

		//- Calculates implicit part of the source term for alphaEqn
		virtual tmp<volScalarField> alphaSourceSp() = 0;

		//- Calculates explicit part of the source term for alphaEqn
		virtual tmp<volScalarField> alphaSourceSu() = 0;

		//- Calculates implicit part of the source term for pEqn
		virtual tmp<volScalarField> pSourceSp() = 0;

		//- Calculates explicit part of the source term for pEqn
		virtual tmp<volScalarField> pSourceSu() = 0;

		//- Calculates implicit part of the source term for TEqn
		virtual tmp<volScalarField> TSourceSp() = 0;

		//- Calculates explicit part of the source term for TEqn
		virtual tmp<volScalarField> TSourceSu() = 0;
		
        ////- Return the volumetric mass condensation (+) and vaporisation (-) rates 
		////  as a coefficient to multiply by (1 - alphal) for the condensation rate
        ////  and a coefficient to multiply by alphal for the vaporisation rate
		////  in alphaEqn. 
		////  Units are kg/(m3 s)
        //virtual Pair<tmp<volScalarField> > mDotAlphal() const = 0;

        ////- Return the volumetric mass condensation (+) and vaporisation (-) rates as coefficients
		////  to multiply in pEqn
		////  units are kg/(m3 s)
        //virtual Pair<tmp<volScalarField> > mDotP() const = 0;

        ////- Return the mass condensation (-) and vaporisation (+) rates as coefficients
        ////  to multiply by (T - TSat) in TEqn
		////  units are kg/(m3 s K)
        //virtual Pair<tmp<volScalarField> > mDotT() const = 0;

        //- Return the volumetric mass condensation (+) and vaporisation (-) rates 
		//  as a coefficient to multiply by (1 - alphal) for the condensation rate
        //  and a coefficient to multiply by alphal for the vaporisation rate
		//  in alphaEqn. 
		//  units are 1/s
        Pair<tmp<volScalarField>> alphaSource();

        //- Return the volumetric mass condensation (+) and vaporisation (-) rates as coefficients
		//  to multiply in pEqn
		//  units are 1/s
        Pair<tmp<volScalarField>> pSource();
		
        //- Return the mass condensation (-) and vaporisation (+) rates as coefficients
        //  to multiply by (T - TSat) in TEqn
		//  units are kg/(m s3 K)
        Pair<tmp<volScalarField>> TSource();

        //- Correct the phaseChange model
		//  In this function mass transfer model should be implemented.
		//  Algorithm used in interBoilCondFoam: 
		//  1) calculate volumetric mass fluxes in alphaEqn.H that are:
		//     i) multiplied by (1-alphal) for condensation and by alphal for
		//        evaporation and store them in mCondNoAlphal_ and mEvapNoAlphal_ 
		//     ii) store also volumetric mass fluxes multiplied by current alphal
		//         in mCondAlphal_ and mEvapAlphal_
		//     iii) store also volumetric mass fluxes that are multiplied by (T-TSat) in TEqn
		//          in mCondNoTmTSat_ and mEvapNoTmTSat_
		//     These all fluxes can be considered as from previous time-step 
		//     because they are calculate based on alphal and T from previous time-step
		//  2) use these fluxes to calculate new alphal field in alphaEqn
		//  3) if momentumPredictor is on in UEqn, calculate new U based on new alphal
		//  4) in TEqn calculate new T based on new alphal but using the stored mass fluxes 
		//  5) in pEqn calculate new p and U based on new alphal but using the stored mass fluxes
		//  6) go to next PIMPLE iteration or time-step and start from 1)
        virtual void correct();

        //- Read the transportProperties dictionary and update
        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
