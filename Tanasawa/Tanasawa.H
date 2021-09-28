/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU Generac License as published by the
    Free Software Foundation; either  2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the ho it will be useful, but WITHOUT
    ANY WARRANTY; without even the imarranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.he GNU General Public License
    for more details.

    You should have received a copy oNU General Public License
    along with OpenFOAM; if not, write Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floorn, MA 02110-1301 USA

Class
    Foam::phaseChangeTwoPhaseMixtures::Tanasawa

Description
    Tanasawa phase change model.
	The model supports condensation and evaporation mass transfer.

    Reference:
    \verbatim
        Tanasawa, I.
        "Advances in condensation heat transfer",
        Advances in Heat Transfer, vol. 21, 1991, pp. 55-139.
    \endverbatim


SourceFiles
    Tanasawa.C

\*--------------------------------------------------------------------*/

#ifndef Tanasawa_H
#define Tanasawa_H

#include "phaseChangeTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{

/*--------------------------------------------------------------------*\
                              Class Tanasawa
\*--------------------------------------------------------------------*/

class Tanasawa
:
    public phaseChangeTwoPhaseMixture
{
    // Private data

    	Switch cond_;

    	Switch evap_;

    	dimensionedScalar gamma_;

        dimensionedScalar mCoeff_;
		
    // Private member functions

        volScalarField calcGradAlphal() const;

public:

    //- Runtime type information
    TypeName("Tanasawa");


    // Constructors

        //- construct from components
        Tanasawa
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~Tanasawa()
        {}


    // Member Functions


        //- Return the volumetric mass condensation (+) and vaporisation (-) rates 
		//  as a coefficient to multiply by (1 - alphal) for the condensation rate
        //  and a coefficient to multiply by alphal for the vaporisation rate
		//  in alphaEqn. 
		//  In this function the mass transfer model should be implemented
		//  without alphal and (1 - alphal) terms included and stored in 
		//  mCondNoAlphal_ and mEvapNoAlphal_ variables.
		//  units are kg/(m3 s)
		//
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
        virtual Pair<tmp<volScalarField> > mDotAlphal();

        //- Return the volumetric mass condensation (+) and vaporisation (-) rates as coefficients
		//  to multiply in pEqn
		//  units are kg/(m3 s)
        virtual Pair<tmp<volScalarField> > mDotP() const;

        //- Return the mass condensation (-) and vaporisation (+) rates as coefficients
        //  to multiply by (T - TSat) in TEqn
		//  units are kg/(m3 s K)
        virtual Pair<tmp<volScalarField> > mDotT() const;

        //- Read the transportProperties dictionary and update
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseChangeTwoPhaseMixtures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
