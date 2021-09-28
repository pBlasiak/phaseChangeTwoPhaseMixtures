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
    Foam::phaseChangeTwoPhaseMixtures::Lee

Description
    Lee phase change model. It is based on saturation temperature.
	For T < TSat => condensation occurs
	and T > TSat => evaporation occurs

    Reference:
    \verbatim
        Lee, W.H.
        "A Pressure Iteration Scheme for Two-Phase Modeling. Technical Report LA-UR 79-975",
        Los Alamos Scientific Laboratory, Los Alamos, New Mexico, 1979.
    \endverbatim

SourceFiles
    Lee.C

\*--------------------------------------------------------------------*/

#ifndef Lee_H
#define Lee_H

#include "phaseChangeTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{

/*--------------------------------------------------------------------*\
                              Class Lee
\*--------------------------------------------------------------------*/

class Lee
:
    public phaseChangeTwoPhaseMixture
{
    // Protected data
protected:

        dimensionedScalar Cc_;
        dimensionedScalar Cv_;

        dimensionedScalar mcCoeff_;
        dimensionedScalar mvCoeff_;


public:

    //- Runtime type information
    TypeName("Lee");


    // Constructors

        //- construct from components
        Lee
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~Lee()
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
