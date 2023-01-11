/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 AUTHOR,AFFILIATION
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
    Foam::HardtWondra

Description
	Algorithm for spreading out source terms for the phase change
	simulations.

References
@article{HARDT20085871,
title = {Evaporation model for interfacial flows based on a continuum-field representation of the source terms},
journal = {Journal of Computational Physics},
volume = {227},
number = {11},
pages = {5871-5895},
year = {2008},
issn = {0021-9991},
doi = {https://doi.org/10.1016/j.jcp.2008.02.020},
url = {https://www.sciencedirect.com/science/article/pii/S0021999108001228},
author = {S. Hardt and F. Wondra},
keywords = {Evaporation model, Interfacial flow, Volume-of-fluid method, Stefan problem, Droplet evaporation, Film boiling},
}

SourceFiles
    HardtWondra.C

\*---------------------------------------------------------------------------*/

#ifndef HardtWondra_H
#define HardtWondra_H

#include "fvCFD.H"
#include "typeInfo.H"
#include "SaturationProperties.H"
//#include "SaturationPropertiesModels/SaturationProperties/SaturationProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
//class Istream;
//class Ostream;
class phaseChangeTwoPhaseMixture;
//
//Istream& operator>>(Istream&, HardtWondra&);
//Ostream& operator<<(Ostream&, const HardtWondra&);


/*---------------------------------------------------------------------------*\
                         Class HardtWondra Declaration
\*---------------------------------------------------------------------------*/

class HardtWondra
{
    // Private Data

		dictionary HWdict_;

        //- Reference to liquid void fraction field
        const volScalarField& alphalRef_;

        //- Cutoff for mass source
		//  default is 1e-3
        scalar cutoff_;

        //- Spread for mass source
		//  default is 3
        scalar spread_;

		//TODO: czy to jest potrzebne? Czy nie wprowadza bledow?
		//      czy nie powinno byc uzywane "raw alpha1()"?
		//- Limited liquid void fraction between 0 and 1 
		volScalarField limitedAlphal_;

		//- Gradient magnitude of limited alphal
		volScalarField magGradLimitedAlphal_;

		//- Mass source term for condensation
		volScalarField rhoSourcel_;
		
		//- Entalpy source term for condensation
		volScalarField hSourcel_;

		//- Access to saturation data
		const SaturationProperties& mixtureSatProps_;

		//- Access to mixture data
		const phaseChangeTwoPhaseMixture& mixture_;

    // Private Member Functions

        //- No copy construct
        HardtWondra(const HardtWondra&) = delete;

        //- No copy assignment
        void operator=(const HardtWondra&) = delete;

		//- Calculates limited alphal
		void calcLimitedAlphal();

		//- Calculates magnitude of gradient of limited alphal
		void calcMagGradLimitedAlphal();


public:

    // Static Data Members

        //- Static data staticData
        //static const dataType staticData;


    // Generated Methods

//        //- Default construct
//        HardtWondra() = default;
//
//        //- Copy construct
//        HardtWondra(const HardtWondra&) = default;
//
//        //- Copy assignment
//        HardtWondra& operator=(const HardtWondra&) = default;


    // Constructors

        //- Default construct
        //HardtWondra();

        //- Construct from components
        HardtWondra
		(
			const volScalarField& alpha1,
			const SaturationProperties& sat,
			const phaseChangeTwoPhaseMixture& mix
		);


    // Selectors

        //- Select default constructed
        //autoPtr<HardtWondra> New();
        //static autoPtr<HardtWondra> New();


    //- Destructor
    ~HardtWondra();


    // Member Functions
		
		const volScalarField& magGradLimitedAlphal() const
		{
			return magGradLimitedAlphal_;
		}

		const volScalarField& rhoSourcel() const
		{
			return rhoSourcel_;
		}

		const volScalarField& hSourcel() const
		{
			return hSourcel_;
		}

	    //- Hardt-Wondra algorithm
		//  OLD implementation WAS:
		//  - no spread mass source terms in alphaEqn
		//  - spread mass source terms in pEqn
		//  - no spread mass source terms in TEqn
		//  Additionally in TEqn mass source term the correction
		//  term h_c is ommited (Eq. (23)) because it gives 
		//  worse results
		//
		//  Ustalenia:
		//  1) przy MULESCorr = false uzywany jest explicitSolve
		//     bez correct
		//  2) wydaje się ze dla 1) czlon zrodlowy nie jest przemnazany 
		//     przez alphal (nie ma nigdzie fvm::Sp)
		//  3) evapVOFHardt wstawia do alphaEqn surowy czlon zrodlowy (bez spread)
		//     Smeared source term jest wstawiany tylko do pEqn i to jawnie bez
		//     fvm::Sp
		//
		//  Wnioski:
		//  1) na poczatek zrob algorytm z uzyciem tylko explicitSolve 
		//     (MULESCorr = false)
		//  2) zrob implementacje jak w evapVOFHardt (spread tylko dla pEqn i explicit)
		//  3) zaimplementuj caly czlon h
		//  4) niech funkcja spread przyjmuje jc i je
		//
		void spread
		(
			const volScalarField& condMassFlux, 
			const volScalarField& evapMassFlux
		);

		void correct();

		void read();

    // Access

    // Check

    // Edit

    // Write


    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators

        //friend Istream& operator>>(Istream&, HardtWondra&);
        //friend Ostream& operator<<(Ostream&, const HardtWondra&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
