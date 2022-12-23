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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Istream;
class Ostream;
class HardtWondra;

Istream& operator>>(Istream&, HardtWondra&);
Ostream& operator<<(Ostream&, const HardtWondra&);


/*---------------------------------------------------------------------------*\
                         Class HardtWondra Declaration
\*---------------------------------------------------------------------------*/

class HardtWondra
{
    // Private Data

		dictionary HWdict_;

        //- Cutoff for mass source
		//  default is 1e-3
        scalar cutoff_;

        //- Spread for mass source
		//  default is 3
        scalar spread_;

		//- Flag for the limited alphal calculation
		Switch limitedAlphalCalculated_;

		//- Flag for the magGradLimitedAlphal calculation
		Switch magGradLimitedAlphalCalculated_;

		//- Limited liquid void fraction between 0 and 1 
		volScalarField limitedAlphal_;

		//- Gradient magnitude of limited alphal
		volScalarField magGradLimitedAlphal_;

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
        HardtWondra(const volScalarField& alpha1);


    // Selectors

        //- Select default constructed
        //static autoPtr<HardtWondra> New();


    //- Destructor
    ~HardtWondra();


    // Member Functions

    // Access

    // Check

    // Edit

    // Write


    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, HardtWondra&);
        friend Ostream& operator<<(Ostream&, const HardtWondra&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
