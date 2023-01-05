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
    Foam::phaseChangeTwoPhaseMixtures::Xu

Description
    Xu phase change model. It is based on saturation temperature.
	It is Lee model with variable Cc and Cv.
	For T < TSat => condensation occurs
	and T > TSat => evaporation occurs

    Reference:
    \verbatim
Zhi Xu, Yaning Zhang, Bingxi Li, Jingqi Huang,
Modeling the phase change process for a two-phase closed thermosyphon by considering transient mass transfer time relaxation parameter,
International Journal of Heat and Mass Transfer,
Volume 101,
2016,
Pages 614-619,
ISSN 0017-9310,
https://doi.org/10.1016/j.ijheatmasstransfer.2016.05.075.
(http://www.sciencedirect.com/science/article/pii/S0017931015312898)
    \endverbatim

SourceFiles
    Xu.C

\*--------------------------------------------------------------------*/

#ifndef Xu_H
#define Xu_H

#include "ChoiZhang.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{

/*--------------------------------------------------------------------*\
                              Class Xu
\*--------------------------------------------------------------------*/

class Xu
:
    public ChoiZhang
{
    // Private data

        dimensionedScalar Cvvar_;


public:

    //- Runtime type information
    TypeName("Xu");


    // Constructors

        //- construct from components
        Xu
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~Xu()
    {}


    // Member Functions
	

        //- Correct the phaseChange model
        virtual void correct() override;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseChangeTwoPhaseMixtures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
