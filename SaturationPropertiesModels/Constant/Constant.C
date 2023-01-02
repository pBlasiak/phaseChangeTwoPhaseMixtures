/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "Constant.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace SaturationPropertiesModels
{
    defineTypeNameAndDebug(Constant, 0);
    addToRunTimeSelectionTable(SaturationProperties, Constant, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SaturationPropertiesModels::Constant::Constant
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    SaturationProperties(typeName, U, phi)

    //Cc_(SaturationPropertiesCoeffs_.subDict(type() + "Coeffs").lookup("Cc")),
    //Cv_(SaturationPropertiesCoeffs_.subDict(type() + "Coeffs").lookup("Cv")),
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::SaturationPropertiesModels::Constant::calcTSat()
{ }

bool Foam::SaturationPropertiesModels::Constant::read()
{
    //if (SaturationProperties::read())
    //{
        //SaturationPropertiesCoeffs_ = subDict(type() + "Coeffs");

        //SaturationPropertiesCoeffs_.lookup("Cc") >> Cc_;
        //SaturationPropertiesCoeffs_.lookup("Cv") >> Cv_;

        //mcCoeff_ = Cc_*rho2();
        //mvCoeff_ = Cv_*rho1();

        return true;
    //}
    //else
    //{
    //    return false;
    //}
}


// ************************************************************************* //
