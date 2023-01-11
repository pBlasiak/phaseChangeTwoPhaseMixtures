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

#include "ClausiusClapeyron.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace SaturationPropertiesModels
{
    defineTypeNameAndDebug(ClausiusClapeyron, 0);
    addToRunTimeSelectionTable(SaturationProperties, ClausiusClapeyron, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SaturationPropertiesModels::ClausiusClapeyron::ClausiusClapeyron
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    SaturationProperties(typeName, U, phi),
    R_("R", dimGasConstant, SaturationPropertiesDict_.subDict(type() + "SatPropModel"))
{ }


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::SaturationPropertiesModels::ClausiusClapeyron::calcTSat()
{
	TSat_ = 1.0/(1.0/TSatG_ - R_/hEvap_*log(max(p_/pSat_,1E-08)));
}

bool Foam::SaturationPropertiesModels::ClausiusClapeyron::read()
{
    if (SaturationProperties::read())
    {
		SaturationPropertiesDict_.subDict(type() + "SatPropModel").lookup("R") >> R_;

    	return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
