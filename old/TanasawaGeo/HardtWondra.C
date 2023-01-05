/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
    Copyright (C) 2020 Henning Scheufler
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

#include "HardtWondra.H"
#include "addToRunTimeSelectionTable.H"
//#include "twoPhaseMixtureEThermo.H"
#include "cutCellIso.H"
#include "volPointInterpolation.H"
#include "calculatedFvPatchFields.H"
#include "wallPolyPatch.H"
#include "fvcSmooth.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
//namespace phaseChangeTwoPhaseMixtures
namespace phaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(HardtWondra, 0);
    addToRunTimeSelectionTable
    (
        //temperaturePhaseChangeTwoPhaseMixture,
        phaseChangeTwoPhaseMixture,
        HardtWondra,
        components
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixtures::HardtWondra::
HardtWondra
(
            const volVectorField& U,
            const surfaceScalarField& phi
    //const thermoIncompressibleTwoPhaseMixture& mixture,
    //const fvMesh& mesh
)
:
    //temperaturePhaseChangeTwoPhaseMixture(mixture, mesh),
    phaseChangeTwoPhaseMixture(typeName, U, phi),
    cond_("condensation", phaseChangeTwoPhaseMixtureCoeffs_.subDict(type() + "Coeffs")),
    evap_("evaporation", phaseChangeTwoPhaseMixtureCoeffs_.subDict(type() + "Coeffs")),
    gamma_
    (
        "gamma",
        //dimless, optionalSubDict(type() + "Coeffs")
        dimless, phaseChangeTwoPhaseMixtureCoeffs_.subDict(type() + "Coeffs")
	),
    HTC_
    (
        "heatTransferCoeff",
		(2.0*gamma_/(2.0 - gamma_)/sqrt(2.0*M_PI*R_)*hEvap_*hEvap_*rho2())
		//heatResistance_( (2-gamma_)/(2*gamma_)*sqrt(2.0*M_PI*R_)/hEvap_/hEvap_/rho2() )
        //dimPower/dimArea/dimTemperature, optionalSubDict(type() + "Coeffs")
    ),

    interfaceArea_
    (
        IOobject
        (
            "interfaceArea",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dimless/dimLength, Zero)
    ),

    mDotc_
    (
        IOobject
        (
            "mDotc",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dimDensity/dimTime, Zero)
    ),

    mDote_
    (
        IOobject
        (
            "mDote",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dimDensity/dimTime, Zero)
    ),

    mDotcSpread_
    (
        IOobject
        (
            "mDotcSpread",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dimDensity/dimTime, Zero)
    ),

    mDoteSpread_
    (
        IOobject
        (
            "mDoteSpread",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dimDensity/dimTime, Zero)
    ),

    spread_
    (
        phaseChangeTwoPhaseMixtureCoeffs_.subDict(type() + "Coeffs").get<scalar>("spread")
    )
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//Foam::Pair<Foam::tmp<Foam::volScalarField>>
//Foam::phaseChangeTwoPhaseMixtures::HardtWondra::
//vDotAlphal() const
//{
//    dimensionedScalar alphalCoeff(1.0/rho1());
//
//    return Pair<tmp<volScalarField>>
//    (
//        (alphalCoeff*mDotc_)/(alpha2() + SMALL),
//       -(alphalCoeff*mDote_)/(alpha1() + SMALL)
//    );
//}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixtures::HardtWondra::
mDotAlphal() const
{
    volScalarField limitedAlpha1
    (
        min(max(alpha1(), scalar(0)), scalar(1))
    );

    volScalarField limitedAlpha2
    (
        min(max(alpha2(), scalar(0)), scalar(1))
    );

    //return Pair<tmp<volScalarField>>
    //(
    //  //  (mDotc_/(limitedAlpha2 + SMALL)),
    //  // -(mDote_/(limitedAlpha1 + SMALL))
    //    (mDotc_*scalar(1)),
    //   -(mDote_*scalar(1))
    //);
	return Pair<tmp<volScalarField>>
	(
		tmp<volScalarField>(mCondNoAlphal_),
		tmp<volScalarField>(mEvapNoAlphal_)
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixtures::HardtWondra::
//mDot() const
mDotP() const
{
    return Pair<tmp<volScalarField>>
    (
        //tmp<volScalarField>(mDotc_),
        //tmp<volScalarField>(mDote_)
        //tmp<volScalarField>(mDotcSpread_),
        //tmp<volScalarField>(-mDoteSpread_)
        mDotcSpread_*pos(p_-pSat_)/max(p_-pSat_,1E-6*pSat_),
	    -mDoteSpread_*neg(p_-pSat_)/max(pSat_-p_,1E-6*pSat_)
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixtures::HardtWondra::
//mDotDeltaT() const
mDotT() const
{
   //const twoPhaseMixtureEThermo& thermo =
   //     refCast<const twoPhaseMixtureEThermo>
   //     (
   //         U.mesh().lookupObject<basicThermo>(basicThermo::dictName)
   //     );

    //const volScalarField& T = U.mesh().lookupObject<volScalarField>("T");

    //const dimensionedScalar& TSat = thermo.TSat();

    //Pair<tmp<volScalarField>> mDotce(mDot());

    //return Pair<tmp<volScalarField>>
    //(
    //    mDotc_*pos(TSat_ - T_.oldTime())/(TSat_ - T_.oldTime()),
    //   -mDote_*pos(T_.oldTime() - TSat_)/(T_.oldTime() - TSat_)
    //);
	return Pair<tmp<volScalarField> >
	(
	    tmp<volScalarField>(mCondNoTmTSat_),
	    tmp<volScalarField>(mEvapNoTmTSat_)
	);
}


//Foam::tmp<Foam::fvScalarMatrix>
//Foam::phaseChangeTwoPhaseMixtures::HardtWondra::
//TSource() const
//{
//    //const volScalarField& T = U.mesh().lookupObject<volScalarField>("T");
//
//    auto tTSource = tmp<fvScalarMatrix>::New(T_, dimEnergy/dimTime);
//    auto& TSource = tTSource.ref();
//
//    //const twoPhaseMixtureEThermo& thermo =
//    //    refCast<const twoPhaseMixtureEThermo>
//    //    (
//    //        U.mesh().lookupObject<basicThermo>(basicThermo::dictName)
//    //    );
//
//    //const dimensionedScalar& TSat = thermo.TSat();
//
//    // interface heat resistance
//    volScalarField IHRcoeff(interfaceArea_*R_);
//
//    TSource = fvm::Sp(IHRcoeff, T_) - IHRcoeff*TSat_;
//
//    return tTSource;
//}


void Foam::phaseChangeTwoPhaseMixtures::HardtWondra::
correct()
{
	phaseChangeTwoPhaseMixture::correct();

    // Update Interface
    updateInterface();

    // Update mDotc_ and mDote_
    //const volScalarField& T = U.mesh().lookupObject<volScalarField>("T");

    //const twoPhaseMixtureEThermo& thermo =
    //    refCast<const twoPhaseMixtureEThermo>
    //    (
    //        U.mesh().lookupObject<basicThermo>(basicThermo::dictName)
    //    );

    //const dimensionedScalar& TSat = thermo.TSat();
    const dimensionedScalar T0(dimTemperature, Zero);
    volScalarField limitedAlpha1 = min(max(alpha1(), scalar(0)), scalar(1));

    //dimensionedScalar L = mixture_.Hf2() - mixture_.Hf1();

    // interface mass fluxes
	// q = HTC*(T-TSat)
	// M = q/hEvap
	// m = M*interfaceArea
	if (cond_)
	{
		//mDotc_ = interfaceArea_*HTC_*max(TSat_ - T_, T0)/hEvap_/sqrt(pow(TSat_,3.0));
		mCondNoAlphal_ = interfaceArea_*HTC_*max(TSat_ - T_, T0)/hEvap_/sqrt(pow(TSat_,3.0));
		mCondAlphal_   = mCondNoAlphal_*(1-limitedAlpha1);
		mDotc_ = mCondAlphal_;
		mCondNoTmTSat_ = -neg(T_ - TSat_)*interfaceArea_*HTC_/hEvap_/sqrt(pow(TSat_,3.0));
	}
	if (evap_)
	{
		//mDote_ = interfaceArea_*HTC_*max(T_ - TSat_, T0)/hEvap_/sqrt(pow(TSat_,3.0));
		mEvapNoAlphal_ = -interfaceArea_*HTC_*max(T_ - TSat_, T0)/hEvap_/sqrt(pow(TSat_,3.0));
		mEvapAlphal_   = mEvapNoAlphal_*limitedAlpha1;
		mDote_ = mEvapAlphal_;
		mEvapNoTmTSat_ = pos(T_ - TSat_)*interfaceArea_*HTC_/hEvap_/sqrt(pow(TSat_,3.0));
	}

    forAll(mDotc_, celli)
    {
        scalar rhobyDt = rho1().value()/T_.mesh().time().deltaTValue();
        scalar maxEvap = alpha1()[celli]*rhobyDt; // positive
        scalar maxCond = -alpha2()[celli]*rhobyDt; // negative
        mDote_[celli] = min(max(mDote_[celli], maxCond), maxEvap);
        mDotc_[celli] = min(max(mDotc_[celli], maxCond), maxEvap);
    }

    // Calculate the spread sources

    dimensionedScalar D
    (
        "D",
        dimArea,
        spread_/sqr(gAverage(T_.mesh().nonOrthDeltaCoeffs()))
    );


    //const volScalarField& alpha1 = alpha1();
    //const volScalarField& alpha2 = alpha2();

    const dimensionedScalar MDotMin("MdotMin", mDotc_.dimensions(), 1e-3);

    if (max(mDotc_) > MDotMin)
    {
        fvc::spreadSource
        (
            mDotcSpread_,
            mDotc_,
            alpha1(),
            alpha2(),
            D,
            1e-3
        );
    }

    if (max(mDote_) > MDotMin)
    {
        fvc::spreadSource
        (
            mDoteSpread_,
            mDote_,
            alpha1(),
            alpha2(),
            D,
            1e-3
        );
    }
}


void Foam::phaseChangeTwoPhaseMixtures::HardtWondra::
updateInterface()
{
    //const volScalarField& T = U.mesh().lookupObject<volScalarField>("T");
    //const twoPhaseMixtureEThermo& thermo =
    //refCast<const twoPhaseMixtureEThermo>
    //(
    //    U.mesh().lookupObject<basicThermo>(basicThermo::dictName)
    //);

    //const dimensionedScalar& TSat = thermo.TSat();

    // interface heat resistance
    // Interpolating alpha1 cell centre values to mesh points (vertices)
    scalarField ap
    (
        volPointInterpolation::New(T_.mesh()).interpolate(alpha1()) // ??
    );

    cutCellIso cutCell(T_.mesh(), ap);                              // ??

    forAll(interfaceArea_, celli)
    {
        label status = cutCell.calcSubCell(celli, 0.5);
        interfaceArea_[celli] = 0;
        if (status == 0) // cell is cut
        {
            interfaceArea_[celli] =
                mag(cutCell.faceArea())/T_.mesh().V()[celli];     // ??
        }
    }

    const polyBoundaryMesh& pbm = T_.mesh().boundaryMesh();     // ??

    forAll(pbm, patchi)
    {
        if (isA<wallPolyPatch>(pbm[patchi]))
        {
            const polyPatch& pp = pbm[patchi];
            forAll(pp.faceCells(),i)
            {
                const label pCelli = pp.faceCells()[i];

                if
                (
                    (TSat_[pCelli] - T_[pCelli]) > 0
                  && alpha1()[pCelli] < 0.9
                )
                {
                    interfaceArea_[pCelli] =
                        mag(pp.faceAreas()[i])/T_.mesh().V()[pCelli];     // ??
                }
            }
        }
    }
}


//Foam::Pair<Foam::tmp<Foam::volScalarField>>
//Foam::phaseChangeTwoPhaseMixtures::HardtWondra::
//vDot() const
//{
//
//    dimensionedScalar pCoeff(1.0/rho1() - 1.0/rho2());
//
//    return Pair<tmp<volScalarField>>
//    (
//        pCoeff*mDotcSpread_,
//       -pCoeff*mDoteSpread_
//    );
//}


bool Foam::phaseChangeTwoPhaseMixtures::HardtWondra::
read()
{
    //if (temperaturePhaseChangeTwoPhaseMixture::read())
    if (phaseChangeTwoPhaseMixture::read())
    {
		// nie wiem czy tu nie czyta teraz z transportProperties
        optionalSubDict(type() + "Coeffs").readEntry("gamma", gamma_);
        optionalSubDict(type() + "Coeffs").readEntry("HTC", HTC_);
        optionalSubDict(type() + "Coeffs").readEntry("spread", spread_);
        optionalSubDict(type() + "Coeffs").readEntry("condensation", cond_);
        optionalSubDict(type() + "Coeffs").readEntry("evaporation", evap_);

        return true;
    }

    return false;
}


// ************************************************************************* //
