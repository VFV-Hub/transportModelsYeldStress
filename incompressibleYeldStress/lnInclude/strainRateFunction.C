/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "strainRateFunction.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModelYeldStresss
{
    defineTypeNameAndDebug(strainRateFunction, 0);

    addToRunTimeSelectionTable
    (
        viscosityModelYeldStress,
        strainRateFunction,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModelYeldStresss::strainRateFunction::strainRateFunction
(
    const word& name,
    const word& nameDois,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    //const volScalarField& T,
    const surfaceScalarField& phi
)
:
    viscosityModelYeldStress(name, nameDois, viscosityProperties, U, phi),
    strainRateFunctionCoeffs_
    (
        viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),
    strainRateFunction_
    (
        Function1<scalar>::New("function", strainRateFunctionCoeffs_)
    ),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar(name, dimViscosity, 0)
    ),
    tauY_
    (
        IOobject
        (
            nameDois,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar(nameDois, dimViscosity/dimTime, 0)
    )
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModelYeldStresss::strainRateFunction::nu() const
{
    return nu_;
}

Foam::tmp<Foam::volScalarField>
Foam::viscosityModelYeldStresss::strainRateFunction::tauY() const
{
    return tauY_;
}

Foam::tmp<Foam::scalarField>
Foam::viscosityModelYeldStresss::strainRateFunction::nu(const label patchi) const
{
    return nu_.boundaryField()[patchi];
}

Foam::tmp<Foam::scalarField>
Foam::viscosityModelYeldStresss::strainRateFunction::tauY(const label patchi) const
{
    return tauY_.boundaryField()[patchi];
}


void Foam::viscosityModelYeldStresss::strainRateFunction::correct()
{
    tmp<volScalarField> tsigma = strainRate();
    const volScalarField& sigma = tsigma();

    nu_.primitiveFieldRef() = strainRateFunction_->value(sigma());

    volScalarField::Boundary& nuBf = nu_.boundaryFieldRef();
    const volScalarField::Boundary& sigmaBf = sigma.boundaryField();

    forAll(nuBf, patchi)
    {
        nuBf[patchi] = strainRateFunction_->value(sigmaBf[patchi]);
    }
}


bool Foam::viscosityModelYeldStresss::strainRateFunction::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModelYeldStress::read(viscosityProperties);

    strainRateFunctionCoeffs_ = viscosityProperties.optionalSubDict
    (
        typeName + "Coeffs"
    );

    strainRateFunction_.clear();
    strainRateFunction_ = Function1<scalar>::New
    (
        "function",
        strainRateFunctionCoeffs_
    );

    return true;
}


// ************************************************************************* //
