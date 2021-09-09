/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "singlePhaseTransportModelYeldStress.H"
#include "viscosityModelYeldStress.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(singlePhaseTransportModelYeldStress, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singlePhaseTransportModelYeldStress::singlePhaseTransportModelYeldStress
(
    const volVectorField& U,
    //const volScalarField& T,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    viscosityModelYeldStressPtr_(viscosityModelYeldStress::New("nu", "tauY",  *this, U, phi))//mexi aqui
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::singlePhaseTransportModelYeldStress::~singlePhaseTransportModelYeldStress()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::singlePhaseTransportModelYeldStress::nu() const
{
    return viscosityModelYeldStressPtr_->nu();
}

Foam::tmp<Foam::volScalarField>
Foam::singlePhaseTransportModelYeldStress::tauY() const
{
    return viscosityModelYeldStressPtr_->tauY();
}


Foam::tmp<Foam::scalarField>
Foam::singlePhaseTransportModelYeldStress::nu(const label patchi) const
{
    return viscosityModelYeldStressPtr_->nu(patchi);
}

Foam::tmp<Foam::scalarField>
Foam::singlePhaseTransportModelYeldStress::tauY(const label patchi) const
{
    return viscosityModelYeldStressPtr_->tauY(patchi);
}


void Foam::singlePhaseTransportModelYeldStress::correct()
{
    viscosityModelYeldStressPtr_->correct();
}


bool Foam::singlePhaseTransportModelYeldStress::read()
{
    if (regIOobject::read())
    {
        return viscosityModelYeldStressPtr_->read(*this);
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
