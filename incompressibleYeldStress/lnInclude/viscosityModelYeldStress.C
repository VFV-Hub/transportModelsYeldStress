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

#include "viscosityModelYeldStress.H"
#include "volFields.H"
#include "fvcGrad.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(viscosityModelYeldStress, 0);
    defineRunTimeSelectionTable(viscosityModelYeldStress, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModelYeldStress::viscosityModelYeldStress
(
    const word& name,
    const word& nameDois,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    //const volScalarField& T,
    const surfaceScalarField& phi
)
:
    name_(name),
    nameDois_(nameDois),
    viscosityProperties_(viscosityProperties),
    U_(U),
    //T_(T),
    phi_(phi)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::viscosityModelYeldStress::strainRate() const
{
    return sqrt(2.0)*mag(symm(fvc::grad(U_)));
}
/*Foam::tmp<Foam::volScalarField> Foam::viscosityModelYeldStress::dTdt() const
{
    return fvc::ddt(T_);
}*/


bool Foam::viscosityModelYeldStress::read(const dictionary& viscosityProperties)
{
    viscosityProperties_ = viscosityProperties;

    return true;
}


// ************************************************************************* //
