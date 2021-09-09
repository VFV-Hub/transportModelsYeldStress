/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "fractalDimension.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModelYeldStresss
{
    defineTypeNameAndDebug(fractalDimension, 0);

    addToRunTimeSelectionTable
    (
        viscosityModelYeldStress,
        fractalDimension,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModelYeldStresss::fractalDimension::calcNu() const
{

	dimensionedScalar TVSMALL("VSMALL", dimTemperature, VSMALL);	
	dimensionedScalar rtVSMALL("VSMALL", dimless/dimTime, VSMALL);
	
    tmp<volScalarField> sr(strainRate());
    //tmp<volScalarField> q(dTdt());
	
	const volScalarField& T=U_.mesh().lookupObject<volScalarField>("T");
	const volScalarField& f=U_.mesh().lookupObject<volScalarField>("f");
	
	return
	(
		min
        	(
			nuMax_,
			(
			nu0_*exp(En_/max(T, TVSMALL))*sr()
			//+CTau_*pow((f0_-f),2.0/(3.0-D0_*exp(D1_*q())))
			)
        	/(
				max(
					sr(), rtVSMALL
					)
		 	)
			)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::viscosityModelYeldStresss::fractalDimension::calcTauY() const
{	
	dimensionedScalar rtVSMALL("VSMALL", dimless/dimTime, VSMALL);
	
    tmp<volScalarField> sr(strainRate());
    //tmp<volScalarField> q(dTdt());
	
	const volScalarField& f=U_.mesh().lookupObject<volScalarField>("f");
	
	return
	(
		min
        	(
			nuMax_,
			(
			0
			//CTau_*pow((f0_-f),2.0/(3.0-D0_*exp(D1_*q())))
			)
        	/(
				max(
					sr(), rtVSMALL
					)
		 	)
			)
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModelYeldStresss::fractalDimension::fractalDimension
(
    const word& name,
    const word& nameDois,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const volScalarField& T,
    const surfaceScalarField& phi
)
:
    viscosityModelYeldStress(name, nameDois, viscosityProperties, U, T, phi),
    fractalDimensionCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
//inicio da edição    
	//CTau_("CTau", dimViscosity/dimTime, fractalDimensionCoeffs_),
	nu0_("nu0", dimViscosity, fractalDimensionCoeffs_),	
	En_("En", dimTemperature, fractalDimensionCoeffs_),
	nuMax_("nuMax", dimViscosity, fractalDimensionCoeffs_),
	CTau_("CTau", dimViscosity/dimTime, fractalDimensionCoeffs_),
	f0_("f0", dimless, fractalDimensionCoeffs_),
	D0_("D0", dimless, fractalDimensionCoeffs_),
	D1_("D1", dimTime/dimTemperature, fractalDimensionCoeffs_),
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
        calcNu()
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
        calcTauY()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModelYeldStresss::fractalDimension::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModelYeldStress::read(viscosityProperties);

    fractalDimensionCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

//inicio da edição
	fractalDimensionCoeffs_.lookup("nu0") >> nu0_;
	fractalDimensionCoeffs_.lookup("En") >> En_;
	fractalDimensionCoeffs_.lookup("nuMax") >> nuMax_;
	fractalDimensionCoeffs_.lookup("CTau") >> CTau_;
	fractalDimensionCoeffs_.lookup("f0") >> f0_;
	fractalDimensionCoeffs_.lookup("D0") >> D0_;
	fractalDimensionCoeffs_.lookup("D1") >> D1_;
//final da edição

    return true;
}


// ************************************************************************* //
