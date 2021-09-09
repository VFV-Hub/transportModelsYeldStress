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

#include "turnOnPlasticityYeldStress.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModelYeldStresss
{
    defineTypeNameAndDebug(turnOnPlasticityYeldStress, 0);

    addToRunTimeSelectionTable
    (
        viscosityModelYeldStress,
        turnOnPlasticityYeldStress,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModelYeldStresss::turnOnPlasticityYeldStress::calcNu() const
{

	dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);
	dimensionedScalar TVSMALL("VSMALL", dimTemperature, VSMALL);

    tmp<volScalarField> sr(strainRate());
	
	const volScalarField& T=U_.mesh().lookupObject<volScalarField>("T");
	
	return
	(
		min
        	(
			nuMax_,
			(tauYRef_*(
					max(
						scalar(0.0),
						exp(
							STau_*(
									scalar(1.0)/max(T, TVSMALL)-scalar(1.0)/TYS_
									)
							)
						-scalar(1.0)
						)
					)  	
			+ kRef_*(
						max(
							scalar(0.0),
							exp(
								Sk_*(
										scalar(1.0)/max(T, TVSMALL)-scalar(1.0)/TPL_
										)
								)
							-scalar(1.0)
							)
						)
					*rtone*pow(tone*sr(), n_)
			+nu0Ref_*(
						exp(
							SNu_*(
									scalar(1.0)/max(T, TVSMALL)-scalar(1.0)/TRef_
									)
							)
						)
					*sr()
			)
        	/(
				max(
					sr(), dimensionedScalar ("VSMALL", dimless/dimTime, VSMALL)
					)
		 	)
			)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::viscosityModelYeldStresss::turnOnPlasticityYeldStress::calcTauY() const
{
	dimensionedScalar TVSMALL("VSMALL", dimTemperature, VSMALL);
	const volScalarField& T=U_.mesh().lookupObject<volScalarField>("T");

	return
	(     	
		tauYRef_*(
					max(
						scalar(0.0),
						exp(
							STau_*(
									scalar(1.0)/max(T, TVSMALL)-scalar(1.0)/TYS_
									)
							)
						-scalar(1.0)
						)
					)
			  	
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModelYeldStresss::turnOnPlasticityYeldStress::turnOnPlasticityYeldStress
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
    turnOnPlasticityYeldStressCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
//inicio da edição    
	n_("n", dimless, turnOnPlasticityYeldStressCoeffs_),
	tauYRef_("tauYRef", dimViscosity/dimTime, turnOnPlasticityYeldStressCoeffs_),
	STau_("STau", dimTemperature, turnOnPlasticityYeldStressCoeffs_),
	TYS_("TYS", dimTemperature, turnOnPlasticityYeldStressCoeffs_),
	TPL_("TPL", dimTemperature, turnOnPlasticityYeldStressCoeffs_),
	kRef_("kRef", dimViscosity, turnOnPlasticityYeldStressCoeffs_),
	Sk_("Sk", dimTemperature, turnOnPlasticityYeldStressCoeffs_),
	nu0Ref_("nu0Ref", dimViscosity, turnOnPlasticityYeldStressCoeffs_),	
	SNu_("SNu", dimTemperature, turnOnPlasticityYeldStressCoeffs_),
	TRef_("TRef", dimTemperature, turnOnPlasticityYeldStressCoeffs_),
	nuMax_("nuMax", dimViscosity, turnOnPlasticityYeldStressCoeffs_),
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

bool Foam::viscosityModelYeldStresss::turnOnPlasticityYeldStress::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModelYeldStress::read(viscosityProperties);

    turnOnPlasticityYeldStressCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

//inicio da edição
	turnOnPlasticityYeldStressCoeffs_.lookup("n") >> n_;
	turnOnPlasticityYeldStressCoeffs_.lookup("tauYRef") >> tauYRef_;
	turnOnPlasticityYeldStressCoeffs_.lookup("STau") >> STau_;
	turnOnPlasticityYeldStressCoeffs_.lookup("TYS") >> TYS_;
	turnOnPlasticityYeldStressCoeffs_.lookup("TPL") >> TPL_;
	turnOnPlasticityYeldStressCoeffs_.lookup("kRef") >> kRef_;
	turnOnPlasticityYeldStressCoeffs_.lookup("Sk") >> Sk_;
	turnOnPlasticityYeldStressCoeffs_.lookup("nu0Ref") >> nu0Ref_;
	turnOnPlasticityYeldStressCoeffs_.lookup("SNu") >> SNu_;
	turnOnPlasticityYeldStressCoeffs_.lookup("TRef") >> TRef_;
	turnOnPlasticityYeldStressCoeffs_.lookup("nuMax") >> nuMax_;
//final da edição

    return true;
}


// ************************************************************************* //
