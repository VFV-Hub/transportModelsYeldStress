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

Class
    Foam::viscosityModels::fractalDimension

Description
     Standard power-law non-Newtonian viscosity model.

SourceFiles
    fractalDimension.C

\*---------------------------------------------------------------------------*/

#ifndef fractalDimension_H
#define fractalDimension_H

#include "viscosityModelYeldStress.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModelYeldStresss
{

/*---------------------------------------------------------------------------*\
                           Class fractalDimension Declaration
\*---------------------------------------------------------------------------*/

class fractalDimension
:
    public viscosityModelYeldStress
{
    // Private data

        dictionary fractalDimensionCoeffs_;

//inicio da edição
		dimensionedScalar nu0_;
		dimensionedScalar En_;
		dimensionedScalar nuMax_;
		dimensionedScalar CTau_;
		dimensionedScalar f0_;
		dimensionedScalar D0_;
		dimensionedScalar D1_;		


        volScalarField nu_;
        volScalarField tauY_;

//final da edição


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() const;
        tmp<volScalarField> calcTauY() const;

public:

    //- Runtime type information
    TypeName("fractalDimension");


    // Constructors

        //- Construct from components
        fractalDimension
        (
            const word& name,
            const word& nameDois,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const volScalarField& T,
            const surfaceScalarField& phi
        );


    //- Destructor
    ~fractalDimension()
    {}


    // Member Functions

        //- Return the laminar viscosity
        tmp<volScalarField> nu() const
        {
            return nu_;
        }
        tmp<volScalarField> tauY() const
        {
            return tauY_;
        }
        
        //- Return the laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }
        
        //- Return the laminar viscosity for patch
        tmp<scalarField> tauY(const label patchi) const
        {
            return tauY_.boundaryField()[patchi];
        }

        //- Correct the laminar viscosity
        void correct()
        {
        	tauY_ = calcTauY();
            nu_ = calcNu();
        }

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModelYeldStresss
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
