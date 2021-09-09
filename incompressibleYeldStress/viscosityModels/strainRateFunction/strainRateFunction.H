/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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
    Foam::viscosityModelYeldStresss::strainRateFunction

Description
    Run-time selected strain-rate function non-Newtonian viscosity model.

    Example linear function of strain-rate:
    \verbatim
        transportModel  strainRateFunction;

        strainRateFunctionCoeffs
        {
            function polynomial ((0 0.1) (1 1.3));
        }
    \endverbatim

See also
    Foam::viscosityModelYeldStress
    Foam::Function1

SourceFiles
    strainRateFunction.C

\*---------------------------------------------------------------------------*/

#ifndef strainRateFunction_H
#define strainRateFunction_H

#include "viscosityModelYeldStress.H"
#include "volFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModelYeldStresss
{

/*---------------------------------------------------------------------------*\
                           Class strainRateFunction Declaration
\*---------------------------------------------------------------------------*/

class strainRateFunction
:
    public viscosityModelYeldStress
{
    // Private Data

        //- Coefficients dictionary
        dictionary strainRateFunctionCoeffs_;

        //- Strain-rate function
        autoPtr<Function1<scalar>> strainRateFunction_;

        //- Current viscosity field
        volScalarField nu_;
        
        //- Current viscosity field
        volScalarField tauY_;


public:

    //- Runtime type information
    TypeName("strainRateFunction");


    // Constructors

        //- Construct from components
        strainRateFunction
        (
            const word& name,
            const word& nameDois,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            //const volScalarField& T,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~strainRateFunction()
    {}


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const;
        
        //- Return the laminar viscosity
        virtual tmp<volScalarField> tauY() const;

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const;
        
        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> tauY(const label patchi) const;

        //- Correct the laminar viscosity
        virtual void correct();

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModelYeldStresss
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
