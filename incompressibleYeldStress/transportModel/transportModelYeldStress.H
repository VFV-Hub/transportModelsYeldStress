/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::transportModelYeldStress

Description
    Base-class for all transport models used by the incompressible turbulence
    models.

SourceFiles
    transportModelYeldStress.C

\*---------------------------------------------------------------------------*/

#ifndef transportModelYeldStress_H
#define transportModelYeldStress_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class transportModelYeldStress Declaration
\*---------------------------------------------------------------------------*/

class transportModelYeldStress
{
public:

    //- Runtime type information
    TypeName("transportModelYeldStress");


    // Constructors

        //- Construct from components
        transportModelYeldStress
        ();

        //- Disallow default bitwise copy construction
        transportModelYeldStress(const transportModelYeldStress&);


    //- Destructor
    virtual ~transportModelYeldStress();


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const = 0;
        
        virtual tmp<volScalarField> tauY() const = 0;

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const = 0;
        
        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> tauY(const label patchi) const = 0;

        //- Correct the laminar viscosity
        virtual void correct() = 0;

        //- Read transportProperties dictionary
        virtual bool read() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const transportModelYeldStress&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
