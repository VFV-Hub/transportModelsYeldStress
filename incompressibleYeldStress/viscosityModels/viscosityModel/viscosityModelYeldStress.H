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

Namespace
    Foam::viscosityModelYeldStresss

Description
    A namespace for various incompressible viscosityModelYeldStress implementations.

Class
    Foam::viscosityModelYeldStress

Description
    An abstract base class for incompressible viscosityModelYeldStresss.

    The strain rate is defined by:

        mag(symm(grad(U)))


SourceFiles
    viscosityModelYeldStress.C
    viscosityModelYeldStressNew.C

\*---------------------------------------------------------------------------*/

#ifndef viscosityModelYeldStress_H
#define viscosityModelYeldStress_H

#include "dictionary.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class viscosityModelYeldStress Declaration
\*---------------------------------------------------------------------------*/

class viscosityModelYeldStress
{

protected:

    // Protected data

        word name_;
        word nameDois_;
        dictionary viscosityProperties_;

        const volVectorField& U_;
        //const volScalarField& T_;
        const surfaceScalarField& phi_;


public:

    //- Runtime type information
    TypeName("viscosityModelYeldStress");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            viscosityModelYeldStress,
            dictionary,
            (
                const word& name,
                const word& nameDois,
                const dictionary& viscosityProperties,
                const volVectorField& U,
                //const volScalarField& T,
                const surfaceScalarField& phi
            ),
            (name, nameDois, viscosityProperties, U, phi)
        );


    // Constructors

        //- Construct from components
        viscosityModelYeldStress
        (
            const word& name,
            const word& nameDois,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            //const volScalarField& T,
            const surfaceScalarField& phi
        );

        //- Disallow default bitwise copy construction
        viscosityModelYeldStress(const viscosityModelYeldStress&);


    // Selectors

        //- Return a reference to the selected viscosity model
        static autoPtr<viscosityModelYeldStress> New
        (
            const word& name,
            const word& nameDois,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            //const volScalarField& T,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~viscosityModelYeldStress()
    {}


    // Member Functions

        //- Return the phase transport properties dictionary
        const dictionary& viscosityProperties() const
        {
            return viscosityProperties_;
        }

        //- Return the strain rate
        tmp<volScalarField> strainRate() const;
        
        //tmp<volScalarField> dTdt() const;

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const = 0;
        
        //- Return the laminar viscosity
        virtual tmp<volScalarField> tauY() const = 0;

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const = 0;
        
        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> tauY(const label patchi) const = 0;

        //- Correct the laminar viscosity
        virtual void correct() = 0;

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties) = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const viscosityModelYeldStress&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
