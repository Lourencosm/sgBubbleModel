/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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
    Foam::bin

Description
    Base class for bubble bins

SourceFiles
    bin.C    
    binNew.C

Author
    Lourenço Sassetti Mendes (lourenco.sassetti@gmail.com)


\*---------------------------------------------------------------------------*/

#ifndef bin_H
#define bin_H

// openfoam libraries
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fvMesh.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class bin Declaration
\*---------------------------------------------------------------------------*/

class bin
{

protected:
        const fvMesh& mesh_;
        const dictionary& bubbleBinDict_;

        word binType_;

        bool debugMode;
        bool debugInfoPos;
        bool debugInfoAux;
        bool debugInfoVar;        
        bool debugInfoCalc; 
    // Protected Member Functions

private:    

    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        bin(const bin&);

        //- Disallow default bitwise assignment
        void operator=(const bin&);


public:

    //- Runtime type information
    TypeName("bin");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        bin,
        dictionary,
        (
            const fvMesh& mesh_,
            const dictionary& bubbleBinDict_            
        ),
        (mesh_, bubbleBinDict_)
    );

    // Constructors

    //- Construct from components
    bin
    (
        const fvMesh& mesh_,
        const dictionary& bubbleBinDict_         
    );

    //- Destructor
    virtual ~bin();

    // Selectors

        //- Return a reference to the selected bubble model
        static autoPtr<bin> New
        (
            const fvMesh& mesh_,
            const dictionary& bubbleBinDict_            
        );

    // Public Member Functions

        //- Return const-access to bin type
        virtual const word& binType() const
        {
            return binType_;
        }

        // display bin name, activation status and debug mode
        virtual void bin_status();

        // Info messages - code position
        void InfoPos(const word& msg);

        // Info messages - auxiliary info
        void InfoAux(const word& msg);

        // Info messages - variable info: min,max,avg, dimensions
        void InfoVar(const volScalarField& sf);
        void InfoVar(const surfaceScalarField& sf);

        virtual const word& name() const = 0;

        const dictionary& dict() const
        {
            return bubbleBinDict_;
        }

        //- Return const-access to bubbleBin diffusion coefficient
        virtual const dimensionedScalar& a_b() const =  0;

        //- Return const-access to bubbleBin diffusion coefficient
        virtual const dimensionedScalar& P_r0_threshold() const =  0;

        //- Return const-access to bubbleBin diffusion coefficient
        virtual const scalar& EnAlphaMin() const =  0;

        //- Return const-access to bubbleBin radius
        virtual const dimensionedScalar& r_b() const =  0;

        //- Return const-access to bubbleBin averaged bubble radius spacing
        virtual const dimensionedScalar& dr_b() const =  0;

        //- Return const-access to bubbleBin averaged bubble radius spacing
        virtual const dimensionedScalar& VOneBubble() const =  0;

        virtual const volScalarField& C() const =  0;

        // field
        virtual const volScalarField& C_evap() const =  0;

        // field
        virtual const volScalarField& Nb() const =  0;
        // field
        virtual const volScalarField& Nb_evap() const =  0;
 
        // field
        virtual const volScalarField& Vb() const =  0;
 
        // field
        virtual const volScalarField& Vb_evap() const =  0;

        // field
        virtual const volScalarField& En() const =  0;

        // field
        virtual const surfaceVectorField& En_faces() const =  0;

        // field
        virtual const volVectorField& Ub() const =  0;

        // field
        virtual const surfaceScalarField& phi_b() const =  0;

        // air bubbles turbulent diffusivity
        virtual const volScalarField& Db() const =  0;

        // term for bubble intergroup transfer in the number of bubble eq.
        virtual const volScalarField& IG() const =  0;

        // update bubble velocity
        virtual void update_Ub() = 0;

        // calculate air entranment in faces - method A
        virtual void entrainment() = 0;

        // solve bubble transport equation - number of bubbles density
        virtual void solve_Nb() = 0;
        
        // calculate bubble volume and bubble fraction       
        virtual void update_Vb_C() = 0;

        // bounding bubble bins fields by a factor in each cell
        virtual void Cb_bounding(const volScalarField& Cb_bounder) = 0;               
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
