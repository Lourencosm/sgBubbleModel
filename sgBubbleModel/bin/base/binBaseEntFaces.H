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
    Foam::bubbleBin

Description
    Base class for bubble bins of type Shi et al. 2010
    Entrainment at faces

SourceFiles
    binBaseEntFaces.C

Author
    Lourenço Sassetti Mendes (lourenco.sassetti@gmail.com)

Date
    2019-02-15

\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
bubble population class for air entrainment and transport model 
belongs to interFoam_sgbm solver

\*---------------------------------------------------------------------------*/

#ifndef bins_binBaseEntFaces_H
#define bins_binBaseEntFaces_H

#include "bin.H"


// openfoam libraries
#include "addToRunTimeSelectionTable.H"
#include "linear.H"
#include "fvm.H"
//#include "fvc.H"
#include "uniformDimensionedFields.H"
#include "fvOptions.H"
#include "bound.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace bins
{

/*---------------------------------------------------------------------------*\
                           Class bubbleBin Declaration
\*---------------------------------------------------------------------------*/


class binBaseEntFaces
:
    public bin
{
    // Private data
protected:
        word name_;

        volScalarField Nb_;

        volVectorField Ub_; 
        surfaceScalarField phi_b_;
        volScalarField Db_;

        surfaceVectorField En_faces_;
        volScalarField En_;
        volScalarField IG_;

        volScalarField Nb_evap_;
        volScalarField Vb_;
        volScalarField Vb_evap_;
        volScalarField C_;
        volScalarField C_evap_;

        /*
        Shi et al 2010 - model variables
        */
        dimensionedScalar a_b_;
        dimensionedScalar Pr0_threshold_;
        
        scalar Cb_max_;
        scalar EnAlphaMin_;
        scalar Sch_;

        dimensionedScalar r_b_;
        dimensionedScalar dr_b_;
        dimensionedScalar VOneBubble_;
        dimensionedScalar Nb_min_;
        dimensionedScalar Di_b_;

        fv::options& fvOptions;
        
public:

    //- Runtime type information
    TypeName("binBaseEntFaces");

    // Constructors

        //- Construct from components
        binBaseEntFaces
        (
            const fvMesh& mesh_,
            const dictionary& bubbleBinDict_   
        );

    // Destructor
        virtual ~binBaseEntFaces();


    // Member Functions

        virtual const word& name() const;

        //- Return const-access to bubbleBin diffusion coefficient
        virtual const dimensionedScalar& a_b() const;

        //- Return const-access to bubbleBin diffusion coefficient
        virtual const dimensionedScalar& P_r0_threshold() const;

        //- Return const-access to bubbleBin diffusion coefficient
        virtual const scalar& EnAlphaMin() const;

        //- Return const-access to bubbleBin radius
        virtual const dimensionedScalar& r_b() const;

        //- Return const-access to bubbleBin averaged bubble radius spacing
        virtual const dimensionedScalar& dr_b() const;

        //- Return const-access to bubbleBin averaged bubble radius spacing
        virtual const dimensionedScalar& VOneBubble() const;

        // field
        virtual const volScalarField& C() const;

        // field
        virtual const volScalarField& C_evap() const;

        // field
        virtual const volScalarField& Nb() const;

        // field
        virtual const volScalarField& Nb_evap() const;
 
        // field
        virtual const volScalarField& Vb() const;
 
        // field
        virtual const volScalarField& Vb_evap() const;

        // field
        virtual const volScalarField& En() const;

        // field
        virtual const surfaceVectorField& En_faces() const;

        // field
        virtual const volVectorField& Ub() const;

        // field
        virtual const surfaceScalarField& phi_b() const;

        // air bubbles turbulent diffusivity
        virtual const volScalarField& Db() const;

        // term for bubble intergroup transfer in the number of bubble eq.
        virtual const volScalarField& IG() const;

        // update bubble velocity
        virtual void update_Ub();

        // calculate air entranment in faces - method A
        virtual void entrainment();

        // solve bubble transport equation - number of bubbles density
        virtual void solve_Nb();
        
        // calculate bubble volume and bubble fraction       
        virtual void update_Vb_C();

        // bounding bubble bins fields by a factor in each cell
        virtual void Cb_bounding(const volScalarField& Cb_bounder);        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace bins
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
