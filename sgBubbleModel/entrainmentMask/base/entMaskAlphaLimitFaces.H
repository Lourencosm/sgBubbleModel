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

Namespace
    Foam::entrainmentMask

Class
    Foam::entrainmentMask::entMaskAlphaLimitFaces

Description
    Base class for entrainment masks
    Faces selected for entrainment based on a sacalar field function approach

SourceFiles
    entMaskAlphaLimitFaces.H

Author
    Lourenço Sassetti Mendes (lourenco.sassetti@gmail.com)

\*---------------------------------------------------------------------------*/
#ifndef entrainmentMask_entMaskAlphaLimitFaces_H
#define entrainmentMask_entMaskAlphaLimitFaces_H

#include "entrainmentMask.H"

// openfoam libraries
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "faceZoneSet.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace entrainmentMasks
{

/*---------------------------------------------------------------------------*\
                           Class entrainmentMask Declaration
\*---------------------------------------------------------------------------*/


class entMaskAlphaLimitFaces
:
    public entrainmentMask
{
    // Private data
protected:
        scalar surfaceAlpha;

        volScalarField      alphaSearchMask;
        surfaceScalarField  entrainmentMask_faces_scalar;
        surfaceVectorField  alpha1_grad_faces;
        surfaceVectorField  entrainmentMask_faces_vector;
        faceSet             entSurfSet;        

  // private Functions    

        
public:

    //- Runtime type information
    TypeName("alphaLimitFaces");

    // Constructors

        //- Construct from components
        entMaskAlphaLimitFaces
        (
            const fvMesh& mesh_,
            const dictionary& bubbleModelDict_
        );

    //- Destructor
    virtual ~entMaskAlphaLimitFaces();

    // Public Member Functions

        // identify entrainment surface faces - method B - an alpha.water 1-0 mask and finding the faces in between
        //                                    - about 200x faster than the other method
        virtual void mask();

        // update the entrainmentMask
        virtual void update();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace entrainmentMask
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //