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
    Foam::entrainmentMask::entMaskIscHits

Description
    Base class for entrainment masks
    Faces selected for bubble entrainment based on the IsoSurfaceCell method

SourceFiles
    entMaskIscHits.H

Author
    Lourenço Sassetti Mendes (lourenco.sassetti@gmail.com)

\*---------------------------------------------------------------------------*/
#ifndef entrainmentMask_entMaskIscHits_H
#define entrainmentMask_entMaskIscHits_H

#include "entrainmentMask.H"

// openfoam libraries
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "volPointInterpolation.H"
#include "isoSurface.H"
#include "isoSurfaceCell.H"
#include "triSurfaceMesh.H"
#include "faceZoneSet.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace entrainmentMasks
{

/*---------------------------------------------------------------------------*\
                           Class entrainmentMask Declaration
\*---------------------------------------------------------------------------*/


class entMaskIscHits
:
    public entrainmentMask
{
    // Private data
protected:
        scalar surfaceAlpha;

        volScalarField      alphaSearchMask;
        surfaceScalarField  entrainmentMask_faces_scalar;
//        volVectorField      alpha1_grad_cells;
        surfaceVectorField  alpha1_grad_faces;
        surfaceVectorField  entrainmentMask_faces_vector;
        faceSet             entSurfSet;        

  // private Functions    

        
public:

    //- Runtime type information
    TypeName("iscHits");

    // Constructors

        //- Construct from components
        entMaskIscHits
        (
            const fvMesh& mesh_,
            const dictionary& bubbleModelDict_
        );

    //- Destructor
    virtual ~entMaskIscHits();

    // Public Member Functions
    
        // get unique and sorted values from labelList
        virtual labelList uniqueOrderList(const labelList& labelList_);

        // create isoSurfaceCell of alpha.water field withe default 0.5 iso
        virtual isoSurfaceCell alphaWater_isoSurfaceCell(const fvMesh& mesh_, const scalar alphaIso = 0.5);

        // update cell mask field - using a list of cell indexes
        virtual void cells_mask(const labelList& cellIdxlabelList, volScalarField& cellsMaskField);

        // identify entrainment surface faces - method A - using the isoSurfaceCell object class and a surface hits detection method 
        //                                    - in this method only isoSurfaceCell Meshed cells are seached for surface intersection: it's faster than all mesh cells
        //                                    - inspired in: Foam::searchableSurfaceToFaceZone::applyToSet
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