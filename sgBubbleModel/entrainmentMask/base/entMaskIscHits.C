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
    Foam::entrainmentMasks

Description
    air bubble entrainment models: entrainment mask
    belongs to interFoam_sgbm solver version 46 and forward

SourceFiles
    entMaskIscHits.H

Author
    Lourenço Sassetti Mendes (lourenco.sassetti@gmail.com)

Date
    2019-02-14

\*---------------------------------------------------------------------------*/

#include "entMaskIscHits.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace entrainmentMasks
{
    defineTypeNameAndDebug(entMaskIscHits, 0);
    addToRunTimeSelectionTable(entrainmentMask, entMaskIscHits, dictionary);    
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::entrainmentMasks::entMaskIscHits::entMaskIscHits
(
    const fvMesh& mesh,
    const dictionary& bubbleModelDict
)
:
    // base class
    entrainmentMask(mesh, bubbleModelDict),
    // entrainment surface alpha value
    surfaceAlpha(readScalar(bubbleModelDict.lookup("surfaceAlpha"))),
    // air entrainment surface search mask at cells
    alphaSearchMask
    (
        IOobject
        (
            "alphaSearchMask",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    // air entrainment production mask at faces 1-0
    entrainmentMask_faces_scalar
    (
        IOobject
        (
            "entrainmentMask_faces_scalar",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    // gradient of alpha @ faces
    alpha1_grad_faces
    (
        IOobject
        (
            "alpha.water_grad_faces",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
//        linearInterpolate(alpha1_grad_cells)
        linearInterpolate(fvc::grad(mesh.lookupObject<volScalarField>("alpha.water")))
    ),
    // air entrainment surface mask with alpha gradient unitary vectors
    entrainmentMask_faces_vector
    (
        IOobject
        (
            "entrainmentMask_faces_vector",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimless, vector( 0,0,0 ))    
    ),
    entSurfSet
    (
        IOobject
        (
            "faceSet_entSurf",
            mesh.time().timeName()+"/polyMesh/sets",    
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    )
{
    InfoPos("\n<---entrainmentMasks::entMaskIscHits::entMaskIscHits(): start"); 

    update();

    InfoPos("\n<---entrainmentMasks::entMaskIscHits::entMaskIscHits(): end");    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::entrainmentMasks::entMaskIscHits::~entMaskIscHits()
{}

// * * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * //  



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// get unique and sorted values from labelList
labelList Foam::entrainmentMasks::entMaskIscHits::uniqueOrderList(const labelList& labelList_)
{

    labelList uniqueIdx;
    Foam::uniqueOrder(labelList_, uniqueIdx);
    
    // create unique indexes list = remove repeated indexes 
    DynamicList<label> uniqueSortedList;
    forAll(uniqueIdx, i)        // loop over all cells in meshed cells                
    {
        uniqueSortedList.append(labelList_[uniqueIdx[i]]);
    }
    return uniqueSortedList;
}

// create isoSurfaceCell of alpha.water field withe default 0.5 iso
isoSurfaceCell Foam::entrainmentMasks::entMaskIscHits::alphaWater_isoSurfaceCell(const fvMesh& mesh_, const scalar alphaIso)
{
    //******* create isoSurfaceCell ********//
    const volScalarField& alphaWater_ = mesh_.lookupObject<volScalarField>("alpha.water");
    const bool regularize_ = true; //false //true

    // volScalarField "alphaWater_" values at cells centers
    const scalarField& alphaWater_cells = alphaWater_.internalField();

    // interpolate the volScalarField "alphaWater_" to get values at points
    volPointInterpolation pMesh(mesh_);
    pointScalarField alphaWaterPTS_(pMesh.interpolate(alphaWater_)());

    // create IsoSurfaceCell
    isoSurfaceCell isoSurfCell_(mesh_, alphaWater_cells, alphaWaterPTS_, alphaIso, regularize_);

    return isoSurfCell_;
}


// update cell mask field - using a list of cell indexes
void Foam::entrainmentMasks::entMaskIscHits::cells_mask(const labelList& cellIdxlabelList, volScalarField& cellsMaskField)
{
    // crete unique indexes list = remove repeated indexes
    labelList uniqueIdx;
    Foam::uniqueOrder(cellIdxlabelList, uniqueIdx);

    // re-define entrainment cells mask field
    cellsMaskField = 0.0;

    forAll(uniqueIdx, i)
    {
        // change value of entrainmentMask_ volScalarField
        cellsMaskField[cellIdxlabelList[uniqueIdx[i]]] = 1.0;
    } 

}


// identify entrainment surface faces - method A - using the isoSurfaceCell object class and a surface hits detection method 
//                                    - in this method only isoSurfaceCell Meshed cells are seached for surface intersection: it's faster than all mesh_ cells
//                                    - inspired in: Foam::searchableSurfaceToFaceZone::applyToSet
void Foam::entrainmentMasks::entMaskIscHits::mask()
{
    //******* create isoSurfaceCell ********//
    isoSurfaceCell isoSurfCell_ = alphaWater_isoSurfaceCell(mesh_, surfaceAlpha);
    
    // write isoSurfaceCell meshed cells to volScalarField
    cells_mask(isoSurfCell_.meshCells(), alphaSearchMask);


    /******* isoSurface to triSurface ********/
    isoSurfCell_.triangulate();

//v1712    pointField  isc_pts(isoSurfCell_.xferPoints ());
//v1712    List<face> isc_faces(isoSurfCell_.xferFaces ()); 

//v1806
    const List<face> isc_faces = isoSurfCell_.surfFaces(); 
    const pointField  isc_pts = isoSurfCell_.points();

    // create triFaceList
    List<labelledTri> triFaces(isc_faces.size());
    forAll(triFaces, facei)
    {
        const face& f = isc_faces[facei];
        labelledTri& tri = triFaces[facei];

        tri[0] = f[0];
        tri[1] = f[1];
        tri[2] = f[2];
        tri.region() = 0;
    }

    // create a triSurfaceMesh
    triSurfaceMesh surfMesh_
    (
        IOobject
        (
            "isoSurfaceCell.vtk",
            mesh_.time().timeName(),
            mesh_.objectRegistry::db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),        
        triSurface(triFaces, isc_pts) 
    );

    /******* isoSurface meshed cells ********/
    labelList meshedCellsIdx = uniqueOrderList(isoSurfCell_.meshCells());

    /******* isoSurface meshed cells all faces ********/
    DynamicList<label> meshedFacesIdx;

    // loop over all cells in meshed cells and get the faces
    forAll(meshedCellsIdx, n)                        
    {
        const label& cellIdx = meshedCellsIdx[n];
        const cell& cell_ = mesh_.cells()[cellIdx];
        forAll(cell_, i)        // loop over all faces in cellID
        {
            meshedFacesIdx.append(cell_[i]);
        }
    }

    labelList uniqueIdx_f;    
    Foam::uniqueOrder(meshedFacesIdx, uniqueIdx_f);

    // loop over all faces of meshed cells: remove duplicate and boundary faces
    DynamicList<label> meshedFaces;
    forAll(uniqueIdx_f, i)                        
    {
        label facei = meshedFacesIdx[uniqueIdx_f[i]];
        if (facei < mesh_.nInternalFaces())
        {
            meshedFaces.append(facei);
        }
    }

    /******* simplify isoSurface merging points********/
    const scalar   mergeTol = 1e-3;

    Info<< "        merging isoSurfaceCell points within " << mergeTol << " metre." << endl;
    while (true)
    {
        label nOldVert = surfMesh_.nPoints();
        triSurfaceTools::mergePoints(surfMesh_, mergeTol);
        if (nOldVert == surfMesh_.nPoints())
        {
            break;
        }
    }

    /******* select faces that match the iso surface the best - only internal faces ********/

    // create faceSeZoneSet of the alpha.water surface faces
    faceSet& surfacesSet = entSurfSet;
    surfacesSet.deleteSet(surfacesSet);

    // check if it is a faceSet
    if (!isA<faceSet>(surfacesSet))
    {
        FatalError
        << "bubbleEntMask----> entSurfSet <-  faceSet doesn't exist!"
        << "    aborting."
        << exit(FatalError);
    }
    else
    {
        // ~~~~ Get faces cell-cell centre vectors ~~~~
        pointField start(meshedFaces.size());
        pointField end(meshedFaces.size());

        const pointField& cc = mesh_.cellCentres();

        forAll(meshedFaces, i)
        {
            const label& facei = meshedFaces[i];
            start[i] = cc[mesh_.faceOwner()[facei]];
            end[i] = cc[mesh_.faceNeighbour()[facei]];
        }

        // ~~~~ Do all intersection tests ~~~~ 
        List<pointIndexHit> hits;
        surfMesh_.findLine(start, end, hits); // note: findLine takes a lot of time

        // ~~~~ Select intersected faces
        forAll(hits, i)
        {
            if (hits[i].hit())
            {
                const label& facei = meshedFaces[i];
                surfacesSet.insert(facei);
                // create entrainmentSurfaceFieldMask to determine where entrainment takes place (@faces)
                entrainmentMask_faces_scalar[facei]=1.0;
            }
        }
    }


}


// update the entrainmentMask
void Foam::entrainmentMasks::entMaskIscHits::update()
{
    InfoPos("\n<---entrainmentMasks::entMaskIscHits::update()---->air entrainment mask @ faces -> method (B): alpha field mask: start");

    const volScalarField& alpha1_ = mesh_.lookupObject<volScalarField>("alpha.water");

    // identify entrainment surface faces - using the isoSurfaceCell object class + select surface cut cells - method A
    // option a:
    mask();

    // gradient of alpha @faces
    alpha1_grad_faces = linearInterpolate(fvc::grad(alpha1_));

    // mag of alpha gradient vectors
    dimensionedScalar smallNum("smallNum", alpha1_grad_faces.dimensions(), 1.112112e-15);
    
    surfaceScalarField fieldMagLimited = max(mag(alpha1_grad_faces),smallNum);

    // ~~~ entrainment surface mask with alpha gradient unitary vectors
    entrainmentMask_faces_vector = entrainmentMask_faces_scalar*alpha1_grad_faces/fieldMagLimited;

    InfoPos("\n---->entrainmentMasks::entMaskIscHits::update()---->air entrainment mask @ faces -> method (B): alpha field mask: end");
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
