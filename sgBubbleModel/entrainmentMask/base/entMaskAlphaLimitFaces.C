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
    entMaskAlphaLimitFaces.H

Author
    Lourenço Sassetti Mendes (lourenco.sassetti@gmail.com)

Date
    2019-02-14

\*---------------------------------------------------------------------------*/

#include "entMaskAlphaLimitFaces.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace entrainmentMasks
{
    defineTypeNameAndDebug(entMaskAlphaLimitFaces, 0);
    addToRunTimeSelectionTable(entrainmentMask, entMaskAlphaLimitFaces, dictionary);    
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::entrainmentMasks::entMaskAlphaLimitFaces::entMaskAlphaLimitFaces
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
    InfoPos("\n<---entrainmentMasks::entMaskAlphaLimitFaces::entMaskAlphaLimitFaces(): start"); 

    update();

    InfoPos("\n<---entrainmentMasks::entMaskAlphaLimitFaces::entMaskAlphaLimitFaces(): end");    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::entrainmentMasks::entMaskAlphaLimitFaces::~entMaskAlphaLimitFaces()
{}

// * * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * //  


// identify entrainment surface faces - method B - an alpha.water 1-0 mask and finding the faces in between
//                                    - about 200x faster than the other method
void Foam::entrainmentMasks::entMaskAlphaLimitFaces::mask()
{
//    clockTime time_ = clockTime();
    const volScalarField& alpha_ = mesh_.lookupObject<volScalarField>("alpha.water");

    faceSet& surfacesSet = entSurfSet;

    //******* update base mask ********//
    alphaSearchMask = pos(surfaceAlpha - alpha_);
    
    // create entrainmentSurfaceFieldMask to determine where entrainment takes place (@faces)
    entrainmentMask_faces_scalar = pos(0.4-mag(linearInterpolate(alphaSearchMask)-0.5));

    /******* select faces that match the surface the best - only internal faces ********/

    // reset faceSeZoneSet of the alpha.water surface faces
    surfacesSet.deleteSet(surfacesSet);

    // check if it is a faceSet
    if (!isA<faceSet>(surfacesSet))
    {
        FatalError
        << "bubbleEntMask----> " << surfacesSet.name() <<" <-  faceSet doesn't exist!"
        << "    aborting."
        << exit(FatalError);
    }
    else
    {
        // ~~~~ Select entrainment surface faces and add to faceSet
        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            if (entrainmentMask_faces_scalar[facei]>0)
            {
                surfacesSet.insert(facei);
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// update the entrainmentMask
void Foam::entrainmentMasks::entMaskAlphaLimitFaces::update()
{
    InfoPos("\n<---entrainmentMasks::entMaskAlphaLimitFaces::update()---->air entrainment mask @ faces -> method (B): alpha field mask: start");

    const volScalarField& alpha1_ = mesh_.lookupObject<volScalarField>("alpha.water");

    // identify entrainment surface faces - using the an alphaLimit maks - method B - 200x+ faster
    // obtion b:
    mask();

    // gradient of alpha @faces
    alpha1_grad_faces = linearInterpolate(fvc::grad(alpha1_));

    // mag of alpha gradient vectors
    dimensionedScalar smallNum("smallNum", alpha1_grad_faces.dimensions(), 1.112112e-15);
    
    tmp<surfaceScalarField> fieldMagLimited = max(mag(alpha1_grad_faces),smallNum);

    // ~~~ entrainment surface mask with alpha gradient unitary vectors
    entrainmentMask_faces_vector = entrainmentMask_faces_scalar*alpha1_grad_faces/fieldMagLimited;

    InfoPos("\n---->entrainmentMasks::entMaskAlphaLimitFaces::update()---->air entrainment mask @ faces -> method (B): alpha field mask: end");
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
