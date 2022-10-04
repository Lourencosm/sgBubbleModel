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
    Base class for entrainment masks
    Faces selected for entrainment based on interface cells between min alpha and max alpha

SourceFiles
    interfaceCells.H

Author
    Lourenço Sassetti Mendes (lourenco.sassetti@gmail.com)

Date
    2019-04-26

\*---------------------------------------------------------------------------*/

#include "interfaceCells.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace entrainmentMasks
{
    defineTypeNameAndDebug(interfaceCells, 0);
    addToRunTimeSelectionTable(entrainmentMask, interfaceCells, dictionary);    
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::entrainmentMasks::interfaceCells::interfaceCells
(
    const fvMesh& mesh,
    const dictionary& bubbleModelDict
)
:
    // base class
    entrainmentMask(mesh, bubbleModelDict),
    // entrainment surface min alpha value
    alphaMin(readScalar(bubbleModelDict.lookup("interfaceAlphaMin"))),
    // entrainment surface max alpha value
    alphaMax(readScalar(bubbleModelDict.lookup("interfaceAlphaMax"))),    

    // air entrainment production mask at cells 1-0
    entrainmentMask_cells_scalar
    (
        IOobject
        (
            "entrainmentMask_cells_scalar",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    )
{
    InfoPos("\n<---entrainmentMasks::interfaceCells::interfaceCells(): start"); 

    update();

    InfoPos("\n<---entrainmentMasks::interfaceCells::interfaceCells(): end");    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::entrainmentMasks::interfaceCells::~interfaceCells()
{}

// * * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * //  


// identify entrainment surface cells - method B - an alpha.water 1-0 mask and finding the faces in between
//                                    - about 200x faster than the other method
void Foam::entrainmentMasks::interfaceCells::mask()
{
//    clockTime time_ = clockTime();
    const volScalarField& alpha_ = mesh_.lookupObject<volScalarField>("alpha.water");

    // reset the entrainment faces mask
    //entrainmentMask_cells_scalar = 0.0;

    //******* update base mask ********//
    // create entrainment cell FieldMask to determine where entrainment takes place (@cells)
    entrainmentMask_cells_scalar = pos(alphaMax - alpha_)*pos(alpha_ - alphaMin);    
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// update the entrainmentMask
void Foam::entrainmentMasks::interfaceCells::update()
{
    InfoPos("\n<---entrainmentMasks::interfaceCells::entMask_update()---->air entrainment mask @ cells -> method (B): alpha field mask: start");

    // identify entrainment surface cells - using the an alphaLimit maks - method B - 200x+ faster
    // obtion b:
    mask();

    InfoPos("\n<---entrainmentMasks::interfaceCells::entMask_update()---->air entrainment mask @ cells -> method (B): alpha field mask: end");
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
