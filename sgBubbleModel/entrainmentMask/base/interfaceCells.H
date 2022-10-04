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
    Foam::entrainmentMask::interfaceCells

Description
    Base class for entrainment masks
    Faces selected for entrainment based on interface cells between min alpha and max alpha

SourceFiles
    interfaceCells.H

Author
    Lourenço Sassetti Mendes (lourenco.sassetti@gmail.com)

\*---------------------------------------------------------------------------*/
#ifndef entrainmentMask_interfaceCells_H
#define entrainmentMask_interfaceCells_H

#include "entrainmentMask.H"

// openfoam libraries
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace entrainmentMasks
{

/*---------------------------------------------------------------------------*\
                           Class entrainmentMask Declaration
\*---------------------------------------------------------------------------*/


class interfaceCells
:
    public entrainmentMask
{
    // Private data
protected:
        scalar alphaMin;
        scalar alphaMax;

        volScalarField      entrainmentMask_cells_scalar;

  // private Functions    
        
public:

    //- Runtime type information
    TypeName("interfaceCells");

    // Constructors

        //- Construct from components
        interfaceCells
        (
            const fvMesh& mesh_,
            const dictionary& bubbleModelDict_
        );

    //- Destructor
    virtual ~interfaceCells();

    // Public Member Functions

        // identify entrainment surface faces - method B - an alpha.water 1-0 mask and finding the faces in between
        //                                    - about 200x faster than the isosurface hits other method
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