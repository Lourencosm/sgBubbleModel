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
    Foam::bins::binBaseEntCells

Description
    Base class for bubble bins of type Shi et al. 2010
    Entrainment at cells center

SourceFiles
    binBaseEntCells.C

Author
    Lourenço Sassetti Mendes (lourenco.sassetti@gmail.com)

Date
    2019-02-15

\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
bubble population class for air entrainment and transport model 
belongs to interFoam_sgbm solver

\*---------------------------------------------------------------------------*/

#ifndef bins_binBaseEntCells_H
#define bins_binBaseEntCells_H

#include "binBaseEntFaces.H"


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


class binBaseEntCells
:
    public binBaseEntFaces
{
    // Private data
protected:
        
public:

    //- Runtime type information
    TypeName("binBaseEntCells");

    // Constructors

        //- Construct from components
        binBaseEntCells
        (
            const fvMesh& mesh_,
            const dictionary& bubbleBinDict_   
        );

    // Destructor
        virtual ~binBaseEntCells();


    // Member Functions

        // calculate air entranment in faces - method B - cells
        virtual void entrainment();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace bins
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
