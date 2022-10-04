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
    Foam::bubbleModels::shiNext

Description
    Base class for bubbleModels of type Shi et al. 2010:
        - initializes and solves a bubble transport model that runs over InterFoam solver    

SourceFiles
    shiNext.C

Author
    Lourenço Sassetti Mendes (lourenco.sassetti@gmail.com)

Date
    2019-10-08

\*---------------------------------------------------------------------------*/

#ifndef bubbleModels_shiNext_H
#define bubbleModels_shiNext_H

#include "shiBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace bubbleModels
{

/*---------------------------------------------------------------------------*\
                          Class shiNext Declaration
\*---------------------------------------------------------------------------*/

class shiNext
:
    public shiBase
{
protected:
    // Protected data
        List<word> availBubbleBins;
        List<word> availEntrainmentMasks;
        List<word> availIntergroupTransfers;

        word binType;
        word entrainmentMaskType;
        word intergroupTransferType;  

//        word Pr_rhoName;

        volScalarField magSqrS;
        volScalarField Pr;

    // Protected Member Functions
 

public:

    //- Runtime type information
    TypeName("shiNext");

    // Constructor

    //- Construct from components
    shiNext
    (
        const fvMesh& mesh
    );

    //- Destructor
    virtual ~shiNext();

        // check the model components types considering the available model components types for this bubble model
        virtual void check_model_components();
        
        // calculate some flow fields needed for the bubble entrainment/transport
        virtual void update_auxFields();

        // sum all bubble bins: number of bubbles, volume of bubbles, entrained bubbles, bubble fraction, evaporated bubble fraction
        virtual void sumBins();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace bubbleModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
