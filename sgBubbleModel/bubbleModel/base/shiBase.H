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
    Foam::bubbleModels::shiBase

Description
    Base class for bubbleModels of type Shi et al. 2010:
        - initializes and solves a bubble transport model that runs over InterFoam solver    

SourceFiles
    shiBase.C

Author
    Lourenço Sassetti Mendes (lourenco.sassetti@gmail.com)

Date
    2019-02-13

\*---------------------------------------------------------------------------*/

#ifndef bubbleModels_shiBase_H
#define bubbleModels_shiBase_H

#include "bubbleModel.H"
#include "bin.H"
#include "entrainmentMask.H"
#include "intergroupTransfer.H"

// openfoam libraries
#include "addToRunTimeSelectionTable.H"
#include "linear.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace bubbleModels
{

/*---------------------------------------------------------------------------*\
                          Class shiBase Declaration
\*---------------------------------------------------------------------------*/

class shiBase
:
    public bubbleModel
{
protected:
    // Protected data

        volScalarField Vcell;
        volScalarField Nb;
        volScalarField Nb_evap;
        volScalarField Cb;
        volScalarField Cb_evap;
        volScalarField En;
        volScalarField Vb;
        volScalarField Vb_evap;
        volScalarField rb_avg;        

        scalarList r_bins;
        scalarList dr_bins;
        scalarList DN_bins;
        scalarList D_bins;

        scalar Cb_max_;

        int n_bins;

    
        // bubble entrainment mask class
        autoPtr<entrainmentMask> bubbleEntMask;

        // list of bubble bins names
        wordList bubbleBinsNameList;
        //  list of pointers to bubbleBin class having the size of the number of
        //  bubble bins retrived from bubbleProperties dict @constant folder        
        PtrList<bin> bubbleBinsList;

        // bubble intergroup transfers class
        autoPtr<intergroupTransfer> bubbleInteractions;

    // Protected Member Functions
 

public:
    // Constructor

    //- Construct from components
    shiBase
    (
        const fvMesh& mesh
    );

    //- Destructor
    virtual ~shiBase();

        // check the model components types considering the available model components types for this bubble model
        virtual void check_model_components() = 0;

        // calculate radius spacing of each bin
        virtual void calc_dr_bins();

        // calculate bubble size probability function
        virtual void calc_D_bins();

        // calculate radius spacing and bubble size probability function  of each bin
        virtual void process_bins();

        // display model variables and parameters
        virtual void model_info();

        // initialize fields and variables needed for the bubble entrainment/transport
        void init_fields();

        // calculate some flow fields needed for the bubble entrainment/transport
        virtual void update_auxFields() = 0;

        // initialize bubble bins
        virtual void init_bubbleBins();

        // intergroup bubble source/sink term for the bubbles of air transport equation
        virtual void update_IG();

        // solve all bubble bins transport equation - number of bubbles density
        virtual void solve();

        // sum all bubble bins: number of bubbles, volume of bubbles, entrained bubbles, bubble fraction, evaporated bubble fraction
        virtual void sumBins() = 0;

        // runtime cycle tasks: update_auxFields() + solve() + sumBins()
        virtual void reCalc();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace bubbleModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
