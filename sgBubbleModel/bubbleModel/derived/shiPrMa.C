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
    Foam::bubbleModels::shiPrMa

Description
    Base class for bubbleModels of type Shi et al. 2010:
        - initializes and solves a bubble transport model that runs over InterFoam solver    

Header File
    shiPrMa.H

Author
    Lourenço Sassetti Mendes (lourenco.sassetti@gmail.com)

Date
    2018-04-16

\*---------------------------------------------------------------------------*/

#include "shiPrMa.H"

#if OPENFOAM >= 1906
    #include "gravityMeshObject.H"
#endif
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace bubbleModels
{
    defineTypeNameAndDebug(shiPrMa, 0);
    addToRunTimeSelectionTable(bubbleModel, shiPrMa, dictionary);    
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bubbleModels::shiPrMa::shiPrMa
(
    const fvMesh& mesh
)
:
    // base class
    shiBase(mesh),
    // available model components types for this bubble model
    availBubbleBins({"binBaseEntCells","binBaseEntCellsUn","binBaseEntFaces"}),
    availEntrainmentMasks({"alphaLimitFaces","alphaLimitCells"}),
    availIntergroupTransfers({"intergroupBase"}),
    // default model components types
    binType(this->lookupOrAddDefault<word>("binType", "binBaseEntCells")),
    entrainmentMaskType(this->lookupOrAddDefault<word>("entrainmentMaskType", "alphaLimitCells")),
    intergroupTransferType(this->lookupOrAddDefault<word>("intergroupTransferType", "intergroupBase")),

    // stress production
    Pr
    (
        IOobject
        (
            "Pr",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        #if OPENFOAM < 1906
            mesh_.lookupObject<volScalarField>("k").dimensions()/mesh_.lookupObject<uniformDimensionedVectorField>("g").dimensions()
        #endif
        #if OPENFOAM >= 1906
            mesh_.lookupObject<volScalarField>("k").dimensions()/meshObjects::gravity::New(mesh_.time()).dimensions()
        #endif

//        dimensionSet(0, 1, 0, 0, 0, 0, 0)
    )    
{
    InfoPos("\n--->bubbleModels::shiPrMa::shiPrMa(): start");


    // check the model components types considering the available model components types for this bubble model
    check_model_components();

    // calculate radius spacing and bubble size probability function  of each bin
    process_bins();

    // display model variables and parameters
    model_info();

    // cell volume initilization
    Vcell.ref()=mesh_.V();
//    Vcell.correctBoundaryConditions();

    // bubble bins initialization
    init_bubbleBins();

    // entrainment mask initialization
    bubbleEntMask = Foam::entrainmentMask::New(mesh, this->dict());

    // write the entrainment surface set to a file
    //entSurfSet.write();

    // update bubble populations fields (sum of all bubble bins)
    sumBins();

    // print bubble model dictionary to screen
    if (debugInfoAux) Info << "\n    bubble model properties:\n" << this->dict() << endl;
    // write the runtime bubble model dictionary 
    this->Foam::regIOobject::write();

    // write these fields for the initial time step
    Vcell.write();    
    Cb.write();
    Nb.write();
    Vb.write();

    // update flow fields used for the production of air-entrainment
    update_auxFields();


    // initialize bubble groups interactions
    if (activeBubbleInter)
    {
        // bubble intergroup transfer initialization
        bubbleInteractions = intergroupTransfer::New(mesh, this->dict(), bubbleBinsList);
    }
    else
    {
//        bubbleInteractions = NULL;
    }

    InfoPos("\n<---bubbleModels::shiPrMa::shiPrMa(): end");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bubbleModels::shiPrMa::~shiPrMa()
{
    this->Foam::regIOobject::write();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// check the model components types considering the available model components types for this bubble model
void Foam::bubbleModels::shiPrMa::check_model_components()
{
    InfoPos("\n    --->bubbleModels::shiPrMa::check_model_components() - updating the bubbleModelDict with the default model components types: start");

    if (!availBubbleBins.found(binType))
    {
        FatalErrorInFunction
            << "Incompatible or non existent bubble bin type: "
            << binType << nl << nl
            << "Valid bubble bin types:" << endl
            << availBubbleBins
            << exit(FatalError);
    }

    if (!availEntrainmentMasks.found(entrainmentMaskType))
    {
        FatalErrorInFunction
            << "Incompatible or non existent entrainment mask type: "
            << entrainmentMaskType << nl << nl
            << "Valid entrainment mask types:" << endl
            << availEntrainmentMasks
            << exit(FatalError);
    }

    if (!availIntergroupTransfers.found(intergroupTransferType))
    {
        FatalErrorInFunction
            << "Incompatible or non existent intergroup transfer type: "
            << intergroupTransferType << nl << nl
            << "Valid intergroup transfer types:" << endl
            << availIntergroupTransfers
            << exit(FatalError);
    }

    InfoPos("    <---bubbleModels::shiPrMa::check_model_components() - updating the bubbleModelDict with the default model components types: end");    
}

void Foam::bubbleModels::shiPrMa::update_auxFields()
{
    InfoPos("\n    bubbleModels::shiPrMa::update_auxFields()---->Flow fields for bubble models: initializing/updating...");

    // read gravitationalaceleration
    #if OPENFOAM < 1906
        const uniformDimensionedVectorField& g_ = mesh_.lookupObject<uniformDimensionedVectorField>("g");
    #endif
    #if OPENFOAM >= 1906
        const uniformDimensionedVectorField& g_ = meshObjects::gravity::New(mesh_.time());
    #endif
    
    // get the turbulence kinetic energy
    const volScalarField& k_ = mesh_.lookupObject<volScalarField>("k");

    // production term
    Pr = k_/(mag(g_));

    if (debug)
    {
        InfoVar(Pr);
    }


    InfoPos("    bubbleModels::shiPrMa::update_auxFields()---->Flow fields for bubble models: initialization/update completed");
}

// modified version with air bubble concentration limit at each cell
// top bounding after solving all bins Nb equation
void Foam::bubbleModels::shiPrMa::sumBins()
{
    InfoPos("\n    bubbleModels::shiPrMa::sumBins()---->bubble population sum fields: initializing/updating...");

    const dimensionedScalar volUnit("volUnits", dimVolume, 1.0);

    // volumetric fraction of bubble in cell
    // initialize with first bin
    Cb = bubbleBinsList[0].C();
    // sum all other bins
    for(int i=1; i<n_bins; i++)
    {
        // volumetric fraction of bubble in cells 
        Cb += bubbleBinsList[i].C();         
    }

    // check if bubble concentration is greater than max and needs top bounding
    if ( gMax(Cb)>Cb_max_ )
    {
            
        Info << "\n        Cb_max     =" << Cb_max_ << endl;
        Info << "        max Cb     =" << gMax(Cb) << endl;
        Info << "********    bounding bubbles    ********\n" << endl;

        // Cb bounder field: Cb<Cb_max : 1, Cb>=Cb_max: Cb_max
        tmp<volScalarField> Cb_bounder = max(Cb, Cb_max_);

        // Cb bounding

        //  bounding bubble bins fields
        for(int i=0; i<n_bins; i++)
        {
            if (debugInfoPos) Info <<  ", " <<i << flush;

            bubbleBinsList[i].Cb_bounding(Cb_bounder);
        }
    }

    // initialization
    if (debugInfoPos) Info << "\n        summing bin: " << 0 << flush;

    En = bubbleBinsList[0].En();
    Cb = bubbleBinsList[0].C();    
    Cb_evap = bubbleBinsList[0].C_evap();
    Nb = bubbleBinsList[0].Nb();
    Nb_evap = bubbleBinsList[0].Nb_evap();
    rb_avg = bubbleBinsList[0].Nb() * r_bins[0] * volUnit;

    //  sum all bubble bins fields
    for(int i=1; i<n_bins; i++)
    {
        if (debugInfoPos) Info <<  ", " <<i << flush;
        // entrainment of bubble number field
        En += bubbleBinsList[i].En();
        // volumetric fraction of bubble in cells 
        Cb += bubbleBinsList[i].C(); 
        // volumetric fraction of evaportated bubble in cells
        Cb_evap += bubbleBinsList[i].C_evap();
        // number of bubbles in each cell
        Nb += bubbleBinsList[i].Nb();
        // number of evaportaed bubbles in each cell
        Nb_evap += bubbleBinsList[i].Nb_evap();
        // average bubble radius
        rb_avg += bubbleBinsList[i].Nb() * r_bins[i] * volUnit;
    }

    // volume of bubbles
    Vb.ref() = Cb.ref() * mesh_.V();
    // volume of evaporated bubbles
    Vb_evap.ref() = Cb_evap.ref() * mesh_.V();

    // average bubble radius
    rb_avg.ref() /= Nb.ref() * volUnit + 1e-15;


    InfoPos("\n        correctBoundaryConditions() of: En, Nb, Nb_evap, Cb, Cb_evap"); //, Vb, Vb_evap, rb_avg:" << endl;

    En.correctBoundaryConditions();
    Nb.correctBoundaryConditions();
    Nb_evap.correctBoundaryConditions();
    Cb.correctBoundaryConditions();
    Cb_evap.correctBoundaryConditions();

    if (debugMode)
    {
        InfoVar(Cb);
        InfoVar(Cb_evap);
        InfoVar(En);
        InfoVar(Vb);
        InfoVar(Nb);

        if (debugInfoVar)
        {
            Info<< "\n              Sum En      : " << gSum(En) << endl;
            Info<< "              Sum Nb      : " << gSum(Nb) << endl;
            Info<< "              Sum Nb_evap : " << gSum(Nb_evap) << endl;
            Info<< "              Sum Vb      : " << gSum(Vb) << endl;
        }
    }

    InfoPos("    bubbleModels::shiPrMa::sumBins()---->bubble population sum fields: initialization/update completed");
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
