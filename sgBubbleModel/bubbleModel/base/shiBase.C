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

Header File
    shiBase.H

Author
    Lourenço Sassetti Mendes (lourenco.sassetti@gmail.com)

Date
    2018-04-16

\*---------------------------------------------------------------------------*/

#include "shiBase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bubbleModels::shiBase::shiBase
(
    const fvMesh& mesh
)
:
    // base class
    bubbleModel(mesh),
    // volume of each cell
    Vcell
    (
        IOobject
        (
            "Vcell",
            mesh.time().constant(),
            polyMesh::meshSubDir,            
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimVolume, 0.0)
    ),
    // number of bubbles field
    Nb
    (
        IOobject
        (
            "Nb",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", inv(dimVolume), 0.0)
    ),
    // number of evaportaed bubbles
    Nb_evap
    (
        IOobject
        (
            "Nb_evap",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", Nb.dimensions(), 0.0)
    ),    
    // bubble semi-"volume fraction": bubble volume divide by cell volume
    Cb
    (
        IOobject
        (
            "Cb",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    // bubble evaporated to air field semi-"volume fraction" 
    Cb_evap
    (
        IOobject
        (
            "Cb_evap",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", Cb.dimensions(), 0.0)
    ),
    // entrainment of bubble number field
    En
    (
        IOobject
        (
            "En",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", Nb.dimensions()/dimTime, 0.0)
    ),   
    // volume of bubble field
    Vb
    (
        IOobject
        (
            "Vb",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimVolume, 0.0)
    ),
    // volume of bubble evaporated to air field
    Vb_evap
    (
        IOobject
        (
            "Vb_evap",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimVolume, 0.0)
    ),
    // average bubble radius
    rb_avg
    (
        IOobject
        (
            "rb_avg",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),

    // bubbles radius list for bin i
    r_bins(this->lookup("bubbleBinsRadius")),
    // radius spacing of list bin i
    //dr_bins,
    dr_bins(r_bins*0.0),
    // bubble density N normalized by the maximum N which is the value at rb,min
    DN_bins(),
    // bubble size probability function list for bin i
    D_bins(),

//    // surface iso alpha value
//    surfaceAlpha(readScalar(this->lookup("surfaceAlpha"))),

    // air bubbles maximum volumetric fraction or entrainment and for bounding
    Cb_max_(this->lookupOrDefault<scalar>("bubbleFractionMax",1.0)),

    // number of bubble bins
    n_bins(r_bins.size())

{
    InfoPos("\n--->bubbleModels::shiBase::shiBase(): start");
    InfoPos("\n<---bubbleModels::shiBase::shiBase(): end");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bubbleModels::shiBase::~shiBase()
{
    this->Foam::regIOobject::write();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// calculate radius spacing of each bin
void Foam::bubbleModels::shiBase::calc_dr_bins()
{

    InfoPos("\n    bubbleModels::shiBase::calc_dr_bins() - Calculating the radius spacing of each bin: start");


    if (r_bins.size() > 1)
    {
        const scalar rb_min = min(r_bins);
        const scalar rb_max = max(r_bins); 

        forAll(r_bins, i)
        {
            scalar rb_i = r_bins[i];

            if ( rb_i == rb_min )
            {
                dr_bins[i] = ( r_bins[i+1] - rb_i ) / 2.0;
            }
            else if ( rb_min != rb_i && rb_i != rb_max )
            {
                dr_bins[i] = ( r_bins[i+1] - r_bins[i-1] ) / 2.0;
            }
            else if ( rb_i == rb_max )
            {
                dr_bins[i] = ( rb_i - r_bins[i-1] ) / 2.0;
            }
        }

        if (debugInfoCalc) Info<< "               sum(dr_bins)=" << sum(dr_bins) << " =? " << "rb_max-rb_min=" <<  rb_max-rb_min << endl;

        if (sum(dr_bins) != rb_max-rb_min)
        {
            FatalError
            << "    Inconsistent bubble radius spacing."
            << exit(FatalError);
        }

    }
    else if (r_bins.size() == 1)
    {
        dr_bins[0]=1.0;
    }

    if (debugInfoVar) Info<< "               dr_bins:" << endl << dr_bins << endl;

    InfoPos("    bubbleModels::shiBase::calc_dr_bins() - Calculating the radius spacing of each bin: end");
}


// calculate bubble density number normalized by the volume of all bubbles
void Foam::bubbleModels::shiBase::calc_D_bins()
{
    InfoPos("\n    bubbleModels::shiBase::calc_D_bins() - Calculating bins bubble density number normalized by the volume of all bubbles: start");

    scalar D_vol = 0.0;

    scalarList DN_vol(r_bins*0.0);

    scalar rb_min = min(r_bins);
    scalar rb_max = max(r_bins);   

    scalarList Vbubble_bins = 4.0*Foam::constant::mathematical::pi/3.0*pow3(r_bins);

    // --------------------------------------------------------------------------------
    // DN - bubble density N normalized by the maximum N which is the value at rb,min
    // according to Shi et al 2010

    // Hinze scale rH (Hinze, 1995)
    // can be calculated directed from the turbulent dissipation, surface tension and the critical Weber number
    // in this study was adopted:
    scalar rH = 1e-3;


    forAll(r_bins, i)
    {
        scalar rb_i = r_bins[i];

        if ( rb_min <= rb_i && rb_i <= rH )
        {
            DN_bins.append(
                            pow(rb_min, 3.0/2.0) * pow(rb_i, -3.0/2.0)
                          );
        }
        else if ( rH < rb_i && rb_i <= rb_max )
        {
            DN_bins.append(
                            pow(rb_min, 3.0/2.0) * pow(rH, 11.0/6.0) * pow( rb_i, -10.0/3.0)
                          );
        }

        if (debugInfoVar) Info<< "            DN_"<<  i  << " : " << DN_bins[i] << endl;
    }

    // --------------------------------------------------------------------------------
    // D_bins - bubble density number normalized by volume of all bubbles, 
    //      assuring that the sum of all DN bubble volumes is 1m3, independently of the number of bins or it's spacing


    DN_vol = DN_bins*dr_bins* Vbubble_bins;

    D_bins = DN_bins * dr_bins / sum(DN_vol);

    if (debugInfoVar) Info<< "            Discrete distibution of bubbles per each bin\n            normalized to a total volume of bubbles of 1m3: \n"    
        << "            D_bins (m-3): \n" << D_bins << endl;

    D_vol = sum(D_bins * Vbubble_bins);


    if ( abs(D_vol-1.0) > 1e-15 )
    {
        FatalError
        << "    Inconsistent normalized bins volume: D_vol"
        << "    D_vol (m3)= " << D_vol
        << "    Must be 1.0."
        << exit(FatalError);
    }

    InfoPos("    bubbleModels::shiBase::calc_D_bins()---->Calculating bins bubble density number normalized by the volume of all bubbles...end.");
}

// calculate radius spacing and bubble radius probability function  of each bin
void Foam::bubbleModels::shiBase::process_bins()
{
    // sort the r_bins list
    sort(r_bins);

    // calculate the radius spacing of each bin
    calc_dr_bins();

    // calculate the bins prbability function
    calc_D_bins();

    // updating the bubbleModelDict
    this->add("bubbleBinsRadiusSpacing", dr_bins);
    this->add("bubbleBinsDensity", D_bins);
}

// display model variables and parameters
void Foam::bubbleModels::shiBase::model_info()
{
    InfoPos("\n    bubbleModels::shiBase::model_info()---->display model variables and parameters: start");
    
    // display model status
    model_status();

    Info<< "\n        - Intergroup coalescence and breakup activated: " << activeBubbleInter << endl;

    // Print averaged bubble bin variables values
    Info << "\n        - Number of bubble bins     : " << n_bins << endl;

    Info << "\n            bins radius (m)               : " << r_bins << endl;
    Info << "            bins radius spacing (m)         : " << dr_bins << endl;
    Info << "            bins density vol Normalized PDF : " << D_bins << endl;

    Info << "\n        - Average bubble bin variables:" << endl;
    Info << "            - avg radius (m)          : " << average(r_bins) << endl;
    Info << "            - avg radius spacing (m)  : " << average(dr_bins) << "\n" << endl;
   
    Info<< "\n        - model type                    : " << this->lookup("modelType") << endl;
    Info<< "        - prod. term constant 'ab'      : " << this->lookup("ab") << endl; 
    Info<< "        - prod. onset threshold 'Pr0'   : " << this->lookup("productionThreshold") << endl;
//    Info<< "        - prod. onset density source    : " << Pr_rhoName << endl;
    Info<< "        - bubbles Schmidt Number        : " << this->lookup("bubbleSchmidtNumber") << endl;

    Info<< "\n        - surface alpha for bubble entrainment:" << this->lookup("surfaceAlpha") << endl;   
    Info<< "\n        - surface alpha threshold for the evaporation:" << this->lookup("EnAlphaMin") << endl;

    Info<< "\n******************************************************************************************\n\n" << endl;


    InfoPos("    bubbleModels::shiBase::model_info()---->display model variables and parameters: end");

}

void Foam::bubbleModels::shiBase::init_bubbleBins()
{

    InfoPos("\n    bubbleModels::shiBase::init_bubbleBins()---->Dynamic bubbleBins initialization: starting...");

    //  dynamically allocating bubbleBins in the list by the "new" operator
    for(int i=0; i<r_bins.size(); i++)
    {
        // creating the bubble bin suffix: "bin_##"
        char binIndex [5];
        sprintf(binIndex, "%02i", i+1);
        
        std::ostringstream oss;
        oss << "bin_" << binIndex;
        const word bubbleBinName = oss.str();

        if (debugInfoPos) Info << "\n         ---->" << bubbleBinName << endl;

        // create a new dict for each bin with some keys of the bubblePropertiesDict
        IOdictionary bubbleBinDict
        (
            IOobject
            (
                "bubbleBinDict_" + bubbleBinName,
//                mesh_time().constant(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            )
        );
        bubbleBinDict += this;

        bubbleBinDict.remove("bubbleModelActive");
        bubbleBinDict.remove("bubbleBinsRadius");
        bubbleBinDict.remove("bubbleBinsRadiusSpacing");
        bubbleBinDict.remove("bubbleBinsDensity");

        bubbleBinDict.add("binName", bubbleBinName);
        bubbleBinDict.add("bubbleRadius", r_bins[i]);
        bubbleBinDict.add("radiusSpacing", dr_bins[i]);
        bubbleBinDict.add("bubbleDensity", D_bins[i]);

        if (debugInfoAux) bubbleBinDict.Foam::regIOobject::write();

        // append model name to bubble models name list 
        bubbleBinsNameList.append(bubbleBinName);

        //append bubble bin to bubble bins pointer list 
        bubbleBinsList.append
        (
            // bin initialization
            autoPtr<bin> 
            (
                bin::New(mesh_, bubbleBinDict)
            )
        );
        if (debugInfoPos) Info << "\n         appended---->" << bubbleBinName << endl;
    }

    InfoPos("    bubbleModels::shiBase::init_bubbleBins()---->Dynamic bubbleBins initialization: completed");
}



// intergroup bubble source/sink term for the bubbles of air transport equation
void Foam::bubbleModels::shiBase::update_IG()
{ 
    InfoPos("\n    bubbleModels::shiBase::update_IG(): start");

    if (activeBubbleInter)
    {
        bubbleInteractions->update();
    }
    else
    {
        InfoAux("    bubbleModels::shiBase::update_IG()----> no intergroup transfer");
    }

    InfoPos("    bubbleModels::shiBase::update_IG(): end");    
}

void Foam::bubbleModels::shiBase::solve()
{
    InfoPos("\n    bubbleModels::shiBase::solve()---->All bubble populations: updating & solving...");

    // ~~~ identify entrainment surface faces/cells mask
    // ~~~ with/without optional alpha gradient unitary vectors
    InfoPos("\n    bubbleModels::shiBase::solve()---->calling 'bubbleEntMask->update()'");
    bubbleEntMask->update();


    InfoPos("\n    bubbleModels::shiBase::solve()---->calling 'update_IG()'");
    // intergroup bubble source/sink term for the bubbles of air transport equation
    update_IG();

    // ~~~ check if the sum of interchange bubble bins volume is near to zero
    if (debugMode)
    {
        InfoPos("\n            bubbleModel::solve()---->All bubble populations: inter-group change bubble volume: ");
        scalar IG_Vb_sum = 0.0;

        forAll(bubbleBinsNameList, i)
        {
            if (debugInfoAux) Info<< i << ", " << flush;

            IG_Vb_sum += gSum(bubbleBinsList[i].IG()) * bubbleBinsList[i].VOneBubble().value();
        }

        if (debugInfoCalc) Info << "\n    IG sum(volumes) = " << IG_Vb_sum << endl;
        if (debugInfoCalc) Info << "    Vb sum(volumes) = " << gSum(Vb) << endl << endl;
    }

    // ~~~ solve bubble populations models information
    forAll(bubbleBinsNameList, i)
    {
        if (debugInfoPos) Info<<"\n    "<< bubbleBinsList[i].name() << endl;

        bubbleBinsList[i].update_Ub();
        bubbleBinsList[i].entrainment();
        bubbleBinsList[i].solve_Nb();
    }

    InfoPos("    bubbleModels::shiBase::solve()---->All bubble populations: updating & solving:end");
}


void Foam::bubbleModels::shiBase::reCalc()
{
    InfoPos("\n--->bubbleModels::shiBase::recalc()---->runTime bubbles model calculations: start...");

    // update flow fields for production of air-entrainment
    update_auxFields();
    
    // solve all bubble models
    solve();
    
    // update bubble populations fields (sum of all bubble bins)
    sumBins();

    InfoPos("\n<---bubbleModels::shiBase::recalc()---->runTime bubbles model calculations: end...\n");
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
