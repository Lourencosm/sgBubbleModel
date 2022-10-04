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

\*---------------------------------------------------------------------------*/

#include "bubbleModel.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


namespace Foam
{
    defineTypeNameAndDebug(bubbleModel, 0);
    defineRunTimeSelectionTable(bubbleModel, dictionary);    
}

const Foam::word Foam::bubbleModel::propertiesName("bubbleProperties");


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bubbleModel::bubbleModel
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            propertiesName + "_runTime",
            mesh.time().constant(),           
            mesh,
            IOobject::NO_READ,
//            IOobject::AUTO_WRITE // <--- replaces the original dict every timestep
            IOobject::NO_WRITE // <--- last working
        ),
        mesh.lookupObject<IOdictionary>(propertiesName)
    ),
    // fvMesh
    mesh_(mesh),
    // bubble model activation status
    activeModel(this->lookupOrDefault<Switch>("bubbleModelActive", false)),
    // bubble model type - author/reference/name
    modelType(this->lookup("modelType")),
    // Is the bubbles intergroup model active?
    activeBubbleInter(readBool(this->lookup("bubbleIntergroupActive"))),    
    // debug mode
    debugMode(this->lookupOrDefault<Switch>("debugMode", false)),
    // debug messages
    debugInfoPos(this->lookupOrDefault<Switch>("debugInfoPositional", false)),
    debugInfoAux(this->lookupOrDefault<Switch>("debugInfoAuxiliary", false)),
    debugInfoVar(this->lookupOrDefault<Switch>("debugInfoVariables", true)),
    debugInfoCalc(this->lookupOrDefault<Switch>("debugInfoCalculations", false))    
{
    InfoPos("--->bubbleModel:bubbleModel(mesh): start");


    InfoPos("<---bubbleModel:bubbleModel(mesh): end");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bubbleModel::~bubbleModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// display model name, activation status and debug mode
void Foam::bubbleModel::model_status()
{
    InfoPos("\n    bubbleModel::model_info()---->display model variables and parameters: start");
    
    if (activeModel)
    {
        Info<< "******************************************************************************************" << endl;
        Info<< "    Bubble model is ACTIVE" << endl;
        Info<< "              model type: " << modelType << endl;
        Info<< "              debug mode is: " << debugMode << endl; 
        Info<< "    Running interFoam with sub-grid bubble model" << endl;
        Info<< "******************************************************************************************" << endl;
    }
    else
    {
        Info<< "******************************************************************************************" << endl;
        Info<< "    Bubble model is NOT active" << endl;
        Info<< "    Running interFoam" << endl;
        Info<< "******************************************************************************************" << endl;          
    }

    InfoPos("    bubbleModel::model_info()---->display model variables and parameters: end");
}

void Foam::bubbleModel::reCalc()
{
    FatalErrorInFunction
        << "\n\n\n****--->bubbleModel::recalc(): "
        << "***empy function**** should not happen!! " 
        << "exiting right now! " << endl
        << exit(FatalError);    
}

void Foam::bubbleModel::InfoPos(const word& msg)
{
    if (debugInfoPos) Info<< msg << endl;
}

void Foam::bubbleModel::InfoAux(const word& msg)
{
    if (debugInfoAux) Info<< msg << endl;
}

//void Foam::bubbleModel::InfoVar(const volScalarField& sf)
//template<class T>
//template<typename T>
//void Foam::bubbleModel::InfoVar(const T& sf)
void Foam::bubbleModel::InfoVar(const volScalarField& sf)
{
    if (debugInfoVar)
    {
        Info<< "        "
        << sf.name()
        << "       min, max, avg, dims:  "
        << gMin(sf) << ",  "
        << gMax(sf) << ",  " 
        << gAverage(sf) << "  " 
        << sf.dimensions()  << endl;
    }
}
void Foam::bubbleModel::InfoVar(const surfaceScalarField& sf)
{
    if (debugInfoVar)
    {
        Info<< "        "
        << sf.name()
        << "       min, max, avg, dims:  "
        << gMin(sf) << ",  "
        << gMax(sf) << ",  " 
        << gAverage(sf) << "  " 
        << sf.dimensions()  << endl;
    }
}
// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
