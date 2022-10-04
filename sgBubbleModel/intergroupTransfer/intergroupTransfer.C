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

#include "intergroupTransfer.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


namespace Foam
{

    defineTypeNameAndDebug(intergroupTransfer, 0);
    defineRunTimeSelectionTable(intergroupTransfer, dictionary);    
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::intergroupTransfer::intergroupTransfer
(
    const fvMesh& mesh,
    const dictionary& bubbleModelDict,
    const PtrList<bin>& bubbleBinsList
)
:
    // fvMesh
    mesh_(mesh),
    // bubble model dictionary
    bubbleModelDict_(bubbleModelDict),
    // bubble bins list
    bubbleBinsList_(bubbleBinsList),
    // intergroup transfer type
    intergroupType_(bubbleModelDict.lookup("intergroupTransferType")),
    // debug mode
    debugMode(bubbleModelDict.lookupOrDefault<Switch>("debugMode", false)),
    // debug messages
    debugInfoPos(bubbleModelDict.lookupOrDefault<Switch>("debugInfoPositional", false)),
    debugInfoAux(bubbleModelDict.lookupOrDefault<Switch>("debugInfoAuxiliary", false)),
    debugInfoVar(bubbleModelDict.lookupOrDefault<Switch>("debugInfoVariables", true)),
    debugInfoCalc(bubbleModelDict.lookupOrDefault<Switch>("debugInfoCalculations", false))    
{
    InfoPos("--->intergroupTransfer::intergroupTransfer(mesh, bubbleModelDict): start");

    intergroupTransfer_status();

    InfoPos("<---intergroupTransfer::intergroupTransfer(mesh, bubbleModelDict): end");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::intergroupTransfer::~intergroupTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// display model name, activation status and debug mode
void Foam::intergroupTransfer::intergroupTransfer_status()
{
    InfoPos("\n    intergroupTransfer::intergroupTransfer_status()---->display intergroup transfer variables and parameters: start");
    

        Info<< "\n\n******************************************************************************************" << endl;
        Info<< "    Bubble intergroup transfer" << endl;
        Info<< "                          type: " << intergroupType_ << endl;
        Info<< "******************************************************************************************\n\n" << endl;


    InfoPos("    intergroupTransfer::intergroupTransfer_status()---->display intergroup transfer variables and parameters: end");
}

void Foam::intergroupTransfer::InfoPos(const word& msg)
{
    if (debugInfoPos) Info<< msg << endl;
}

void Foam::intergroupTransfer::InfoAux(const word& msg)
{
    if (debugInfoAux) Info<< msg << endl;
}

void Foam::intergroupTransfer::InfoVar(const volScalarField& sf)
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
void Foam::intergroupTransfer::InfoVar(const surfaceScalarField& sf)
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
