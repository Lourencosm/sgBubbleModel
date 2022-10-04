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

#include "bin.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


namespace Foam
{

    defineTypeNameAndDebug(bin, 0);
    defineRunTimeSelectionTable(bin, dictionary);    
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bin::bin
(
    const fvMesh& mesh,
    const dictionary& bubbleBinDict    
)
:
    // fvMesh
    mesh_(mesh),
    // bin dictionary
    bubbleBinDict_(bubbleBinDict),
    // bubble population name
    binType_(bubbleBinDict.lookup("binType")),
    // debug mode
    debugMode(bubbleBinDict.lookupOrDefault<Switch>("debugMode", false)),
    // debug messages
    debugInfoPos(bubbleBinDict.lookupOrDefault<Switch>("debugInfoPositional", false)),
    debugInfoAux(bubbleBinDict.lookupOrDefault<Switch>("debugInfoAuxiliary", false)),
    debugInfoVar(bubbleBinDict.lookupOrDefault<Switch>("debugInfoVariables", true)),
    debugInfoCalc(bubbleBinDict.lookupOrDefault<Switch>("debugInfoCalculations", false))    
{
    InfoPos("--->bin::bin(mesh): start");

    bin_status();

    InfoPos("<---bin::bin(mesh): end");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bin::~bin()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// display model name, activation status and debug mode
void Foam::bin::bin_status()
{
    InfoPos("\n    bin::bin_status()---->display bin variables and parameters: start");
    

        Info<< "\n\n******************************************************************************************" << endl;
        Info<< "    Bubble bin" << endl;
        Info<< "          type: " << binType_ << endl;
        Info<< "******************************************************************************************\n\n" << endl;


    InfoPos("    bin::bin_status()---->display bin variables and parameters: end");
}

void Foam::bin::InfoPos(const word& msg)
{
    if (debugInfoPos) Info<< msg << endl;
}

void Foam::bin::InfoAux(const word& msg)
{
    if (debugInfoAux) Info<< msg << endl;
}

void Foam::bin::InfoVar(const volScalarField& sf)
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
void Foam::bin::InfoVar(const surfaceScalarField& sf)
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
