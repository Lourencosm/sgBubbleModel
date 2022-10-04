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

#include "entrainmentMask.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


namespace Foam
{

    defineTypeNameAndDebug(entrainmentMask, 0);
    defineRunTimeSelectionTable(entrainmentMask, dictionary);    
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::entrainmentMask::entrainmentMask
(
    const fvMesh& mesh,
    const dictionary& bubbleModelDict
)
:
    // fvMesh
    mesh_(mesh),
    // bubble population dictionary
    bubbleModelDict_(bubbleModelDict),
    // bubble population name
    entMaskType_(bubbleModelDict.lookup("entrainmentMaskType")),
    // debug mode
    debugMode(bubbleModelDict.lookupOrDefault<Switch>("debugMode", false)),
    // debug messages
    debugInfoPos(bubbleModelDict.lookupOrDefault<Switch>("debugInfoPositional", false)),
    debugInfoAux(bubbleModelDict.lookupOrDefault<Switch>("debugInfoAuxiliary", false)),
    debugInfoVar(bubbleModelDict.lookupOrDefault<Switch>("debugInfoVariables", true)),
    debugInfoCalc(bubbleModelDict.lookupOrDefault<Switch>("debugInfoCalculations", false))    
{
    InfoPos("--->entrainmentMask::entrainmentMask(mesh): start");

    status();

    InfoPos("<---entrainmentMask::entrainmentMask(mesh): end");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::entrainmentMask::~entrainmentMask()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// display model name, activation status and debug mode
void Foam::entrainmentMask::status()
{
    InfoPos("\n    entrainmentMask::status()---->display entrainment masks variables and parameters: start");
    

        Info<< "\n\n******************************************************************************************" << endl;
        Info<< "    Bubble entrainment mask" << endl;
        Info<< "                       type: " << entMaskType_ << endl;
        Info<< "******************************************************************************************\n" << endl;


    InfoPos("    entrainmentMask::status()---->display entrainment masks variables and parameters: end");
}

void Foam::entrainmentMask::InfoPos(const word& msg)
{
    if (debugInfoPos) Info<< msg << endl;
}

void Foam::entrainmentMask::InfoAux(const word& msg)
{
    if (debugInfoAux) Info<< msg << endl;
}

void Foam::entrainmentMask::InfoVar(const volScalarField& sf)
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
void Foam::entrainmentMask::InfoVar(const surfaceScalarField& sf)
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
