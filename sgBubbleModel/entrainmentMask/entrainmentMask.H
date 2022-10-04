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
    Foam::entrainmentMask

Description
    Base class for entrainment masks

SourceFiles
    entrainmentMask.C    
    entrainmentMaskNew.C

Author
    Lourenço Sassetti Mendes (lourenco.sassetti@gmail.com)


\*---------------------------------------------------------------------------*/

#ifndef entrainmentMask_H
#define entrainmentMask_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "IOdictionary.H"

#include "fvMesh.H"

// header for InfoVar()
#include "fvc.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class entrainmentMask Declaration
\*---------------------------------------------------------------------------*/

class entrainmentMask
{

protected:
        const fvMesh& mesh_;
        const dictionary& bubbleModelDict_;

        word entMaskType_;

        bool debugMode;
        bool debugInfoPos;
        bool debugInfoAux;
        bool debugInfoVar;        
        bool debugInfoCalc; 
    // Protected Member Functions

private:    

    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        entrainmentMask(const entrainmentMask&);

        //- Disallow default bitwise assignment
        void operator=(const entrainmentMask&);


public:

    //- Runtime type information
    TypeName("entrainmentMask");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        entrainmentMask,
        dictionary,
        (
            const fvMesh& mesh_,
            const dictionary& bubbleModelDict_
        ),
        (mesh_, bubbleModelDict_)
    );

    // Constructors

    //- Construct from components
    entrainmentMask
    (
        const fvMesh& mesh_,
        const dictionary& bubbleModelDict_
    );

    //- Destructor
    virtual ~entrainmentMask();

    // Selectors

        //- Return a reference to the selected bubble model
        static autoPtr<entrainmentMask> New
        (
            const fvMesh& mesh_,
            const dictionary& bubbleModelDict_
        );

    // Public Member Functions

        //- Return const-access to entrainmentMask type
        virtual const word& entMasktype() const
        {
            return entMaskType_;
        }

        // display model name, activation status and debug mode
        virtual void status();

        // update the entrainmentMask
        virtual void update() = 0;

        // Info messages - code position
        void InfoPos(const word& msg);

        // Info messages - auxiliary info
        void InfoAux(const word& msg);

        // Info messages - variable info: min,max,avg, dimensions
        void InfoVar(const volScalarField& sf);
        void InfoVar(const surfaceScalarField& sf);        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
