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
    Foam::intergroupTransfer

Description
    Base class for intergroup bubble transfers by coalescence and breakup

SourceFiles
    intergroupTransfer.C    
    intergroupTransferNew.C

Author
    Lourenço Sassetti Mendes (lourenco.sassetti@gmail.com)


\*---------------------------------------------------------------------------*/

#ifndef intergroupTransfer_H
#define intergroupTransfer_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fvMesh.H"
#include "bin.H"
#include "fvc.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class intergroupTransfer Declaration
\*---------------------------------------------------------------------------*/

class intergroupTransfer
{

protected:
        const fvMesh& mesh_;
        const dictionary& bubbleModelDict_;
        const PtrList<bin>& bubbleBinsList_;

        word intergroupType_;

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
        intergroupTransfer(const intergroupTransfer&);

        //- Disallow default bitwise assignment
        void operator=(const intergroupTransfer&);


public:

    //- Runtime type information
    TypeName("intergroupTransfer");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        intergroupTransfer,
        dictionary,
        (
            const fvMesh& mesh_,
            const dictionary& bubbleModelDict_,
            const PtrList<bin>& bubbleBinsList_
        ),
        (mesh_, bubbleModelDict_, bubbleBinsList_)
    );

    // Constructors

    //- Construct from components
    intergroupTransfer
    (
        const fvMesh& mesh_,
        const dictionary& bubbleModelDict_,
        const PtrList<bin>& bubbleBinsList_       
    );

    //- Destructor
    virtual ~intergroupTransfer();

    // Selectors

        //- Return a reference to the selected bubble model
        static autoPtr<intergroupTransfer> New
        (
            const fvMesh& mesh_,
            const dictionary& bubbleModelDict_,
            const PtrList<bin>& bubbleBinsList_
        );

    // Public Member Functions

        //- Return const-access to intergroupTransfer type
        virtual const word& intergroupType() const
        {
            return intergroupType_;
        }

        // display model name, activation status and debug mode
        virtual void intergroupTransfer_status();

        // update the intergroupTransfer
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
