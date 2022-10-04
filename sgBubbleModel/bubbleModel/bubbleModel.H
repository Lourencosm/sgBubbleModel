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
    Foam::bubbleModel

Description
    Base class for bubbleModels

SourceFiles
    bubbleModel.C    
    bubbleModelNew.C

Author
    Lourenço Sassetti Mendes (lourenco.sassetti@gmail.com)


\*---------------------------------------------------------------------------*/

#ifndef bubbleModel_H
#define bubbleModel_H

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
                          Class bubbleModel Declaration
\*---------------------------------------------------------------------------*/

class bubbleModel
:
    public IOdictionary
{

protected:

    // Protected data
        const fvMesh& mesh_;
//        const fileName inputDictFile;
//        const dictionary& inputDict;
//        const IOdictionary& modelDict;
        bool activeModel;
        word modelType;
        bool activeBubbleInter;        
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
        bubbleModel(const bubbleModel&);

        //- Disallow default bitwise assignment
        void operator=(const bubbleModel&);


public:

    //- Runtime type information
    TypeName("bubbleModel");

    //- Default name of the bubble properties dictionary
    static const word propertiesName;

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        bubbleModel,
        dictionary,
        (
            const fvMesh& mesh_
        ),
        (mesh_)
    );

    // Constructors

    //- Construct from components
    bubbleModel
    (
        const fvMesh& mesh_
    );

    //- Destructor
    virtual ~bubbleModel();

    // Selectors

        //- Return a reference to the selected bubble model
        static autoPtr<bubbleModel> New
        (
            const fvMesh& mesh_
        );

    // Public Member Functions

        //- Return the bubbleModel dictionary
//        inline const dictionary& dict() const;

        virtual const IOdictionary& dict() const
        {
//            const IOdictionary& dict_(this);
            return mesh_.lookupObject<IOdictionary>(propertiesName + "_runTime");
//            return dict_;
        }

        // check bubble model activation status
        virtual const bool& isActive() const
        {
            return activeModel;
        }

        //- Return const-access to bubbleBin model type
        virtual const word& modeltype() const
        {
            return modelType;
        }

        // display model name, activation status and debug mode
        virtual void model_status();

        // runtime cycle tasks. e.g: update_auxFields() + solve() + sumBins()
        virtual void reCalc();

        // Info messages - code position
        void InfoPos(const word& msg);

        // Info messages - auxiliary info
        void InfoAux(const word& msg);

        // Info messages - variable info: min,max,avg, dimensions
//        template<class T>
//        template<typename T>
//        void InfoVar(const T& sf);
        void InfoVar(const volScalarField& sf);
        void InfoVar(const surfaceScalarField& sf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
