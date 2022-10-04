/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::bubbleModel> Foam::bubbleModel::New
(
    const fvMesh& mesh
)
{
//    Info << nl << "--->bubbleModel::New(mesh): start" << endl;

    // read the bubble model properties dictionary
    IOdictionary inputDict
    (
        IOobject
        (
            propertiesName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // check if the bubble model is active or not
    const bool activeModel_(inputDict.lookupOrDefault<Switch>("bubbleModelActive", false));

//    Info << nl << "    bubbleModel is ACTIVE?: " << activeModel_ << endl;

    if (activeModel_)
    {
        // get the bubble model type
        const word modelType_(inputDict.lookup("modelType"));

//        Info<< "    bubbleModel::New(mesh): Selecting bubble model: " << modelType_ << endl;

        auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType_);

        if (!cstrIter.found())
        {
            FatalErrorInFunction
                << "Unknown BubbleModel type: "
                << modelType_ << nl << nl
                << "Valid BubbleModel types:" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }
/*        // verbose
        {
            Info<< "    bubbleModel::New(dictionary): Selected bubble model: " << modelType_ << endl;
        } 
*/
        return autoPtr<bubbleModel>(cstrIter()(mesh));
    }
    else
    {
        return autoPtr<bubbleModel>(new bubbleModel(mesh));
    }

//    Info << "--->bubbleModel::New(bubbleModelDict): end" << endl;        
}



// ************************************************************************* //
