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

#include "intergroupTransfer.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::intergroupTransfer> Foam::intergroupTransfer::New
(
    const fvMesh& mesh_,
    const dictionary& bubbleModelDict_,
    const PtrList<bin>& bubbleBinsList_
)
{
    // get the intergroup transfers type
    const word intergroupType_(bubbleModelDict_.lookup("intergroupTransferType"));

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(intergroupType_);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown intergroup transfers type: "
            << intergroupType_ << nl << nl
            << "Valid entrainment intergroup transfers types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    return autoPtr<intergroupTransfer>(cstrIter()(mesh_, bubbleModelDict_, bubbleBinsList_));
}



// ************************************************************************* //
