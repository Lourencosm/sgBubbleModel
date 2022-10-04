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

#include "bin.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::bin> Foam::bin::New
(
    const fvMesh& mesh_,
    const dictionary& bubbleBinDict_
)
{

    // get the bin type
    const word binType_(bubbleBinDict_.lookup("binType"));

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(binType_);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown bubble bin type: "
            << binType_ << nl << nl
            << "Valid bubble bin types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<bin>(cstrIter()(mesh_, bubbleBinDict_));
    
}



// ************************************************************************* //
