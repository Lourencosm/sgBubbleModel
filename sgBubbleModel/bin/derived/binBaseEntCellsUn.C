/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "binBaseEntCellsUn.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace bins
{
    defineTypeNameAndDebug(binBaseEntCellsUn, 0);
    addToRunTimeSelectionTable(bin, binBaseEntCellsUn, dictionary);    
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bins::binBaseEntCellsUn::binBaseEntCellsUn
(
    const fvMesh& mesh,
    const dictionary& bubbleBinDict   
)
:
    // base class
    binBaseEntFaces(mesh, bubbleBinDict)

{
    InfoPos("\n        bins::binBaseEntCellsUn::binBaseEntCellsUn()---->initializated ");
}    


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bins::binBaseEntCellsUn::~binBaseEntCellsUn()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// calculate air entranment in cells - method B - cells
void Foam::bins::binBaseEntCellsUn::entrainment()
{

    const volScalarField& entrainmentMask_cells = mesh_.lookupObject<volScalarField>("entrainmentMask_cells_scalar");
    const volScalarField& Pr_ = mesh_.lookupObject<volScalarField>("Pr");

    const volVectorField& U_ = mesh_.lookupObject<volVectorField>("U");
    const volScalarField& alpha_ = mesh_.lookupObject<volScalarField>("alpha.water");

    // gradient of alpha @cells
    tmp<volVectorField> alpha_grad = fvc::grad(alpha_);

    // air entranment flux sum of all cell faces @ cell center
    En_ = entrainmentMask_cells * pos(Pr_ - Pr0_threshold_) * pos( U_ & alpha_grad ) * Pr_ * a_b_ * Di_b_ ;
    // En_.dimensions=[0 -3 -1 ...]     

    // entrainment, En_unlimited
    if (debugInfoCalc)
    {
        Info<< "            En_ (unlimited) :" << endl; 
        InfoVar(En_);
    }
    
    // bound the entrainment bubble volume to only Cb_max*VCell
    // by definition the sum(Di_b*Vbubble_i)=1.0 m3
    // so: En <= Di_b / m3 x Cb_max * dTime
    // so: En / Di_b <= Cb_max * dtime

    // _perUnitVol per time
    volScalarField Nb_remaining_capacity = max( (Cb_max_ * Di_b_ - Nb_), Nb_*0.0 );


    if (debugInfoCalc) Info<< "            Nb_remaining_capacity                 :            " << gMin(Nb_remaining_capacity) << " " << gMax(Nb_remaining_capacity) << endl;

    En_ = min( En_, Nb_remaining_capacity / dimensionedScalar("dTime", dimTime, mesh_.time().deltaTValue()) );


    if (debugInfoCalc)
    {
        Info<< "\n        min, max, avg of: " << name_ << ":" << endl;
        Info<< "            prod @ cells, Pr:            " << gMin(Pr_) << " " << gMax(Pr_) << " " << gAverage(Pr_) << endl;
        Info<< "            entrainment, En_       :\n            " << gMin(En_) << " " << gMax(En_) << " " << gAverage(En_) << endl;

        InfoVar(En_);
    }

}
// ************************************************************************* //
