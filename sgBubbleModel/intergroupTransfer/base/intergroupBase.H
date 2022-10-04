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

Class
    Foam::intergroupTransfers::intergroupBase

Description
    Base class for intergroup bubble transfers by coalescence and breakup
    according to Shi et al. (2010)

SourceFiles
    intergroupBase.C

Author
    Lourenço Sassetti Mendes (lourenco.sassetti@gmail.com)

\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
bubble Group Interaction class for air bubble entrainment and transport model 
belongs to interFoam_sgbm solver version 64 and forward

\*---------------------------------------------------------------------------*/
#ifndef intergroupBase_H
#define intergroupBase_H

//#include "bubbleModel.H"
#include "bin.H"
#include "intergroupTransfer.H"

// openfoam libraries
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace intergroupTransfers
{

/*---------------------------------------------------------------------------*\
                           Class bubbleInter Declaration
\*---------------------------------------------------------------------------*/


class intergroupBase
:
    public intergroupTransfer
{
// Private data
//private:

// Protected data
protected:
        scalarList  r_bins_;
        int  n_bins_;
        labelList  mi;
        labelList  mj;

        dimensionedScalar timeUnit;

        const dictionary& transportPropertiesDict_;

        const dimensionedScalar rho_ref_;
        const dimensionedScalar rho_water_;
        const dimensionedScalar rho_air_;
        const dimensionedScalar nu_water_;
        const dimensionedScalar nu_air_;        
        const dimensionedScalar sigma_;

        const scalar beta_;

        dimensionSet dimIG;
        dimensionedScalar IGUnit;

        scalar***   Xikl_;
        scalar**    BrXik_;
        

        PtrList<dimensionedScalar> r_b_;

        PtrList<volScalarField> IG_;

        PtrList<volScalarField> Xplus_;
        PtrList<volScalarField> Xminus_;
        PtrList<volScalarField> Bplus_;
        PtrList<volScalarField> Bminus_;

        PtrList<volScalarField> T_;
        PtrList<volScalarField> W_;
        PtrList<volScalarField> CE_;
        PtrList<volScalarField> BrKernal_;

public:

    //- Runtime type information
    TypeName("intergroupBase")

    // Constructors

        //- Construct from components
        intergroupBase
        (
            const fvMesh& mesh,
            const dictionary& bubbleBinDict,
            const PtrList<bin>& bubbleBinsList
        );

    // Destructor
        virtual ~intergroupBase();

//        //- Return const-access to source term for bubble intergroup transfer in the number of bubble eq.
        const PtrList<volScalarField>& IG() const
        {
            return IG_;
        }

        // update the source term for bubble intergroup transfer in the number of bubble eq.
//        void update_Sn(const int& i, const bubbleInter& bubbleInteractions_);

        // update all fields
        virtual void update();

//protected:
    // index function
    // converts two indexes i,j into a single index
    // i,j should return the same as j,i
    int idx( int i, int j )
    {
        int m=0;

        if (i < j)
        {
            m = i * n_bins_ + j - i*(i+1)/2;
        }
        else
        {
            m = j * n_bins_ + i - j*(j+1)/2;            
        }

        return m;
    }

    // initialize the correspondence between ij and m indexes
    void ij_init( )
    {
        int m = 0;

        InfoAux("\n\n            i,j  indexes in function of m\n");

        forAll(r_bins_,i)
        {    
            forAll(r_bins_,j)
            {
                // create m as function funtion of m
                if (i <= j)
                {
                    mi[m] = i;
                    mj[m] = j;
                    if (debugInfoAux) Info<< "            ij[" << m << "]= (" << mi[m] << "," << mj[m] << ")" << endl;
                    m++;
                }
            }
        }

        Info<< endl << endl << endl;

    }

    // initiate BrKernal_ , BrXik_, Xikl_ scalar arrays
    virtual void init_arrays();

    // calculate the collision rate of bubble groups
    virtual void init_field(PtrList<volScalarField>& F_, const word& fieldName_, const dimensionSet& dimset_);

    // calculate the collision rate of bubble groups
    virtual void calc_T();

    // calculate the Weber number
    virtual void calc_W();

    // calculate the number of bubble transfered from the calescence of two bubble from group k and l to group i
    virtual void calc_Xikl();

    // calculate the coalescence efficiency
    virtual void calc_CE();

    // calculate the coalencence source
    virtual void calc_Xplus();

    // calculate the coalencence sink
    virtual void calc_Xminus();

    // calculate the bubble breakup source
    virtual void calc_Bplus();

    // calculate the bubble breakup sink
    virtual void calc_Bminus();

    // calcutate the breakup kernal function
    virtual void calc_BrKernal();

    // indefinite integral of the integral in the breakage kernal expression of Luo and Svendsen (1996)
    Foam::scalar ksi_indef_integral(const scalar ksi_, const scalar a_);

    // function inside the integral in the breakage kernal expression of Luo and Svendsen (1996)
    Foam::scalar ksi_func(const scalar ksi_, const scalar a_);

    // function inside the integral in the breakage kernal expression of Luo and Svendsen (1996)
    Foam::scalar ksi_func_RiemannIntegral(const scalar a_, const scalar xmin_, const scalar xmax_, const int N_);

    // calculate the number of bubble transfered from the breakup of bubble into two identical daughter bubbles
    virtual void calc_BrXik();

    // calculate the source term for bubble intergroup transfer in the number of bubble eq.
    virtual void calc_IG();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace intergroupTransfers 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
