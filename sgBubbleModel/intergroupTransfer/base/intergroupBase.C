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

#include "intergroupBase.H"

#include <boost/config.hpp>
#include <boost/math/special_functions/gamma.hpp>

using namespace boost;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace intergroupTransfers
{
    defineTypeNameAndDebug(intergroupBase, 0);
    addToRunTimeSelectionTable(intergroupTransfer, intergroupBase, dictionary);    
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::intergroupTransfers::intergroupBase::intergroupBase
(
    const fvMesh& mesh,
    const dictionary& bubbleModelDict,
    const PtrList<bin>& bubbleBinsList
)
:
    // base class
    intergroupTransfer(mesh, bubbleModelDict, bubbleBinsList),
    // bubbles radius list for bin i
    r_bins_(bubbleModelDict_.lookup("bubbleBinsRadius")),
    // number of bubble bins
    n_bins_(r_bins_.size()),
    mi(n_bins_*(n_bins_+1)/2),
    mj(n_bins_*(n_bins_+1)/2),

    // unitary dimensioned scalar of time
    timeUnit("timeUnit", dimTime, 1.0),

    // transport propertis mesh dict
    transportPropertiesDict_(mesh.lookupObject<dictionary>("transportProperties")),

    rho_ref_("rho_ref", dimDensity, transportPropertiesDict_.subDict("water").lookup("rho")),
    rho_water_("rho_water", dimDensity, transportPropertiesDict_.subDict("water").lookup("rho")),
    rho_air_("rho_air", dimDensity, transportPropertiesDict_.subDict("air").lookup("rho")),
    nu_water_("nu_water", dimViscosity, transportPropertiesDict_.subDict("water").lookup("nu")),
    nu_air_("nu_air", dimViscosity, transportPropertiesDict_.subDict("water").lookup("nu")),
    sigma_("sigma", dimForce/dimLength, transportPropertiesDict_.lookup("sigma")),

    beta_(2.0), // Kuboi et al. (1972a)=2.0 or Luo and Svendsen (1996)=2.046

    dimIG(bubbleBinsList_[0].IG().dimensions()),
    // unitary dimensioned scalar of IG    
    IGUnit("IGUnit", dimIG, 1.0)
{

    InfoPos("\n--->intergroupTransfers::intergroupBase()---->initialization: bubble intergroup transfer class: start");

    if (debugInfoAux) Info << "\n            - Number of bubble bins     : " << n_bins_ << endl;
    if (debugInfoAux) Info<< "\n            - bubble population radius (m)        : " << r_bins_ << endl;


    InfoPos("\n\n    intergroupTransfers::intergroupBase()---->FIELDS declaration/initialization-coalencence and breakup terms\n\n");


    // initialize the correspondence between ij and m indexes
    ij_init();

    // initialize intergroup transfer class fields
    forAll(r_bins_,i)
    { 
        char binIndex [5];
        sprintf(binIndex, "%02i", i+1);
        
        std::ostringstream iss;
        iss << "bin_" << binIndex;

        // intergroup source/sink term for each bin transport equation
        IG_.append( mesh_.lookupObjectRefPtr<volScalarField>("IG_"+ iss.str()) );

        // intermediate fields
        init_field(Xplus_, "Xplus_"+ iss.str(), dimIG);      // coalescence source
        init_field(Xminus_, "Xminus_"+ iss.str(), dimIG);    // coalescence sink
        init_field(Bplus_, "Bplus_"+ iss.str(), dimIG);      // breakup source
        init_field(Bminus_, "Bminus_"+ iss.str(), dimIG);    // breakup sink

        // averaged bubble radius - dimensioned scalar pointer list
        r_b_.append
        (
            new dimensionedScalar
            (
                "r_b_",
                dimLength,
                r_bins_[i]
            )
        );                        
    }

    // initiate BrKernal_ , BrXik_, Xikl_ scalar arrays
    init_arrays();

    // initialization of coalencence symmetric fields
    forAll(mi,m)
    {    
        int i = mi[m];
        int j = mj[m];

        std::ostringstream ijss;
        ijss << i << "_" << j;
        const string ij = ijss.str();

        init_field(T_, "T_"+ ij, dimIG);
        init_field(W_, "W_"+ ij, dimless);
        init_field(CE_, "CE_"+ ij, dimless);
    }
    InfoPos("\n\n    intergroupBase::intergroupBase()---->    calling update();");

    update();

    InfoPos("\n<---intergroupTransfers::intergroupBase()---->initialization: bubble intergroup transfer class: end");
}    


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::intergroupTransfers::intergroupBase::~intergroupBase()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// initiate BrKernal_ , BrXik_, Xikl_ scalar arrays
void Foam::intergroupTransfers::intergroupBase::init_arrays()
{
    // initiate BrKernal_
    InfoPos("\n\n    intergroupTransfers::init_arrays()---->    initiate BrKernal_");

    forAll(r_bins_,k)
    {
        if (debugInfoAux) Info << "                        k: " << k << endl;

        init_field(BrKernal_, "BrKernal_"+ k, IG_[0].dimensions());    
    }    



    // initiate BrXik_
    InfoPos("\n\n    intergroupTransfers::init_arrays()---->    initiate BrXik_");
    
    BrXik_ = new scalar*[n_bins_];

    forAll(r_bins_,i)
    {
        if (debugInfoAux) Info << "                        i: " << i << " k x " << n_bins_ << endl; 
        BrXik_[i] = new scalar[n_bins_];
    }
    // calculate the number of bubble transfered from the breakup of bubble into two identical daughter bubbles
    calc_BrXik();


    InfoPos("\n\n    intergroupTransfers::init_arrays()---->    initiate Xikl_");

    Xikl_ = new scalar**[n_bins_];

    forAll(r_bins_,i)
    {
        Xikl_[i] = new scalar*[n_bins_];

        forAll(r_bins_,k)
        {   
            if (debugInfoAux) Info << "                        i: " << i << " k: " << k << "  l x " << n_bins_<< endl;
            Xikl_[i][k] = new scalar[n_bins_];
        }
    }

    // calculate the number of bubble transfered from the calescence of two bubble from group k and l to group i
    calc_Xikl();
}


// initiate field
void Foam::intergroupTransfers::intergroupBase::init_field(PtrList<volScalarField>& F_, const word& fieldName_, const dimensionSet& dimset_)
{

    F_.append
    (
        new volScalarField
        (
            IOobject
            (
                fieldName_,                
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimset_
        )
    );
    // info index
    if (debugInfoAux) Info<< "            " << fieldName_ << "   initialized as volScalarField " <<
           "with dimensions: " << dimset_ << endl;

}

// calculate the collision rate of bubble groups k and l
void Foam::intergroupTransfers::intergroupBase::calc_T()
{

    const scalar& pi_ = Foam::constant::mathematical::pi;
    const volScalarField& epsilon_ = mesh_.lookupObject<volScalarField>("epsilon");

    forAll(T_,m)
    {
        if (debugInfoPos) Info<< "            calc of:  " << T_[m].name() << endl;

        int k = mi[m];
        int l = mj[m];

        // only fill the upper matrix and the diagonal because: T_ij=T_ji

        T_[m] = sqrt(2.0)/4.0 * pi_ * sqr(2.0*r_b_[k] + 2.0*r_b_[l])
                * pow(epsilon_, 1.0/3)
                * sqrt( pow(2*r_b_[k], 2.0/3) + pow(2*r_b_[l], 2.0/3) )
                * bubbleBinsList_[k].Nb() * bubbleBinsList_[l].Nb();
//                /  dimensionedScalar("volumeUnit", dimVolume, 1.0); <--- PLEASE VALIDATE
    }
}

// calculate the Weber number
// according to Luo (1993) and mainly Chen et al. (2005)
void Foam::intergroupTransfers::intergroupBase::calc_W()
{
    const volScalarField& epsilon_ = mesh_.lookupObject<volScalarField>("epsilon");

    forAll(W_,m)
    {
        if (debugInfoPos) Info<< "            calc of:  " << W_[m].name() << endl;

        int k = mi[m];
        int l = mj[m];

        // only fill the upper matrix and the diagonal because: W_ij=W_ji
        W_[m] = rho_water_ * 2 * r_b_[k]
                * ( beta_ * pow(epsilon_ * 2 * r_b_[k], 2.0/3 ) )
                * ( 1 + pow(r_b_[k]/r_b_[l], -2.0/3) )
                / sigma_;
    }
}

// calculate the number of bubble transfered from the calescence of two bubble from group k and l to group i
void Foam::intergroupTransfers::intergroupBase::calc_Xikl()
{
    InfoPos("\n\n    intergroupTransfers::calc_Xikl()---->");

    forAll(r_bins_,i)
    {
        forAll(r_bins_,k)
        {   
            forAll(r_bins_,l)
            {

                if ( i > 0 && k < i && l < i )
                {
                    // bubble volumes
                    scalar vi_minus_1 = 4*Foam::constant::mathematical::pi/3*pow3(r_bins_[i-1]);

                    scalar vi_plus_1 = 1e-15;  
                    if ( i < r_bins_.size()-1 )
                    {
                        vi_plus_1 = 4*Foam::constant::mathematical::pi/3*pow3(r_bins_[i+1]);
                    }

                    scalar vi = 4*Foam::constant::mathematical::pi/3*pow3(r_bins_[i]);
                    scalar vl = 4*Foam::constant::mathematical::pi/3*pow3(r_bins_[l]);    
                    scalar vk = 4*Foam::constant::mathematical::pi/3*pow3(r_bins_[k]);

                    if ( vi_minus_1 < vk+vl && vk+vl < vi )
                    {
                        Xikl_[i][k][l] = ( vk + vl - vi_minus_1)
                                         / (vi - vi_minus_1)
                                         * vi_minus_1;

                        if (debugInfoCalc) Info<< "            calc of:  Xikl_" << i << "_" << k << "_" << l << " = " << Xikl_[i][k][l] << endl;                                                 
                    }
                    else if ( i < r_bins_.size()-1 && vi < vk+vl && vk+vl < vi_plus_1 )
                    {
                        Xikl_[i][k][l] = ( vi_plus_1 - (vk-vl) )
                                         / (vi_plus_1 - vi)
                                         * vi;

                        if (debugInfoCalc) Info<< "            calc of:  Xikl_" << i << "_" << k << "_" << l << " = " << Xikl_[i][k][l] << endl;
                    }
                    else
                    {
                        Xikl_[i][k][l] = 0.0;
                        if (debugInfoCalc) Info<< "            calc of:  Xikl_" << i << "_" << k << "_" << l << " = " << Xikl_[i][k][l] << endl;                        
                    }
                }
                else
                {
                    Xikl_[i][k][l] = 0.0;
                    if (debugInfoCalc) Info<< "            calc of:  Xikl_" << i << "_" << k << "_" << l << " = " << Xikl_[i][k][l] << endl;
                }
            }    
        }
    }

}

// calculate the coalescence efficiency
//    represents the probability of coalescence when colision occurs
// based on Luo (1993)
void Foam::intergroupTransfers::intergroupBase::calc_CE()
{
    forAll(CE_,m)
    {
        if (debugInfoPos) Info<< "            calc of:  " << CE_[m].name() << endl;

        int k = mi[m];
        int l = mj[m];

        // only fill the upper matrix and the diagonal because: CE_ij=CE_ji
        CE_[m] = exp(
                    -sqrt(  0.75 
                            * ( 1 + sqr(r_b_[k]/r_b_[l]) )
                            * ( 1 + pow(r_b_[k]/r_b_[l],3) )
                         )
                    / (
                        sqrt( rho_air_/rho_ref_ + 0.5)
                        * pow( 1 + r_b_[k]/r_b_[l], 3)
                      )
                    * sqrt(W_[m])
                    );
    }
}


// calculate the coalencence source:
//     represent the gain in bubble group "i" due to coalescence of smaller bubbles.
// Prince and Blanch (1990)
void Foam::intergroupTransfers::intergroupBase::calc_Xplus()
{
    InfoPos("\n\n    intergroupTransfers::calc_Xplus()");

    forAll(r_bins_,i)
    {
        Xplus_[i] = 0.0*IGUnit;  //<--- please VALIDATE next line!

        forAll(r_bins_,k)
        {    
            forAll(r_bins_,l)
            {
                if (k < i && l < i )
                {

                    // symmetric matrix single index of T_ and CE_;
                    const int m = idx(k,l);
                    
                    Xplus_[i] += 0.5 * T_[m]*CE_[m]*Xikl_[i][k][l];
                }
            }
        }

        if (debugInfoCalc) Info<< "            calc of:  " << Xplus_[i].name() << "  sum= " << gSum(Xplus_[i]) << endl;
    }
}

// calculate the coalencence sink
void Foam::intergroupTransfers::intergroupBase::calc_Xminus( )
{
    InfoPos("\n    intergroupTransfers::calc_Xminus");

    forAll(r_bins_,i)
    {
        Xminus_[i] = 0.0*IGUnit;  //<--- please VALIDATE next line!

        forAll(r_bins_,k)
        { 
            // symmetric matrix single index of T_ and CE_;
            int m = idx(i,k);

            Xminus_[i] += T_[m]*CE_[m];
        }
        if (debugInfoCalc) Info<< "            calc of:  " << Xminus_[i].name() << "  sum= " << gSum(Xminus_[i]) << endl;
    }
}


// calculate the bubble breakup source
void Foam::intergroupTransfers::intergroupBase::calc_Bplus( )
{
    InfoPos("\n    intergroupTransfers::calc_Bplus");

    forAll(r_bins_,i)
    {

        Bplus_[i] = 0.0*IGUnit;  //<--- please VALIDATE next line!

        forAll(r_bins_,k)
        {
            Bplus_[i] += BrKernal_[k] * BrXik_[i][k];
        }

        if (debugInfoCalc) Info<< "            calc of:  " << Bplus_[i].name() << "  sum= " << gSum(Bplus_[i]) << endl;
    }
}


// calculate the bubble breakup sink
void Foam::intergroupTransfers::intergroupBase::calc_Bminus()
{
    InfoPos("\n    intergroupTransfers::calc_Bminus");

    forAll(r_bins_,i)
    {
        Bminus_[i] = BrKernal_[i];

        if (debugInfoCalc) Info<< "            calc of:  " << Bminus_[i].name() << "  sum= " << gSum(Bminus_[i]) << endl;
    }
}


// calcutate the breakup kernal function
// by Luo and Svendsen (1996)
// Luo, H., Svendsen, H.F., 1996. Theorethical model for drop and bubble breakup in turbulent dispersions. AIChE J. 42, 1225-1233.
void Foam::intergroupTransfers::intergroupBase::calc_BrKernal()
{
    InfoPos("\n\n    intergroupTransfers::calc_BrKernal()---->");

    const scalar cb = 0.923;
    const scalar cf = 0.2599;

    const volScalarField& epsilon_ = mesh_.lookupObject<volScalarField>("epsilon");
    const volScalarField& alpha_b_ = mesh_.lookupObject<volScalarField>("Cb");
    const volScalarField& alphaWater = mesh_.lookupObject<volScalarField>("alpha.water");


    forAll(r_bins_,k)
    {
        if (debugInfoPos) Info<< "            calc of:  BrKernal_" << k << flush;
        
        // dimensionless minimum eddy size
        // by van den Hengel et al. (2005)
        volScalarField ksi_min_ = 11.4 * pow( pow3(nu_water_)/epsilon_, 1.0/4.0) / (2*r_b_[k]);

            // Check "l_min_" dimensions consistency
        if ( ksi_min_.dimensions() != dimless )
        {
            FatalError
            << "    intergroupTransfers::calc_BrKernal()"
            << "    Check the dimesions of l_min_"
            << exit(FatalError);
        }

        scalarField ksi_min_smaller1 = pos(ksi_min_ - 1.0);

        if (debugInfoCalc) Info << "            # cells where 'ksi_min' >= 1: " << gSum(ksi_min_smaller1) << " of total " << mesh_.cells().size() << endl;


        volScalarField ksi_integral
        (
            IOobject
            (
                "ksi_integral_bin_" + k,                
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimless //IG_[0].dimensions()
        );

        // ksi_integral constant A
        volScalarField A_ = 12 * cf * sigma_
                            /   (
                                    beta_ * rho_water_ * pow(epsilon_, 2.0/3.0)
                                     * pow(2*r_b_[k], 5.0/3.0)
                                 );

        forAll(mesh_.cells(), celli)
        {
            scalar ksi_min_cell = ksi_min_[celli];

            if(ksi_min_cell < 1.0)
            {
                scalar a_ = A_[celli];

                // ksi_integral by with incomplete gamma functions
                // difference of the indef integral of 1, minus the indef integral of ksi_min
                ksi_integral[celli] =  ksi_indef_integral(1.0,a_) - ksi_indef_integral(ksi_min_cell,a_);  
            }
            else
            {
                ksi_integral[celli] = 0.0;
            }
        }

        ksi_integral.correctBoundaryConditions();


        BrKernal_[k] = cb
                       * ( alphaWater - alpha_b_ ) // option 3: consdider the vof and the bubble phase
                       * bubbleBinsList_[k].Nb()
                       * pow( epsilon_ / ( 4 * sqr(r_b_[k]) ), 1.0/3.0)
                       * ksi_integral;
    }

}

// indefinite integral of the integral in the breakage kernal expression of Luo and Svendsen (1996)
Foam::scalar Foam::intergroupTransfers::intergroupBase::ksi_indef_integral(const scalar ksi_, const scalar a_)
{
    // according to http://www.wolframalpha.com
    scalar ksi_indef_integral_ = 3*ksi_*pow(a_/pow(ksi_,11.0/3),3.0/11)
                                    *(
                                        2*ksi_*pow(a_/pow(ksi_,11.0/3),3.0/11)*boost::math::tgamma(5.0/11,a_/pow(ksi_,11.0/3))
                                        + boost::math::tgamma(8.0/11,a_/pow(ksi_,11.0/3))
                                        + pow(ksi_,2.0)*pow(a_/pow(ksi_,11.0/3),6.0/11)*boost::math::tgamma(2.0/11,a_/pow(ksi_,11.0/3))
                                    )
                                 /(11*a_);

    return ksi_indef_integral_;
}

// function inside the integral in the breakage kernal expression of Luo and Svendsen (1996)
Foam::scalar Foam::intergroupTransfers::intergroupBase::ksi_func(const scalar ksi_, const scalar a_)
{
    scalar ksi_func_val = sqr(1+ksi_) / pow(ksi_, 11.0/3.0)
                       * exp( - a_ * pow(ksi_, -11.0/3.0) );

    return ksi_func_val;
}

// function inside the integral in the breakage kernal expression of Luo and Svendsen (1996)
Foam::scalar Foam::intergroupTransfers::intergroupBase::ksi_func_RiemannIntegral(const scalar a_, const scalar xmin_, const scalar xmax_, const int N_)
{
    // Riemann integral of the following function:
    const scalar step_ = ( xmax_ - xmin_) / N_;
    scalar integral_ = 0.0;

    if (debugInfoCalc) Info<< "                xmin = " << xmin_ << "  to xmax_ = " << xmax_ << endl;

    if( xmax_ > xmin_)
    {
        if (debugInfoCalc) Info<< "                interval lenght = " << xmax_ - xmin_ << "                dscretization intervals = " << N_ << endl;

        // x coordinate
        scalar x_ = xmin_ + step_/2.0;

        for(int i=0; i<N_; i++)
        {
            // x coordinate update
            if (debugInfoPos) Info << "------>intergroupTransfers::ksi_func_RiemannIntegral(): (i-1)=(" << i << "-1)=" << (i-1) << endl;
            x_ += (i-1) * step_;

            integral_ += ksi_func(x_, a_) * step_;
        }        

        if (debugInfoCalc) Info<< "                Riemann integral of ksi_func(x,a) with " << N_ << " intervals = " << integral_ << endl;


        if (debugInfoCalc) Info<< "\n                ksi_func(xmin) = " << ksi_func(xmin_, a_) << "      ksi_func(xmax) = " << ksi_func(xmax_, a_) << endl;
        if (debugInfoCalc) Info<< "                (ksi_func(xmin)+ksi_func(xmax))/2 * (xmax - xmin) = " << (ksi_func(xmin_, a_)+ksi_func(xmax_, a_))/2 * (xmax_ - xmin_) << endl;

    }
    else
    {
        Info<< "                Riemann integral NOT CALCULATED: ksi_min > 1.0 " << endl;
    }
    return  integral_;
}

// calculate the number of bubble transfered from the breakup of bubble into two identical daughter bubbles
void Foam::intergroupTransfers::intergroupBase::calc_BrXik()
{
    InfoPos("\n\n    intergroupTransfers::calc_BrXik()---->");

    forAll(r_bins_,i)
    {
        forAll(r_bins_,k)
        {   
            if ( i > 0 )
            {
                // bubble volumes
                scalar vi_minus_1 = 4*Foam::constant::mathematical::pi/3*pow3(r_bins_[i-1]);

                scalar vi_plus_1 = 1e-15;  
                if ( i < r_bins_.size()-1 )
                {
                    vi_plus_1 = 4*Foam::constant::mathematical::pi/3*pow3(r_bins_[i+1]);
                }

                scalar vi = 4*Foam::constant::mathematical::pi/3*pow3(r_bins_[i]);
                scalar vk = 4*Foam::constant::mathematical::pi/3*pow3(r_bins_[k]);


                if ( vi_minus_1 < vk/2.0 && vk/2.0 < vi )
                {
                    BrXik_[i][k] = 2 * ( vk/2.0 - vi_minus_1)
                                     / (vi - vi_minus_1);

                    if (debugInfoCalc) Info<< "            calc of:  BrXik_" << i << "_" << k << " = " << BrXik_[i][k] << endl;                                                 
                }
                else if ( i < r_bins_.size()-1 && vi_plus_1 > vk/2.0 && vk/2.0 > vi )
                {
                    BrXik_[i][k] = 2 * ( vi_plus_1 - vk/2.0)
                                     / (vi_plus_1 - vi_minus_1);

                    if (debugInfoCalc) Info<< "            calc of:  BrXik_" << i << "_" << k << " = " << BrXik_[i][k] << endl;
                }
                else
                {
                    BrXik_[i][k] = 0.0;
                    if (debugInfoCalc) Info<< "            calc of:  BrXik_" << i << "_" << k << " = " << BrXik_[i][k] << endl;                        
                }
            }
            else
            {
                BrXik_[i][k] = 0.0;
                if (debugInfoCalc) Info<< "            calc of:  BrXik_" << i << "_" << k << " = " << BrXik_[i][k] << endl;
            }
        }    
    }


}

// calculate the source term for bubble intergroup transfer in the number of bubble eq.
void Foam::intergroupTransfers::intergroupBase::calc_IG()
{
    InfoPos("\n\n    intergroupTransfers::calc_IG_");

    forAll(r_bins_,i)
    {

        IG_[i] = Xplus_[i] - Xminus_[i] + Bplus_[i] - Bminus_[i];

        if (debugInfoCalc) Info<< "            calc of:  " << IG_[i].name() << "  = " << gSum(IG_[i]) << endl;
    }
}

// update all fields
void Foam::intergroupTransfers::intergroupBase::update()
{
    // calculate the collision rate of bubble groups
    calc_T();
    // calculate the coalescence efficiency
    calc_CE();
    // calculate the Weber number
    calc_W();
    // calculate the coalencence source
    calc_Xplus( );
    // calculate the coalencence sink
    calc_Xminus( );
    // calcutate the breakup kernal function
    calc_BrKernal();
    // calculate the bubble breakup source
    calc_Bplus( );
    // calculate the bubble breakup sink
    calc_Bminus( );
    // calculate the bubble number equation source term
    calc_IG();
}


// ************************************************************************* //