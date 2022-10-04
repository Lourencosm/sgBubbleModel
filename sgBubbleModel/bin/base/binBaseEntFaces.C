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

#include "binBaseEntFaces.H"

#if OPENFOAM >= 1812
    #include "gravityMeshObject.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace bins
{
    defineTypeNameAndDebug(binBaseEntFaces, 0);
    addToRunTimeSelectionTable(bin, binBaseEntFaces, dictionary);    
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bins::binBaseEntFaces::binBaseEntFaces
(
            const fvMesh& mesh,
            const dictionary& bubbleBinDict   
)
:
    // base class
    bin(mesh, bubbleBinDict),
    // bubble population name
    name_(bubbleBinDict.lookup("binName")),
    /*
    // fields declaration
    */
    // number of bubbles in each cell
    Nb_
    (
        IOobject
        (
            "Nb_" + name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    // bubble velocity field
    Ub_
    (
        IOobject
        (
            "Ub_" + name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE 
        ),
        mesh_,
        dimVelocity
    ),
    // face flux field phi_b of air bubbles in water
    phi_b_
    (
        IOobject
        (
            "phi_b_" + name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(Ub_) & mesh_.Sf()
    ),
    // bubble difusion constant
    Db_
    (
        IOobject
        (
            "Db_" + name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE 
        ),
        mesh_,
        dimensionSet(0, 2, -1, 0, 0, 0, 0)
    ),
    // surfaceVectorField: entrainment of bubble number field @ faces
    En_faces_
    (
        IOobject
        (
            "En_faces_" + name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE 
        ),
        mesh_,
        dimensionedVector("zero", Nb_.dimensions()/(dimTime*dimArea), Zero)                  
    ),
    // volScalarField: entrainment of bubble number field
    En_
    (
        IOobject
        (
            "En_" + name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE 
        ),
        mesh_,
        dimensionedScalar("zero", Nb_.dimensions()/dimTime, Zero)       
    ),    
    // source term for bubble intergroup transfer in the number of bubble eq.
    IG_
    (
        IOobject
        (
            "IG_" + name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE 
        ),
        mesh_,
        dimensionedScalar("zero", Nb_.dimensions()/dimTime, 0.0)        
    ),
    // number of bubbles bubble evaporated to air field
    Nb_evap_
    (
        IOobject
        (
            "Nb_evap_" + name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE 
        ),
        mesh_,
        dimensionedScalar("zero", Nb_.dimensions(), 0.0)
    ),      
    // volume of bubble field
    Vb_
    (
        IOobject
        (
            "Vb_" + name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE 
        ),
        mesh_,
        dimensionedScalar("zero", dimVolume, 0.0)
    ),
    // volume of evaporated bubble field
    Vb_evap_
    (
        IOobject
        (
            "Vb_evap_" + name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE 
        ),
        mesh_,
        dimensionedScalar("zero", dimVolume, 0.0)
    ),
    // air concentraion in cells
    C_
    (
        IOobject
        (
            "C_" + name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
//        dimless
    ),
    // bubble evaporated to air field
    C_evap_
    (
        IOobject
        (
            "C_evap_" + name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE            
        ),
        mesh_,
        dimensionedScalar("zero", C_.dimensions(), 0.0)
    ),
    /*
    Shi et al 2010 - model variables
    */
    // air bubbles turbulent diffusivity
    // production term constant to be determined
    a_b_
    (
        "a_b",
//        dimensionSet(-1, -1, 2, 0, 0, 0, 0),<- for Shi production
        En_.dimensions()/(mesh_.lookupObject<volScalarField>("Pr").dimensions()*inv(dimVolume)),
        bubbleBinDict.lookup("ab")
    ),
    // production threshold for the onset air entrainment
    Pr0_threshold_
    (
        "Pr0_threshold_",
        mesh_.lookupObject<volScalarField>("Pr").dimensions(),
        bubbleBinDict.lookup("productionThreshold")
    ),
    // air bubbles maximum volumetric fraction or entrainment and for bounding
    Cb_max_(readScalar(bubbleBinDict.lookup("bubbleFractionMax"))),
    // surface alpha thresholds for the air entrainment evaportation region
    EnAlphaMin_(readScalar(bubbleBinDict.lookup("EnAlphaMin"))),
    // Schmidt number of air bubbles in water (Buscaglia et al. 2002)    
    Sch_(readScalar(bubbleBinDict.lookup("bubbleSchmidtNumber"))),
    // averaged bubble radius
    r_b_
    (
        "r_b_" + name_,
        dimLength,
        bubbleBinDict.lookup("bubbleRadius")
    ),
    // averaged bubble radius spacing    
    dr_b_
    (
        "dr_b_" + name_,
        dimLength,
        bubbleBinDict.lookup("radiusSpacing")
    ),    
    // volume of one bubble of this bin
    VOneBubble_ 
    (
        "VBubble_" + name_,
        4.0*Foam::constant::mathematical::pi/3.0*pow3(r_b_)
    ),
    // minimum boundary for Nb_
    Nb_min_
    (
        "Nb_min",
        Nb_.dimensions(),
        0.0
    ),
    // averaged bubble distribution    
    Di_b_
    (
        "Di_b_" + name_,
        inv(dimVolume),
        bubbleBinDict.lookup("bubbleDensity")
    ),
    // fvOptions
    fvOptions(fv::options::New(mesh_))

{
    Info<< "\n        Init bubble " << name_ << endl;
    Info<< "            - radius (m)                        : " << r_b_ << endl;
    Info<< "            - distribution radius spacing (m)   : " << dr_b_ << endl;
    Info<< "            - density number PF (m^-3)          : " << Di_b_ << endl; 

    // update bubble volume and bubble number 
    update_Vb_C();

    InfoPos("\n        bins::binBaseEntFaces::bubbleBin()---->initializated ");
}    


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bins::binBaseEntFaces::~binBaseEntFaces()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
const Foam::word& Foam::bins::binBaseEntFaces::name() const
{
    return name_;
}

//- Return const-access to bubbleBin diffusion coefficient
const Foam::dimensionedScalar& Foam::bins::binBaseEntFaces::P_r0_threshold() const
{
    return Pr0_threshold_;
}

//- Return const-access to bubbleBin diffusion coefficient
const Foam::scalar& Foam::bins::binBaseEntFaces::EnAlphaMin() const
{
    return EnAlphaMin_;
}

//- Return const-access to bubbleBin radius
const Foam::dimensionedScalar& Foam::bins::binBaseEntFaces::r_b() const
{
    return r_b_;
}

//- Return const-access to bubbleBin averaged bubble radius spacing
const Foam::dimensionedScalar& Foam::bins::binBaseEntFaces::dr_b() const
{
    return dr_b_;
}        

//- Return const-access to bubbleBin averaged bubble radius spacing
const Foam::dimensionedScalar& Foam::bins::binBaseEntFaces::VOneBubble() const
{
    return VOneBubble_;
} 

const Foam::volScalarField& Foam::bins::binBaseEntFaces::C() const
{
    return C_;
}

// field
const Foam::volScalarField& Foam::bins::binBaseEntFaces::C_evap() const  // <-- to remove
{
    return C_evap_;
}

// field
const Foam::volScalarField& Foam::bins::binBaseEntFaces::Nb() const
{
    return Nb_;
}

// field
const Foam::volScalarField& Foam::bins::binBaseEntFaces::Nb_evap() const
{
    return Nb_evap_;
}

// field
const Foam::volScalarField& Foam::bins::binBaseEntFaces::Vb() const
{
    return Vb_;
}

// field
const Foam::volScalarField& Foam::bins::binBaseEntFaces::Vb_evap() const
{
    return Vb_evap_;
}
// field
const Foam::volScalarField& Foam::bins::binBaseEntFaces::En() const
{
    return En_;
}

// field
const Foam::surfaceVectorField& Foam::bins::binBaseEntFaces::En_faces() const
{
    return En_faces_;
}

// field
const Foam::volVectorField& Foam::bins::binBaseEntFaces::Ub() const
{
    return Ub_;
}

// field
const Foam::surfaceScalarField& Foam::bins::binBaseEntFaces::phi_b() const
{
    return phi_b_;
}                                
// air bubbles turbulent diffusivity
const Foam::volScalarField& Foam::bins::binBaseEntFaces::Db() const
{
    return Db_;
}  
// term for bubble intergroup transfer in the number of bubble eq.
const Foam::volScalarField& Foam::bins::binBaseEntFaces::IG() const
{
    return IG_;
} 

//- Return const-access to bubbleBin diffusion coefficient
const Foam::dimensionedScalar& Foam::bins::binBaseEntFaces::a_b() const
{
    return a_b_;
}

// update bubble velocity
void Foam::bins::binBaseEntFaces::update_Ub()
{

    const volVectorField& U_ = mesh_.lookupObject<volVectorField>("U");

    const scalar& r_b__ = r_b_.value();

    // read gravitationalaceleration
    #if OPENFOAM < 1812
        const uniformDimensionedVectorField& g_ = mesh_.lookupObject<uniformDimensionedVectorField>("g");
    #endif
    #if OPENFOAM >= 1812
        const uniformDimensionedVectorField& g_ = meshObjects::gravity::New(mesh_.time());
    #endif

    // gravity unit opposite vector (pointing the sky)
    const dimensionedVector g_b_unit = (-1.0*g_/mag(g_));

    // mesh_: openfoam solver mesh
    // gv_: gravity unit opposite vector (pointing the sky)
    scalar ws_ = 0.0;    // bubble-slip velocity, Clift et al. (1978)

    // bubble-slip velocity, Clift et al. (1978)
    if ( 0.0 <= r_b__ && r_b__ <= 7e-4 ) {
    ws_ = 4474.0*Foam::pow(r_b__,1.357);
    //Info<< "ws 1" << endl;
    }
    else if ( 7e-4 < r_b__ && r_b__ <= 5.1e-3 ) {
    ws_ = 0.23;
    //Info<< "ws_ 2" << endl;
    }
    else if ( 5.1e-3 < r_b__ ) {
    ws_ = 4.202*Foam::pow(r_b__,0.547);
    //Info<< "ws_ 3" << endl;
    }

    // magnitude of bubble extra velocity parallel to the gravitacion vector
    dimensionedScalar Ub_extraUvert
    (
            "Ub_extraUvert",
            dimensionSet(0, 1, -1, 0, 0, 0, 0),
            ws_
    );

    // bubble velocity
    Ub_ =  U_ + Ub_extraUvert*g_b_unit;

    // face flux field phi_b of air bubbles in water
    phi_b_ = linearInterpolate(Ub_) & mesh_.Sf();

    if (debugInfoVar) Info<<"        Ub updated with " << "        ws_=" << ws_ << endl;
}


// calculate air entranment in faces - method B
void Foam::bins::binBaseEntFaces::entrainment()
{

    const surfaceVectorField& entrainmentMask_faces_vector = mesh_.lookupObject<surfaceVectorField>("entrainmentMask_faces_vector");
    const volScalarField& Pr_ = mesh_.lookupObject<volScalarField>("Pr");

    // production term field magnitude @ faces dimensions=[1 -1 -3]
    surfaceScalarField Pr_faces = linearInterpolate(Pr_);
    // Pr_faces.dimensions=[1 -1 -3 ...]
    if (debugInfoCalc) Info<< "            Pr_faces dimensions     : " << Pr_faces.dimensions() << endl;

    // production term field VECTOR @ faces
    surfaceVectorField PrSurf_faces = Pr_faces*entrainmentMask_faces_vector;
    // PrSurf_faces.dimensions=[1 -1 -3 ...]
    if (debugInfoCalc) Info<< "            PrSurf_faces dimensions     : " << PrSurf_faces.dimensions() << endl;

    // dimension scalar to fix the a_b_ dimensions from base class
    const dimensionedScalar a_b_dim_fix("a_b_dim_fix", pow(dimLength,-2) , 1.0);

    // surfaceVectorField: air entrainment vector field @ faces
    En_faces_ = pos(Pr_faces - Pr0_threshold_)*PrSurf_faces*a_b_*a_b_dim_fix*Di_b_;
    // En_faces_.dimensions=[0 -5 -1 ...]    
    if (debugInfoCalc) Info<< "            En_faces_ dimensions     : " << En_faces_.dimensions() << endl;

    // air entranment flux sum in each face
	surfaceScalarField En_faces_flux = (-1)*(En_faces_ & mesh_.Sf());
    // En_faces_flux.dimensions=[0 -3 -1 ...]     
    if (debugInfoCalc) Info<< "            En_faces_flux dimensions : " << En_faces_flux.dimensions() << endl; 

    // air entranment flux sum of all cell faces @ cell center
    En_ = fvc::surfaceSum(En_faces_flux);
    // En_.dimensions=[0 -3 -1 ...]     

    // mask: keep only the entrainment increment
    En_ *= pos(En_);

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
        Info<< "            prod @ surface, mag(PrSurf_faces):            " << mag(gMin(PrSurf_faces)) << " " << mag(gMax(PrSurf_faces)) << " " << mag(gAverage(PrSurf_faces)) << endl;
        Info<< "            entrainment, En_faces_       :\n            " << gMin(En_faces_) << " " << gMax(En_faces_) << " " << gAverage(En_faces_) << endl;

        InfoVar(En_faces_flux);
        InfoVar(En_);
    }

}


// solve bubble transport equation
void Foam::bins::binBaseEntFaces::solve_Nb()
{
    if (debugInfoPos) Info<< "\n    bins::binBaseEntFaces::solve_Nb()---->solving transport...     " << name_ << endl << endl;

    const volScalarField& alphaWater = mesh_.lookupObject<volScalarField>("alpha.water");
    const volScalarField& nut_ = mesh_.lookupObject<const volScalarField>("nut");
    
    // Schmidt number of air bubbles in water (Buscaglia et al. 2002)
    //1.0 - theoric // 0.7 - shi experiments

    // air bubbles turbulent diffusivity    
    Db_ = nut_/Sch_; 


    InfoPos("\n    bins::binBaseEntFaces::solve_Nb()---->before Nb_eqn");

    // Number of bubbles of air transport equation
    fvScalarMatrix Nb_Eqn
    (
        fvm::ddt(Nb_)
     +  fvm::div(phi_b_, Nb_)
     -  fvm::laplacian(Db_, Nb_)
    ==
        fvOptions(Nb_)
    );

    Nb_Eqn.relax();

    fvOptions.constrain(Nb_Eqn);

    InfoPos("\n    bins::binBaseEntFaces::solve_Nb()---->before solve");
    
    solve
    (
        Nb_Eqn
     ==
        En_
     +  IG_
    );

    if (debugInfoCalc)
    {
        Info<< "\n        after Nb_Eqn:" << endl;
        InfoVar(Nb_);
        Info<< "                  Sum bin Nb_      : " << gSum(Nb_) << endl;       
    }

    fvOptions.correct(Nb_);

    if (debugInfoCalc)
    {
        Info<< "\n        after correct Nb_:" << endl;
        InfoVar(Nb_);
        Info<< "                  Sum bin Nb_      : " << gSum(Nb_) << endl;
    }

    bound(Nb_, Nb_min_);

    if (debugInfoCalc)
    {
        Info<< "\n        after bound Nb_:" << endl;
        InfoVar(Nb_);
        Info<< "                  Sum bin Nb_      : " << gSum(Nb_) << endl;   
    }
    // updating other fields

    // correction of Nb_ - evaporation
    if (debugInfoCalc)
    {
        Info<< "\n        before Nb_evap calc:" << endl;
        InfoVar(alphaWater);
        InfoVar(neg(alphaWater - EnAlphaMin_));        
        Info<< "                  EnAlphaMin_      : " << EnAlphaMin_ << endl;   
    }    

    Nb_evap_ = Nb_ * neg(alphaWater - EnAlphaMin_);
    Nb_ -= Nb_evap_;

    InfoVar(Nb_);
    InfoVar(Nb_evap_);

    if (debugInfoCalc)
    {
        Info<< "                  Sum bin Nb_      : " << gSum(Nb_) << endl; 
        Info<< "                  Sum bin Nb_evap_ : " << gSum(Nb_evap_) << endl;      
    }

    // update bubble volume and bubble number 
    update_Vb_C();
}

// calculate bubble volume and bubble fraction
void Foam::bins::binBaseEntFaces::update_Vb_C()
{
    // bubble fraction @ cell 
    C_.ref() = Nb_() * VOneBubble_;
//    C_.correctBoundaryConditions();

    // bubble fraction @ cell 
    C_evap_.ref() = Nb_evap_() * VOneBubble_;
//    C_evap_.correctBoundaryConditions();

    // bubble volume
    Vb_.ref() = C_.ref() * mesh_.V(); // / dimensionedScalar("volUnits", dimVolume, 1.0);

    // evaporated bubble volume
    Vb_evap_.ref() = C_evap_.ref() * mesh_.V(); // / dimensionedScalar("volUnits", dimVolume, 1.0);


    InfoVar(C_);
    InfoVar(C_evap_);
    InfoVar(Vb_);
    InfoVar(Vb_evap_);
}

// bounding bubble bins fields by a factor in each cell
void Foam::bins::binBaseEntFaces::Cb_bounding(const volScalarField& Cb_bounder)
{
    // volumetric fraction of bubble in cells 
    C_ /= Cb_bounder; 
    // volumetric fraction of evaportated bubble in cells
    C_evap_ /= Cb_bounder;
    // number of bubbles in each cell
    Nb_ /= Cb_bounder;
    // number of evaportaed bubbles in each cell
    Nb_evap_ /= Cb_bounder;
}

// ************************************************************************* //