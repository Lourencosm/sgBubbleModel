/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField
::variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    flowRate_(),
    alphaName_("none"),
    air_bl_height_(0),
    air_bl_zMax_(0),
    UairRef_(0),
    buffer_height_(0),
    water_bl_height_(0)
{}


Foam::variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField
::variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    flowRate_(Function1<scalar>::New("flowRate", dict)),
    alphaName_(dict.lookup("alpha")),
    air_bl_height_(readScalar(dict.lookup("airBlHeight"))),
    air_bl_zMax_(readScalar(dict.lookup("airBlTopLevel"))),
    UairRef_(readScalar(dict.lookup("airFreestreamUmag"))),
    buffer_height_(readScalar(dict.lookup("bufferLayerHeight"))),
    water_bl_height_(readScalar(dict.lookup("water_bl_height")))
{}


Foam::variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField
::variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField
(
    const variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    flowRate_(ptf.flowRate_.clone()),
    alphaName_(ptf.alphaName_),
    air_bl_height_(ptf.air_bl_height_),
    air_bl_zMax_(ptf.air_bl_zMax_),
    UairRef_(ptf.UairRef_),
    buffer_height_(ptf.buffer_height_),
    water_bl_height_(ptf.water_bl_height_)
{}


Foam::variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField
::variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField
(
    const variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    flowRate_(ptf.flowRate_.clone()),
    alphaName_(ptf.alphaName_),
    air_bl_height_(ptf.air_bl_height_),
    air_bl_zMax_(ptf.air_bl_zMax_),
    UairRef_(ptf.UairRef_),
    buffer_height_(ptf.buffer_height_),
    water_bl_height_(ptf.water_bl_height_)
{}


Foam::variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField
::variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField
(
    const variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    flowRate_(ptf.flowRate_.clone()),
    alphaName_(ptf.alphaName_),
    air_bl_height_(ptf.air_bl_height_),
    air_bl_zMax_(ptf.air_bl_zMax_),
    UairRef_(ptf.UairRef_),
    buffer_height_(ptf.buffer_height_),
    water_bl_height_(ptf.water_bl_height_) 
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField
::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalarField alphap =
        patch().lookupPatchField<volScalarField, scalar>(alphaName_);

    alphap = max(alphap, scalar(0));
    alphap = min(alphap, scalar(1));

    const scalar t = db().time().timeOutputValue();
    scalar flowRate = flowRate_->value(t);
    
    scalar UairFS = -UairRef_;

    // patch centre faces z coordinate
    scalarField z = patch().Cf().component(vector::Z);

    // patch boundaries - OPTION A
    const vectorField& patchPoints = patch().patch().localPoints();    
    boundBox bb(patchPoints, true);
    const scalar bbzMin = bb.min().z();
    const scalar bbzMax = bb.max().z();


    // water bottom boundary layer factor
    // initialization
    scalarField wbl_factor = alphap*0.0;

    forAll(wbl_factor, facei)
    {
        if (z[facei]-bbzMin <= water_bl_height_)
        {
            wbl_factor[facei] = pow((z[facei]-bbzMin)/water_bl_height_, (1.0/7.0));            
        }
        else
        {	
        	wbl_factor[facei] = 1.0;
        }
    }


    // a simpler way of doing this would be nice
    scalar UwaterFS = -flowRate/gSum(patch().magSf()*alphap*wbl_factor);

    vectorField n(patch().nf());

    // water level - approximate
    scalar wetArea = gSum(patch().magSf()*alphap);
    scalar dryArea = gSum(patch().magSf());

    scalar wl = wetArea/dryArea*(bbzMax-bbzMin)+bbzMin;

    // air boundary layer height = max( input height, water level - input level)
    scalar abl_h = max(air_bl_height_,air_bl_zMax_-wl);

    // air boundary layer factor
    // initialization
    scalarField Umag = alphap*0.0;

    forAll(Umag, facei)
    {
    	// water boundary layer
    	if (z[facei]-bbzMin <= water_bl_height_)
    	{
    		Umag[facei] = UwaterFS*wbl_factor[facei];
    	}
    	// water freestream
    	else if ( (z[facei]-bbzMin > water_bl_height_) && (z[facei] <= wl) )
    	{
    		Umag[facei] = UwaterFS;
    	}
    	// air buffer layer above water    	    	
        else if ( (z[facei] > wl) && (z[facei] < wl+buffer_height_) )
        {
            Umag[facei] = UwaterFS-((UwaterFS-UairFS)/2)*pow((z[facei]-wl)/buffer_height_, (3.0)); //7.0
        }
    	// air boundary layer above water    	    	
        else if ( (z[facei] > wl+buffer_height_) && (z[facei] < wl+buffer_height_+abl_h) )
        {
            Umag[facei] = (UwaterFS-UairFS)/2-((UwaterFS-3*UairFS)/2)*pow((z[facei]-(wl+buffer_height_))/abl_h, (1.0/7.0));
        }
        else if (z[facei] >= wl+buffer_height_+abl_h)
        {
            Umag[facei] = UairFS;
        }
    }
    operator==(n*Umag);

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchField<vector>::write(os);
    flowRate_->writeData(os);
    os.writeEntry("alpha", alphaName_);
    os.writeEntry("airBlHeight", air_bl_height_);
    os.writeEntry("airBlTopLevel", air_bl_zMax_);
    os.writeEntry("airFreestreamUmag", UairRef_);
    os.writeEntry("bufferLayerHeight", buffer_height_);        
    os.writeEntry("water_bl_height", water_bl_height_);
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField
   );
}


// ************************************************************************* //
