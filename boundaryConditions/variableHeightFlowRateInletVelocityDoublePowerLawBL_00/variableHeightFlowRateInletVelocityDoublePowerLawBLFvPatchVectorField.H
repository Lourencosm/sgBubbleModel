/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    This boundary condition provides a velocity boundary condition for
    multphase flow based on a user-specified volumetric flow rate.

    The flow rate is made proportional to the phase fraction alpha at each
    face of the patch and alpha is ensured to be bound between 0 and 1.

    modified by Lourenco Sassetti

Usage
    \table
        Property     | Description             | Required    | Default value
        flowRate     | volumetric flow rate [m3/s] | yes |
        alpha        | phase-fraction field    | yes |
        airBlHeight     | boundary layer height [m] | yes | 10
        airBlTopLevel   | boundary top level [m] | yes |        
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            variableHeightFlowRateInletVelocity;
        flowRate        0.2;
        alpha           alpha.water;
        value           uniform (0 0 0); // placeholder
    }
    \endverbatim

    The \c flowRate entry is a \c Function1 of time, see Foam::Function1Types.

Note
    - the value is positive into the domain
    - may not work correctly for transonic inlets
    - strange behaviour with potentialFoam since the momentum equation is
      not solved

    The air boundary layer height is the maximum of \c blHeight entry \c and \c  blTopLevel entry \c - water level.
Note
    - at one of this must be specified. 

See also
    Foam::fixedValueFvPatchField
    Foam::Function1Types

SourceFiles
    variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField_H
#define variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
    Class variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Inlet integral flow rate
        autoPtr<Function1<scalar>> flowRate_;

        //- Name of the phase-fraction field
        word alphaName_;

        //- air boundary layer height
        scalar air_bl_height_;

        //- air boundary layer top level
        scalar air_bl_zMax_;

        //- air freestream velocity
        scalar UairRef_;
        
        //- buffer layer height above water surface
        scalar buffer_height_;
        
        //- water boundary layer height
        scalar water_bl_height_;
        
public:

   //- Runtime type information
   TypeName("variableHeightFlowRateInletVelocityDoublePowerLawBL");


   // Constructors

        //- Construct from patch and internal field
        variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField
        //  onto a new patch
        variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField
        (
            const variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField
        (
            const variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField
        (
            const variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new variableHeightFlowRateInletVelocityDoublePowerLawBLFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
