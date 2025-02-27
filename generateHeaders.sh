#!/bin/bash

headers=(
    "liftFoam.C"
    "createAdditionalFields.C"
    "LIFTModel.C"
    "twoTemperatureModel.C"
    "femtosecondLaserModel.C"
    "ultraFastShockWaveModel.C"
    "dropletModel.C"
    "extremeConditionMaterialProperties.C"
    "advancedInterfaceCapturing.C"
    "phaseChangeModel.C"
    "LIFTConfig.C"

)

# Create include directory
mkdir -p include

cd include

# Generate header files
for header in "${HEADERS[@]}"; do
    header_file="${header}.H"
    
    if [ ! -f "$header_file" ]; then
        echo "Creating $header_file..."
        
        # Generate header content
        cat > "$header_file" << EOF
/*---------------------------------------------------------------------------*\\
  =========                 |
  \\\\      /  F ield        | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration    |
    \\\\  /    A nd          | www.openfoam.com
     \\\\/     M anipulation |
-------------------------------------------------------------------------------
    Description
        Header file for LIFT (Laser-Induced Forward Transfer) process solver
        Part of the $header module
-------------------------------------------------------------------------------*/

#ifndef ${header}_H
#define ${header}_H

#include "fvCFD.H"

namespace Foam
{
    // Class declaration will go here
}

#endif
EOF
        echo "Created $header_file"
    fi
done

# Make script executable
chmod +x generateHeaders.sh

echo "Header generation complete"
