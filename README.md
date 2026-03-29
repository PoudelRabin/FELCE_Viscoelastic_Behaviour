# FELCE_Viscoelastic_Behaviour
The repository contains finite element simulation code assocaited with following paper:

"A numerical model for the viscoelastic behaviour of liquid crystal elastomers"

Authors: Rabin Poudel, Yasemin Sengul, Angela Mihai

**Note**: This repository is for hosting code of the published paper and is not intended for regular changes and updates.
## Folders
### src
This folder contains the sorce code for the finite element simualtion. The simulations are done in the open software [FEBio](https://febio.org/) ([FEBio Github](https://github.com/febiosoftware/FEBio)). FEBio 4.4 is used for the simualtion. Most of the codes are the FEBio 4.4 with user codes written by authors starts with the 'FELCE' and are in the folder src/FEBioMech. Some of the original files from src/FEBioMech are edited to register the FELCE codes.

### FEMaterialModel
This folder contains material code for the viscoelastic phenomema described in the paper. The Folders contain the c++ file "FELCEMaterial444" which need to be copied, moved to src/FEBioMech replacing the original according the phenomena that need to be simulated. The c++ file "FELCEInitialAngle.cpp" needs to be edited for the required initial angle and moved to src/FEBioMech.

### Examples
This folder contains the input file (.feb) and output (.xpt) for simualting the LCE material. Read the README file in the folder.

## Building and running the model
Building and running is similar to the one building the original FEBio excutables. Please follow the README.md and BUILD.md in the src folder as well as read the documentation provided. More resources are availabe in [FEBio](https://febio.org/) ([FEBio Github](https://github.com/febiosoftware/FEBio)). Please note that software requres the pardiso solver for solving the equations. Please install intel oneAPI math kernel library.
