## Molecular Dynamics and Monte Carlo simulations of semiconductor samples at the nanoscale
This project is part of my dissertation for the doctorate diploma in Physics.
Its scope is to simulate and model the operation of a virtual Atomic Force Microscope nanotip
when interacting with semiconducting nanosamples (e.g. Si, Ge etc...)

### helper_scripts
This folder contains some pre- or post-processing utilities in order to manipulate 
the nanostructures atomic positions.

### MolDyn
This folder contains the methodology for the thermal relaxation of nanostructures.

### MonteCarlo
Alternatively one can use the Monte Carlo methodology for thermal relaxation.

### MolDyn_tip_surface
This folder contains the methodology for calculating the tip-sample interactions of 
an AFM nano-tip with a nanosurface.

We have used:
- Andersen thermostat
- Simulated anealing
- Lennard Jones potential
- Tersoff potential
- Linked cells optimization method
- Jmol and Paraview visualization

### How-To Build the projects
These projects run in all platforms (Windows, MacOS, Linux) so for a unified makefile solution we use cmake.
In order to build e.g. MolDyn project follow these steps:
- Download cmake and follow installation guide (https://cmake.org/install/)
- Create a bin/ directory in the same location as CMakeLists.txt
- Go to bin/ and run cmake .. in order to perform an out-of-source build
- All produced files and folders  will be placed at bin/
