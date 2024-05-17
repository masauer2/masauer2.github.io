layout: page
title: "GUIDE"
permalink: /docs
# NOTE ON 05.30.2023
# NEED TO TEST PUSHING PERMISSIONS
Here is how you can push changes to the main branch. 
First clone the repo with `git clone` and then follow the steps below.
(`git add` will add changes to the "stage" which are then committed. Pushing changes updates the repo)
```
git add .
git commit -m "Commit message"
git push -u origin main
```

# FRESEAN-CP2K
## Analysis of Anharmonic Vibrations in CP2K Simulations
## Required software

Required Software/Libraries
 - GNU make
 - C/c++ compiler (openmp enabled)
 - lfftw3f (http://www.fftw.org/download.html)
 - lgsl (https://www.gnu.org/software/gsl/)
 - lgslcblas (https://www.gnu.org/software/gsl/)

## Compiling Code
```
git clone https://github.com/masauer2/FRESEANCP2K.git
cd FRESEANCP2K
make
```

## Codebase tour
 - `FRESEAN` contains source code to perform the analysis
 - `example` contains sample data (stored in `example/data`), input files (stored in `example/input`) and output (stored in `example/output`)
 - `src` and `include` contains backend source code needed to compile code in `FRESEAN`
 - `bin` stores executables for analysis
 - `build` stores object files for analysis

# Running the code
A sample script is provided called `RUNME.sh` which will perform the analysis, which is as follows.
1. Generate a `.gro` reference structure for the system given CP2K position and velocity file (`bin/genREF.exe`).
2. Generate a `.mtop` custom topology file from CP2K position file (`bin/genMTOP.exe`).
3. Generate the frequency dependent cross correlation matrices (`bin/gen-modes_omp.exe`).
4. Generate the frequency dependent vibrational modes (`bin/eigen.exe`).
5. Extract the binary eigenvectors to `.xyz` format at a given frequency (`bin/extract.exe`).

## Code Usage
1. bin/genREF.exe <br>
- `./bin/genREF.exe input_position_file.xyz input_velocity_file.xyz output_reference_file.gro`
2.  bin/genMTOP.exe <br>
- `./bin/genMTOP.exe input_position_file.xyz output_topology.mtop`
3.  bin/gen-modes_omp.exe <br>
- `./bin/gen-modes_omp.exe input_file.inp`
- 
4.  bin/eigen.exe <br>
- `./bin/eigen.exe input_covar.mmat nCorr`
- nCorr is an integer representing length of the cross-correlation function. This should be the same as in the input file in the previous step.
5.  bin/extract.exe <br>
- `./bin/extract.exe input_file.inp`
- `extractMode` is used to select whether we want to extract a particlar frequency (mode 0) or a particular index (mode 1). For mode 0, an integer representing the desired extraction frequency in wavenumbers should be provided. The code will automatically determine the closest frequency available to be extracted. For example, if the frequency resolution is 4.2 cm^-1 and the `freqSel` variable is set to 8, the vibrational modes for 8.4 cm^-1 will be extracted (this correspons to index 3 if mode 1 is used). 

## Code Output

- `covar_CP2K.mmat` -> Cross correlation matrices in binary format
- `eval_covar_CP2K.mmat.dat` -> Eigenvalues for eigenvectors at each frequency (Each line is a different frequency)
- `evec_covar_CP2K.mmat` -> Eigenvectors for each frequency in binary format
- `evec_freqX_modeM-N_eigen.xyz` -> Human readable eigenvectors for frequency index X for modes between M and N
# masauer2.github.io
