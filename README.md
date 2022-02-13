# File organization
- configure.sh: script for changing file permissions and compiling the program
- demo: files for demo
    - demo.sh: script for runing the demo
    - input.pqr: an example input file of protein

- package: PES models and related script
    - ALA et al.: scripts for calling the models
    - model: PES models
    - lib: DeepMD lib
    - bin: script for trnafer input file

- pesmpi.sh: script for running the program
- README.md: document
- src: source code
    - makefile: makefile for compiling the program
    - fragment.f90: souce code for cutting proteins into fragments.
    - runfrag-pesmpi.f90: souce code for calling DeepMD and PES models to calculating fragments.
    - cal_energy_pes.f90: souce code for calculating the final energy and atomic forces with the output files of runfrag-pesmpi.f90

# Environment configuration
## System
The scripts for running this program were written by Bash shell.
All the tests were done on Centos 8.6 system.

## Compiler and mpi library
The program was written by Fortran 90 and requre fortran compiler and mpi library.
In general, the code can be compiled by any fortran compiler.
In our test, we used ifort and mpif90 to compile the program and used mpirun to run the program inparallel. The ifort, mpif90, and mpirun were all from intel ONEAPI.
The simplest way for preparing compiler and mpilibray is to install intel ONEAPI (free software).

## DeepMD
We used DeepMD (version 2.0.0) to train the pes models.
Download and install DeepMD version 2.0.0 or later. (Python 3.x and conda are provided by DeepMD)

# Building program
Run the script `configure.sh` by `./configure.sh`.

# Instructions for use
After environment configuration and building program, we can use this program to calculate the energy and atomic forces of proteins in pqr format.
The program read the protein structure from one input file (input.pqr).
After preparing the `input.pqr` file, run the program by `./run.sh`
The energy is located in `energy.dat` (Hartree), force information is in `force.dat` (Hartree/Bohr).

# Demo
Run the demon:

```
cd demo
./demo
cat energy.dat
```

We should fine the following output:

```
 Total energy =     -6179.4632034124
```
