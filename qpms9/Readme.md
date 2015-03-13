# qPMS9 - Parallel Algorithm for Quorum Planted Motif Search

Copyright Marius Nicolae 
marius.nicolae (at) engr.uconn.edu 

## Prerequisites: 
- unix or cygwin environment
- openmpi library (mpi and mpi_cxx libraries are required for parallel execution)
- g++
- make

## To compile 
For the parallel version (recommended) type: make
The executable will be found in the Release folder.

To compile without mpi (only for single core execution) type: make nompi
The executable will be found in the NoMpi folder.

## To clean 
For the parallel version: make clean
To clean version without mpi: make clean-nompi

## To run
Release/qPMS9 or NoMpi/qPMS9 - will print the following:

Arguments: inputFile -l <motifLength> -d <#Mutations> [-q <quorumPercent>] [-s <stackSz>] [-n nPrime] [-o <outputFile>] [-h]

Where:

  -l: length of motif

  -d: max changes for planted instances of the motif

  -o: outputFile [stdout if not specified]

  -q: min percentage of strings that have motif [default 100]

  -s: (optional) number of lmers for which we generate common neighborhoods

  -n: (optional) number of sequences for which we initially find motifs

  -h: print this help
