# qPMS9 - Parallel Algorithm for Quorum Planted Motif Search

Copyright Marius Nicolae 
marius.nicolae (at) engr.uconn.edu 

## Prerequisites: 
- unix or cygwin environment
- openmpi library (mpiCC required for compiling parallel version)
- g++
- make

## To compile 
For the parallel version (recommended) type: 

```sh
make
```

The executable, named qpms9, will be found in the Release folder.

To compile without mpi (only for single core execution) type: 

```
make nompi
```

The executable will be found in the NoMpi folder.

To compile with debug info and without optimizations type:

```
make debug
```
The executable will be found in the Debug folder.

## To clean 
For the parallel version: 

```
make clean
```

To clean version without mpi: 

```
make clean-nompi
```

To clean debug version: 

```
make clean-debug
```

## To run
```
Release/qpms9
```
 
or 

```
NoMpi/qpms9
``` 

or

```
Debug/qpms9
``` 

It will print the following:

```
Arguments: inputFile -l <motifLength> -d <#Mutations> [-q <quorumPercent>] [-s <stackSz>] [-n nPrime] [-o <outputFile>] [-h]

Where:

  -l: length of motif

  -d: max changes for planted instances of the motif

  -o: outputFile [stdout if not specified]

  -q: min percentage of strings that have motif [default 100]

  -s: (optional) number of lmers for which we generate common neighborhoods

  -n: (optional) number of sequences for which we initially find motifs

  -h: print this help
```
