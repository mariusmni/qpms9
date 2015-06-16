# qPMS9

This project includes an efficient parallel algorithm for quorum Planted Motif Search (qPMS). Planted Motif Search (PMS) receives as input n biological sequences and two integers L and D. It returns all possible sequences M of length L such that M occurs in each of the input strings, and each occurrence of M differs from M in no more than D positions. A more practically usable formulation is called quorum PMS (qPMS), which allows the motif to appear in at least q% of the input strings. Therefore, PMS is a special case of qPMS, where q=100%.

This algorithm is described in the following publication

Marius Nicolae & Sanguthevar Rajasekaran, "qPMS9: An Efficient Algorithm for Quorum Planted Motif Search"

http://www.nature.com/srep/2015/150115/srep07813/full/srep07813.html

## Source code

See the [qpms9](qpms9) folder for qPMS9 source code and details on how to compile/run.

See the [qpms9-data](qpms9-data) folder for dataset generator source code. 

## Precompiled binaries

Copyright Marius Nicolae 
marius.nicolae (at) engr.uconn.edu 
