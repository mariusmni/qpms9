# Motif dataset generator

## To compile
```
make
```

## To run
```
Release/qpms9-data
```

```
Arguments: -o <outputFile> -n <nStrings> -m <stringLen> -l <motifLen> -d <motifMutations> -q <quorumPercent> -s <aSize> -a <alphabet> -r <randSeed>

Where:

  -o: outputFile [stdout if not specified]

  -n: number of strings to generate [default 20]

  -m: length of each string [default 600]

  -l: length of motif

  -d: max changes for planted instances of the motif

  -q: quorum percent - in what percentage of strings to plant motif [default 100]

  -s: alphabet size - 4 for DNA, 20 for protein, or specify arbitrary alphabet using -alphabet [default 4]

  -a: alphabet - specify arbitrary alphabet

  -r: seed for the random number generator; use the same seed for repeatable output
```
