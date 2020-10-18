# LoMeX

LoMeX is a long k-mer extraction program that utilizes spaced seeds. This program consists of two separate programs: a modified version of an existing k-mer extraction program Squeakr (https://github.com/splatlab/squeakr), and the actual LoMeX program itself.

## Requirements
- SpacedSeedSqueakr has the same dependencies as the original Squeakr program (which can also be found here https://github.com/splatlab/squeakr).
  - libboost-dev 1.58.0.1ubuntu1
  - libssl-dev 1.0.2g-1ubuntu4.6
  - zlib1g-dev 1:1.2.8.dfsg-2ubuntu4
  - bzip2 1.0.6-8

- C++11 compiler
- OpenMP (version X)

## Installation

1. Download the project directory.
2. Install modified Squeakr.
- Go to ```SpacedSeedSqueakr``` directory and run command ```make squeakr```
3. Install LoMeX.
- Go to ```LoMeX``` directory and run command ```make```

## How to use

1. Run SpacedSeedSqueakr exact k-mer count with

```
./SpacedSeedSqueakr/squeakr count -e -c $M -k $K -s $S -t $T -p $SEED -o $OUTPUT $READS
  
  Required parameters:
  -e : use exact counting (must be used)
  -c : M = minimum abundance
  -k : K = number of fixed positions
  -s : S = lognumslots (refer to the instructions of the original Squeakr)
  -p : SEED = spaced seed pattern e.g. "5-10-5-10-5".
  -o : OUTPUT = output file path
  READS = path to the read file
  
  Optional parameters:
  -t : T = number of threads
```
  Spaced seed requirements:

  Blocks of fixed and "don't care" characters are separated with "-". The number is equal to the block length. Blocks in odd positions are fixed, even positions have the "don't care" blocks. Total length of the spaced seed must be odd and less than. Total length of the fixed blocks must be less than 32. Block length must be at least 1. Spaced seed can have any number of blocks. Middle character of the spaced seed must be in a fixed block. Spaced seed must be symmetrical.


2. Produce readable Squeakr output with

```
./SpacedSeedSqueakr/squeakr list -f $OUTPUT -o $COUNTS

  Required paramters:
  -f : OUTPUT = path to the output of the previous step
  -o : COUNTS = path to the file with extracted k-mers and their counts
```

3. Run LoMeX with
```
./LoMeX/lomex-main.out -k $COUNTS -r $READS -s $SEED -o $OUTPUT_PATH -q $OUTPUT_NAME -w $TMP -m $M -n $N -l $L -b $B -i $I -t $T
  
  Required parameters:
  -k : COUNTS = path to the output of the previous step
  -r : READS = path to the read file
  -s : S = spaced seed pattern
  -o : OUTPUT_PATH = path to a directory where output is created
  -q : OUTPUT_NAME = name of the output file
  -w : TMP = path to an empty directory where temporary files are created (must exist and be empty, inconvenient -> fix this later)
  
  Optional parameters:
  -m : M = minumum spaced k-mer abundance (default=2)
  -n : N = minimum solid character abundance (default=2)
  -l : L = relative minimum solid character abundance (default=0.1)
  -b : B = number of spaced k-mers occurrences that can be stored in the buffer before disk writing (default=5000000)
  -i : I = number of iterations (default=1)
  -t : T = number of threads for search step (default=1) 
```








