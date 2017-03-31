Hieralign- Symmetric Word Alignment using Hierarchial Sub-sentential Alignment Method with Various IBM1 Models
==========

`Hieralign` is yet another symmetric word alignment tool simple, fast and unsupervised, 
it has been implemented with multi-threading in all stages, which makes it run a little faster than fast_align. 

The source code in this repository is provided under the terms of the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0.html).

## Characters:
1. Pure C++ implementation, 
   included the new features of C++11,  no additional library required.

2. multi-threading,
    It has been implemented with multi-threading in all stages (same as  fast_align, using OpenMP) and well-optimized.

3. cross-platform CMake file, easy and one-step compiling and installation.

4. replaced the best-1 version of hierarchical subsentential alignment with beam search version.

5. a little faster than fast_align, less than 2 minutes for 300,000 lines.

6. well-formed codes and simple class structure.

## Input format

The input to `Hieralign` must be tokenized and aligned into parallel sentences. Each line is a source language sentence and its target language translation, separated by a triple pipe symbol with leading and trailing white space (`\t`). An example 2-sentence English-Japanese parallel corpus is:

 
      I have a book . 私 は 本 を 持って いた 。
      I do not know the answer    私 は その 答え を 知ら ない

## Compiling and Installation

Building `Hieralign` requires a modern C++ compiler and the [CMake]() build system. Additionally, the OpenMP library is required here to obtain better performance.
Please make sure your complier version (gcc/clang) supports OpenMP before compiling. 

To compile, do the following 

    cmake . 
    make
    
## Usage
Run `Hieralign -h` to see a list of command line options.
There are 3 ways to use `Hieralign`: 

1. as a normal word aligner:
`Hieralign` generates symmetric alignments (source–target, i.e., left language–right language) alignments is:
 
    ./Hieralign --train train.txt --model vbibm1 > aligned.hier

2. train the model in advance:
The usually and recommended way to train and store a model is to just add the `--model_file` option:

    ./Hieralign --train train.txt --model vbibm1 --model_file path_to_store_model

3. as an online aligner:
Given the trained model, `Hieralign` can be ran as an online aligner:

   ./Hieralign --input test.txt --model vbibm1 --model_file path_to_store_model > aligned.hier

## Output

`Hieralign` produces outputs in the widely-used `i-j` “Pharaoh format,” where a pair `i-j` indicates that the <i>i</i>th word (zero-indexed) of the left language (source) is aligned to the <i>j</i>th word of the right sentence (target) but have been symmetrized. For example, a good alignment for first sentence would be:

      0-0 0-1 1-4 1-5 3-2 4-6
      0-0 1-6 2-6 3-5 4-2 5-3

