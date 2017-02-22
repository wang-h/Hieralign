Hieralign- Word Alignment using Hierarchial Sub-sentential Alignment method with various IBM1 models
==========
Copyright @Hao Wang
email: oko_ips@ruri.waseda.jp

`Hieralign` is yet another word alignment tool simple, fast and unsupervised, it has been implemented with 
multi-threading in all stages, which makes it run faster than fast_align. 

The source code in this repository is provided under the terms of the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0.html).

## Characters:
1. Pure C++ implementation, 
   included the new features of C++11,  no additional library required.

2. multi-threading,
    It has been implemented with multi-threading in all stages (same as  fast_align, using OpenMP) and well-optimized.

3. cross-platform CMake file, easy and one-step compiling and installation.

4. replaced best-1 version of hierarchical subsentential alignment with beam search version.

5. a little faster than fast_align, less than 2 minutes for 300,000 lines.

6. well-formed codes and simple class structure.

## Input format

Input to `Hieralign` must be tokenized and aligned into parallel sentences. Each line is a source language sentence and its target language translation, separated by a triple pipe symbol with leading and trailing white space (`\t`). An example 2-sentence English-Japanese parallel corpus is:

 
   I think X is just a matter of time.	X は 時間 の 問題 と 思い ます	
   I think that X will become an issue in the future.	X は 今後 の 課題 と 思い ます

## Compiling and Installation

Building `Hieralign` requires a modern C++ compiler and the [CMake]() build system. Additionally, the OpenMP library is required here to obtain better performance. 

To compile, do the following 

    cmake . 
    make
    
## Usage
Run `Hieralign` to see a list of command line options.
There are 3 ways to use `Hieralign`: 

1. as a normal word aligner:
`Hieralign` generates symmetric alignments (source–target, i.e., left language–right language) alignments is:
 
    ./Hieralign --train train.txt --model ibm1_vbh > aligned.hier

2. train the model in advance:
The usually and recommended way to train and store a model is to just add the `--model_file` option:

    ./Hieralign --train train.txt --model ibm1_vbh --model_file path_to_store_model

3. as an online aligner:
Given the trained model, `Hieralign` can be ran as an online aligner:

   ./Hieralign --input test.txt --model ibm1_vbh --model_file path_to_store_model > aligned.hier

## Output

`Hieralign` produces outputs in the widely-used `i-j` “Pharaoh format,” where a pair `i-j` indicates that the <i>i</i>th word (zero-indexed) of the left language (source) is aligned to the <i>j</i>th word of the right sentence (target). For example, a good alignment for first sentence would be:

    0-0 1-1 2-2 2-3 3-4 4-5
