#!/bin/bash
# Hao WANG

############################################################################
# train a ngram-based SMT model using MAIRE v1.1
############################################################################
# Written by Hao WANG, 2016. 11. 13
 
 
MOSES_DIR=/itigo/files/Tools/TranslationEngines/mosesdecoder-RELEASE-2.1.1 
MOSES_SCRIPT_DIR=${MOSES_DIR}/scripts
MOSES_BIN_DIR=${MOSES_DIR}/bin 
 
############################################################################
##########################  finish_config  #################################
############################################################################


WORK_DIR=model
mkdir -p ${WORK_DIR}

${MOSES_SCRIPT_DIR}/training/train-model.perl \
    -cores 4 \
    --first-step 4 \
    --last-step 4 \
    --corpus data/encode.train \
    --f source \
    --e target \
    --parallel \
    --alignment grow-diag-final-and
#--score-options "--GoodTuring"
