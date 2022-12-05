#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Download Aotus annotations
# ----------------------------------------------------------------------------------------

# --- Download annotations for Aotus genome, Anon2.0
module load kent/385

cd genomes/ 

AOTUS_URL=ftp.ensembl.org/pub/release-107/gtf

ANNO1=aotus_nancymaae/Aotus_nancymaae.Anan_2.0.107.gtf

for ANNO in $ANNO1; do
    wget $AOTUS_URL/$ANNO.gz \
        -O `basename $ANNO`.gz
    gunzip `basename $ANNO`.gz
done

mv Aotus_nancymaae.Anan_2.0.107.gtf aotus.gtf #rename gtf file


cd ..
