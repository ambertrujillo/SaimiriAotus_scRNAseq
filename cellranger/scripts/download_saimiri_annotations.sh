#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Download Saimiri boliviensis annotations
# ----------------------------------------------------------------------------------------

# --- Download annotations for Saimiri genome, SaiBol1.0

cd genomes/ 

SAIMIRI_URL=ftp.ensembl.org/pub/release-107/gtf

ANNO1=saimiri_boliviensis_boliviensis/Saimiri_boliviensis_boliviensis.SaiBol1.0.107.gtf

for ANNO in $ANNO1; do
    wget $SAIMIRI_URL/$ANNO.gz \
        -O `basename $ANNO`.gz
    gunzip `basename $ANNO`.gz
done

mv Saimiri_boliviensis_boliviensis.SaiBol1.0.107.gtf saimiri.gtf #rename gtf file

cd ..
