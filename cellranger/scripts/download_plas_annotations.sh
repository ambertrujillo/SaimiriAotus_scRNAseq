#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Download Plasmodium vivax annotations
# ----------------------------------------------------------------------------------------

# --- Download annotations for Plasmodium vivax genome, vivax (ASM241v2)


cd genomes/

PLASMODIUM_URL=ftp.ensemblgenomes.org/pub/protists/release-54/gtf

ANNO1=plasmodium_vivax/Plasmodium_vivax.ASM241v2.54.gtf

for ANNO in $ANNO1; do
    wget $PLASMODIUM_URL/$ANNO.gz \
        -O `basename $ANNO`.gz
    gunzip `basename $ANNO`.gz
done

mv Plasmodium_vivax.ASM241v2.54.gtf vivax.gtf #rename gtf file

cd ..
