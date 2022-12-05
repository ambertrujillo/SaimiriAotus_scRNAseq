#!/bin/sh

# Script to download Plasmodium vivax genome, vivax (ASM241v2)

GENOME_FA=vivax.fa

wget 'ftp.ensemblgenomes.org/pub/protists/release-54/fasta/plasmodium_vivax/dna/Plasmodium_vivax.ASM241v2.dna.toplevel.fa.gz' \
    -O ${GENOME_FA}.gz

gunzip -c ${GENOME_FA}.gz > $GENOME_FA  

LAST_OK_LINE=$((`grep -n "^>[^c]" $GENOME_FA | head -n 1 | cut -d":" -f 1` - 1))
if [ $LAST_OK_LINE -gt 0 ]; then
    mv $GENOME_FA ${GENOME_FA}.backup
    head -n $LAST_OK_LINE ${GENOME_FA}.backup > ${GENOME_FA}
    rm ${GENOME_FA}.backup   
fi

mkdir tmp_for_sort
faSplit byname ${GENOME_FA} tmp_for_sort/
cd tmp_for_sort/;
ls -v | xargs cat > ../${GENOME_FA}
cd ..
rm -r tmp_for_sort

exit
