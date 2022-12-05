#!/bin/sh

# Script to download Saimiri genome, SaiBol1.0

GENOME_FA=saimiri.fa

wget 'ftp.ensembl.org/pub/release-107/fasta/saimiri_boliviensis_boliviensis/dna/Saimiri_boliviensis_boliviensis.SaiBol1.0.dna.toplevel.fa.gz' \
    -O ${GENOME_FA}.gz

gunzip -c ${GENOME_FA}.gz > $GENOME_FA  

LAST_OK_LINE=$((`grep -n "^>[^c]" $GENOME_FA | head -n 1 | cut -d":" -f 1` - 1))
if [ $LAST_OK_LINE -gt 0 ]; then
    mv $GENOME_FA ${GENOME_FA}.backup
    head -n $LAST_OK_LINE ${GENOME_FA}.backup > ${GENOME_FA}
    rm ${GENOME_FA}.backup   
fi

#Getting rid of things we don't want -- like alternative haplotypes
mkdir tmp_for_sort
faSplit byname ${GENOME_FA} tmp_for_sort/ 
cd tmp_for_sort/;
ls -v | xargs cat > ../${GENOME_FA}
cd ..
rm -r tmp_for_sort

exit
