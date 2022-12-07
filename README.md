# SaimiriAotus_scRNAseq
scRNAseq pipeline for experimental infection of Aotus and Saimiri with P. vivax
For additional notes, look at [this](https://bioinformaticsworkbook.org/dataAnalysis/RNA-Seq/Single_Cell_RNAseq/Chromium_Cell_Ranger.html#gsc.tab=0) link.


# Pipeline

```bash
mkdir genomes
cd genomes
```

## -> Cellranger
Cell ranger is for mapping sequenced reads to concatenated reference genomes of host and pathogen.
###  Prepare GTF, genome, and reads

**1. Download annotation and reference genome files** 
   + NOTE: *Cell ranger requires a gtf file, not a gff file* 
> Necessary module(s): kent/385, samtools/intel/1.11

  * _Aotus nancymaae_

```bash
scripts/cellranger/download_aotus.sh
```
```bash
scripts/cellranger/download_aotus_annotations.sh
```

  * _Saimiri boliviensis_
  
```bash
scripts/cellranger/download_saimiri.sh
```
```bash
scripts/cellranger/download_saimiri_annotations.sh
```

  * _Plasmodium vivax_
  
```bash
scripts/cellranger/download_plas.sh
```
```bash
scripts/cellranger/download_plas_annotations.sh
```

**2. Prepare each reference genome and annotation file for the creation of a concatenated reference genome and annotation file.** 
> Necessary module(s): cellranger/7.0.0, bcl2fastq/2.20.0.422

  * _Aotus nancymaae_
 
 ```bash
cellranger mkgtf aotus.gtf aotus.filtered.gtf \
    --attribute=gene_biotype:protein_coding \
    --attribute=gene_biotype:lincRNA \
    --attribute=gene_biotype:antisense \
    --attribute=gene_biotype:IG_LV_gene \
    --attribute=gene_biotype:IG_V_gene \
    --attribute=gene_biotype:IG_V_pseudogene \
    --attribute=gene_biotype:IG_D_gene \
    --attribute=gene_biotype:IG_J_gene \
    --attribute=gene_biotype:IG_J_pseudogene \
    --attribute=gene_biotype:IG_C_gene \
    --attribute=gene_biotype:IG_C_pseudogene \
    --attribute=gene_biotype:TR_V_gene \
    --attribute=gene_biotype:TR_V_pseudogene \
    --attribute=gene_biotype:TR_D_gene \
    --attribute=gene_biotype:TR_J_gene \
    --attribute=gene_biotype:TR_J_pseudogene \
    --attribute=gene_biotype:TR_C_gene
```

  * _Saimiri boliviensis_
  
```bash
cellranger mkgtf saimiri.gtf saimiri.filtered.gtf \
    --attribute=gene_biotype:protein_coding \
    --attribute=gene_biotype:lincRNA \
    --attribute=gene_biotype:antisense \
    --attribute=gene_biotype:IG_LV_gene \
    --attribute=gene_biotype:IG_V_gene \
    --attribute=gene_biotype:IG_V_pseudogene \
    --attribute=gene_biotype:IG_D_gene \
    --attribute=gene_biotype:IG_J_gene \
    --attribute=gene_biotype:IG_J_pseudogene \
    --attribute=gene_biotype:IG_C_gene \
    --attribute=gene_biotype:IG_C_pseudogene \
    --attribute=gene_biotype:TR_V_gene \
    --attribute=gene_biotype:TR_V_pseudogene \
    --attribute=gene_biotype:TR_D_gene \
    --attribute=gene_biotype:TR_J_gene \
    --attribute=gene_biotype:TR_J_pseudogene \
    --attribute=gene_biotype:TR_C_gene
```

  * _Plasmodium vivax_

```bash
cellranger mkgtf vivax.gtf vivax.filtered.gtf \
    --attribute=gene_biotype:protein_coding \
    --attribute=gene_biotype:lincRNA \
    --attribute=gene_biotype:antisense \
    --attribute=gene_biotype:IG_LV_gene \
    --attribute=gene_biotype:IG_V_gene \
    --attribute=gene_biotype:IG_V_pseudogene \
    --attribute=gene_biotype:IG_D_gene \
    --attribute=gene_biotype:IG_J_gene \
    --attribute=gene_biotype:IG_J_pseudogene \
    --attribute=gene_biotype:IG_C_gene \
    --attribute=gene_biotype:IG_C_pseudogene \
    --attribute=gene_biotype:TR_V_gene \
    --attribute=gene_biotype:TR_V_pseudogene \
    --attribute=gene_biotype:TR_D_gene \
    --attribute=gene_biotype:TR_J_gene \
    --attribute=gene_biotype:TR_J_pseudogene \
    --attribute=gene_biotype:TR_C_gene
```

**3. Make concatenated host-pathogen reference genome. (SBATCH JOB)** 

```bash
sbatch cellranger/prepare_ref_genome.sbatch
```

**4. Download FASTQ reads**
+ NOTE: *Hand make list of SRRs as SRR.numbers by pasting SRRs into text file called "SRRnumbers.txt"*

```bash
mkdir data/

FILENAME="SRRnumbers.txt"

LINES=$(cat $FILENAME)

for LINE in $LINES
do
    echo "$LINE"
    fasterq-dump --split-files $LINE
done


mv data/SRR*/*.sra ..
rm -r data/SRR*
mv *.sra data/


tar -zcvf singlecell_fastqs.tar.gz data/
```
Seperate Aotus and Saimiri reads into their own directories:

```bash
mkdir data/saimiri
mkdir data/aotus

#saimiri
mv data/SRR11008269* data/saimiri/
mv data/SRR11008272* data/saimiri/
mv data/SRR11008273* data/saimiri/
mv data/SRR11008276* data/saimiri/

#aotus
mv data/SRR* data/aotus
```
Rename FASTQ files for cellranger:

```bash
cd data/saimiri

mv SRR11008269_1.fastq SRR11008269_S1_R1_001.fastq 
mv SRR11008269_2.fastq SRR11008269_S1_R2_001.fastq
mv SRR11008272_1.fastq SRR11008272_S2_R1_001.fastq
mv SRR11008272_2.fastq SRR11008272_S2_R2_001.fastq
mv SRR11008273_1.fastq SRR11008273_S3_R1_001.fastq
mv SRR11008273_2.fastq SRR11008273_S3_R2_001.fastq
mv SRR11008276_1.fastq SRR11008276_S4_R1_001.fastq
mv SRR11008276_2.fastq SRR11008276_S4_R2_001.fastq

cd data/aotus

mv SRR11008270_1.fastq SRR11008270_S1_R1_001.fastq
mv SRR11008270_2.fastq SRR11008270_S1_R2_001.fastq
mv SRR11008271_1.fastq SRR11008271_S2_R1_001.fastq
mv SRR11008271_2.fastq SRR11008271_S2_R2_001.fastq
mv SRR11008274_1.fastq SRR11008274_S3_R1_001.fastq
mv SRR11008274_2.fastq SRR11008274_S3_R2_001.fastq
mv SRR11008275_1.fastq SRR11008275_S4_R1_001.fastq
mv SRR11008275_2.fastq SRR11008275_S4_R2_001.fastq
mv SRR11008277_1.fastq SRR11008277_S5_R1_001.fastq
mv SRR11008277_2.fastq SRR11008277_S5_R2_001.fastq
mv SRR11008278_1.fastq SRR11008278_S6_R1_001.fastq
mv SRR11008278_2.fastq SRR11008278_S6_R2_001.fastq
```
+ NOTE: *In order to parallelize it (array job) need to make separate folders for each individual*

```bash
cd data/saimiri

mkdir SRR11008269_S1_R1
mv SRR11008269_S* SRR11008269_S1_R1/
mkdir SRR11008272_S2_R1
mv SRR11008272_S* SRR11008272_S2_R1/
mkdir SRR11008273_S3_R1
mv SRR11008273_S* SRR11008273_S3_R1/
mkdir SRR11008276_S4_R1
mv *.fastq SRR11008276_S4_R1

cd data/aotus
mkdir SRR11008270_S1_R1
mv SRR11008270_S* SRR11008270_S1_R1/


mkdir SRR11008271_S2_R1
mv SRR11008271_S* SRR11008271_S2_R1/
mkdir SRR11008274_S3_R1
mv SRR11008274_S* SRR11008274_S3_R1/
mkdir SRR11008275_S4_R1
mv SRR11008275_S* SRR11008275_S4_R1/
mkdir SRR11008277_S5_R1
mv SRR11008277_S* SRR11008277_S5_R1/
mkdir SRR11008278_S6_R1
mv SRR11008278_S* SRR11008278_S6_R1/
```

+ NOTE: *Each individual has it's own cellranger_count.sbatch in scripts/aotus  (or saimiri)*
Run the following for each individual

Example:
```bash
sbatch cellranger/scripts/aotus/cellranger_count_70.sbatch
```


