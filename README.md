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
      + *Cell ranger requires a gtf file, not a gff file* 
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

**2. Prepare each reference genome and annotation file for the creation of a concatenated reference genome and annotation file** 
> Necessary module(s): cellranger/7.0.0, bcl2fastq/2.20.0.422



