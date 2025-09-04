# ROH-pipeline

**ROH-pipeline** is a Snakemake workflow for detecting **Runs of Homozygosity (ROH)**, estimating heterozygosity, and generating summary statistics and plots from short-read sequencing data.  
It integrates standard preprocessing (QC, trimming, alignment), variant calling, and downstream ROH/heterozygosity analyses.

---

## Features
- Quality control of raw reads with **FastQC**
- Adapter and quality trimming with **Trimmomatic**
- Alignment to a reference genome with **BWA**
- BAM sorting and indexing with **Samtools**
- Variant calling (planned: **GATK HaplotypeCaller**)
- Frequency-based genome statistics
- Estimation of heterozygosity with **ANGSD**
- ROH estimation with **PLINK**
- Generation of summary tables and plots:
  - ROH length distributions
  - ROH SNP distributions
  - Relationships between heterozygosity and ROH statistics (SROH, NROH, FROH)

---

## Dependencies

The pipeline relies on [Snakemake wrappers](https://snakemake-wrappers.readthedocs.io), as well as external tools.  
You should have the following installed and available in your environment:

- [Snakemake](https://snakemake.readthedocs.io) (≥7.0 recommended)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [BWA](http://bio-bwa.sourceforge.net/)
- [Samtools](http://www.htslib.org/)
- [GATK](https://gatk.broadinstitute.org/) (optional, for variant calling)
- [PLINK](https://www.cog-genomics.org/plink/)
- [ANGSD](http://www.popgen.dk/angsd/)
- **Python ≥3.8** with:
  - pandas
  - matplotlib / seaborn (for plotting scripts)

---

## Input files

1. **Reference genome**  
   Defined in `config/config.yaml`, e.g.:

   ```yaml
   reference: "genome.fasta"
