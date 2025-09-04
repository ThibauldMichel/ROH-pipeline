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

The pipeline relies on [Snakemake wrappers](https://snakemake-wrappers.readthedocs.io), as well as external tools.You should have the following installed and available in your environment:

- [Snakemake](https://snakemake.readthedocs.io) (≥7.0 recommended)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [BWA](http://bio-bwa.sourceforge.net/)
- [Samtools](http://www.htslib.org/)
- [GATK](https://gatk.broadinstitute.org/)
- [PLINK](https://www.cog-genomics.org/plink/)
- [ANGSD](http://www.popgen.dk/angsd/)
- **Python ≥3.8** with:
  - pandas
  - matplotlib / seaborn (for plotting scripts)

---

## Input files

1. **Reference genome**Defined in `config/config.yaml`, e.g.:

   ```yaml
   reference: "genome.fasta"

   ```
2. **Sample table**
   Located at `config/table_reads.tsv`, tab-delimited with at least the following columns:

   - `read_ID` : unique identifier for each sample
   - `Population` : population assignment

The pipeline assumes raw reads are in `reads/{sample}_1.fq.gz` and `reads/{sample}_2.fq.gz`.

Example `table_reads.tsv`:

```
read_ID    Population
sample1    PopA
sample2    PopB
sample3    PopA
```

---

## Usage

1. **Clone the repository** and move into the project directory:
   git clone <repo_url>
   cd ROH-pipeline


2. **Prepare config files** :

* Place your reference genome path in `config/config.yaml`
* Add sample metadata in `config/table_reads.tsv`
* Place raw reads in `reads/` as `sample_1.fq.gz`, `sample_2.fq.gz`

3. Run the pipeline with Snakemake:

   ```
   snakemake --cores 8
   ```
   Adjust `--cores` according to available resources.
4. Outputs will be written to:

   * `trimmed/` : cleaned reads
   * `mapped/` : aligned BAM files
   * `calls/` : variant calls
   * `heterozygosity/` : heterozygosity matrices
   * `genofreq/` : genome frequency data
   * `plink/` : PLINK ROH results
   * `plots/` : summary figures
   * `output/` : combined ROH statistics

   ## Example outputs


   * `plots/plot_ROH_length_hist.png` : histogram of ROH lengths
   * `plots/plot_ROH_SNPs_hist.png` : histogram of SNP counts per ROH
   * `plots/plot_SROH_NROH_FROH.png` : summary of ROH statistics across samples
   * `plots/plot_heterozygosity_FROH_regline.png` : regression of heterozygosity vs FROH









