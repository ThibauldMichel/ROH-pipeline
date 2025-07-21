# ROH-pipeline

The pipeline is for the moment split between different BASH and Python scripts. A Snakemake wrapper is in preparation.

Put the script in the same directory than the common VCF file containing all the samples.

1. Run script_plink.sh
2. Run script_geno_freq.sh 
3. Run calculation_minimal_number_SNP_for_ROH_loop.py
4. Finally run hom_SROH_NROH.py on the outputs.
