#!/bin/bash

# make directory
mkdir genofreq

# list all samples:
awk '{print $3}' config/table_reads.tsv | tail -n +2 > genofreq/list_samples.tsv
# list populations:
awk '{print $2}' config/table_reads.tsv | tail -n +2 | sort | uniq  > genofreq/list_populations.tsv
# table samples:
awk '{print $2 "\t" $3}' config/table_reads.tsv  | tail -n +2 | sort > genofreq/list_specimen_population.tsv


# Make the populations VCF files:
while read f; 
do 
	grep -w "$f" genofreq/list_specimen_population.tsv | cut -f2 > genofreq/specimens.tsv  ;
	bcftools view -S genofreq/specimens.tsv calls/all.vcf -O z -o genofreq/population_"$f".vcf.gz 
done < genofreq/list_populations.tsv



while read f;
do
	# print chrom pos format fields
	# Solution for genotypes heres:
	# https://www.biostars.org/p/312304/
	bcftools query -f '%CHROM %POS  %REF  %ALT [ %GT]\n' genofreq/population_"$f".vcf.gz > genofreq/table_"$f".tsv ;



	# We will attribute a "0" to homozygous loci, and "1" to heterozygous loci.
	# sum("0")/ sum("1" & "0") will be the rate of homozygous in the group. 

	# Expected frequencies:
	# P 100% A/A -> freq 0
	# F1  100%  A/a -> freq 1
	# BC genotypes output AA aA AA aA = 50% aA 50% AA
	# 50% heterozygous  50% homozygous -> freq 0.5

	# replace heterozygous by 1 and homozygous by 0:
	sed 's%\./\.%NA%g; s%\.[|]\.%NA%g; s%0/0%0%g; s%0[|]0%0%g; s%0/1%1%g; s%0[|]1%1%g; s%1/1%0%g; s%1[|]1%0%g; s%1/2%1%g; s%1[|]2%1%g; s%0/2%1%g; s%0[|]2%1%g; s%2/2%0%g; s%2[|]2%0%g; s%0/3%1%g; s%0[|]3%1%g; s%3/3%0%g; s%3[|]3%0%g; s%1/3%1%g; s%1[|]3%1%g ; s%2/3%1%g; s%2[|]3%1%g; s%0/4%1%g; s%0[|]4%1%g; s%1/4%1%g; s%1[|]4%1%g; s%2/4%1%g; s%2[|]4%1%g; s%3/4%1%g; s%3[|]4%1%g; s%4/4%0%g; s%4[|]4%0%g' genofreq/table_"$f".tsv > genofreq/table_geno_"$f".tsv ;

	# make the headers of the matrix for pandas processing: 
	#	zgrep "#CHROM" ./"$f"/GVCF_SNPs.vcf.gz | awk '{$1=$2=$3=$4=$5=$6=$7=$8=$9=""; print $0}' | sed 's/.conch_baits//g' | sed 's/ /" "/g'  | sed 's/ "/, "/g'  > headers_"$f".tsv ;



	#zgrep "#CHROM" "$f".vcf.gz  | sed 's/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t//' | sed 's/\t/\n/'  > samples_list_"$f".tsv ;

	zgrep "#CHROM" genofreq/population_"$f".vcf.gz   | sed 's/\tID//g' | sed 's/QUAL\tFILTER\tINFO\tFORMAT\t//' > genofreq/samples_list_"$f".tsv ;

	cat genofreq/samples_list_"$f".tsv > genofreq/geno_freq_"$f".tsv ;
	cat genofreq/table_geno_"$f".tsv >> genofreq/geno_freq_"$f".tsv ;
	sed 's/\t/,/g' genofreq/geno_freq_"$f".tsv | sed 's/ /,/g' | sed 's/,,/,/g' > genofreq/geno_freq_"$f".csv ;

done < genofreq/list_populations.tsv

