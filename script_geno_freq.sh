#!/bin/bash

VCF=subB.vcf.gz 

# Get names of the samples:
zgrep "#CHROM" "$VCF"  | sed 's/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t//' | sed 's/\t/\n/g'  > samples_list.tsv

# Make individual VCF files:
#while read f;
#do
#	bcftools view -s "$f" GVCF_SNPs.vcf.gz -O z -o "$f".vcf.gz;
#done < samples_list.tsv

# Make the populations VCF files:
while read f; 
do 
	#specimen=$(echo "$f" | cut -f1) ;
	#population=$(echo "$f" | cut -f2);
	grep -w "$f" list_specimen_population.tsv | cut -f1 > specimens.tsv  ;
	bcftools view -S specimens.tsv "$VCF" -O z -o "$f".vcf.gz 
done < list_populations.tsv



while read f;
do
	# print chrom pos format fields
	# Solution for genotypes heres:
	# https://www.biostars.org/p/312304/
	bcftools query -f '%CHROM %POS  %REF  %ALT [ %GT]\n' "$f".vcf.gz > table_"$f".tsv ;



	# We will attribute a "0" to homozygous loci, and "1" to heterozygous loci.
	# sum("0")/ sum("1" & "0") will be the rate of homozygous in the group. 

	# Expected frequencies:
	# P 100% A/A -> freq 0
	# F1  100%  A/a -> freq 1
	# BC genotypes output AA aA AA aA = 50% aA 50% AA
	# 50% heterozygous  50% homozygous -> freq 0.5

	# replace heterozygous by 1 and homozygous by 0:
	sed 's%\./\.%NA%g; s%\.[|]\.%NA%g; s%0/0%0%g; s%0[|]0%0%g; s%0/1%1%g; s%0[|]1%1%g; s%1/1%0%g; s%1[|]1%0%g; s%1/2%1%g; s%1[|]2%1%g; s%0/2%1%g; s%0[|]2%1%g; s%2/2%0%g; s%2[|]2%0%g; s%0/3%1%g; s%0[|]3%1%g; s%3/3%0%g; s%3[|]3%0%g; s%1/3%1%g; s%1[|]3%1%g ; s%2/3%1%g; s%2[|]3%1%g; s%0/4%1%g; s%0[|]4%1%g; s%1/4%1%g; s%1[|]4%1%g; s%2/4%1%g; s%2[|]4%1%g; s%3/4%1%g; s%3[|]4%1%g; s%4/4%0%g; s%4[|]4%0%g' table_"$f".tsv > table_geno_"$f".tsv ;

	# make the headers of the matrix for pandas processing: 
	#	zgrep "#CHROM" ./"$f"/GVCF_SNPs.vcf.gz | awk '{$1=$2=$3=$4=$5=$6=$7=$8=$9=""; print $0}' | sed 's/.conch_baits//g' | sed 's/ /" "/g'  | sed 's/ "/, "/g'  > headers_"$f".tsv ;



	#zgrep "#CHROM" "$f".vcf.gz  | sed 's/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t//' | sed 's/\t/\n/'  > samples_list_"$f".tsv ;

	zgrep "#CHROM" "$f".vcf.gz   | sed 's/\tID//g' | sed 's/QUAL\tFILTER\tINFO\tFORMAT\t//' > samples_list_"$f".tsv ;

	cat samples_list_"$f".tsv > geno_freq_"$f".tsv ;
	cat table_geno_"$f".tsv >> geno_freq_"$f".tsv ;
	sed 's/\t/,/g' geno_freq_"$f".tsv | sed 's/ /,/g' | sed 's/,,/,/g' > geno_freq_"$f".csv ;

done < list_populations.tsv

