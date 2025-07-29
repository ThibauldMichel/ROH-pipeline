#!/bin/bash

#SBATCH --job-name="impreads"
#SBATCH --export=ALL



rsync -azvh /home/tmichel/projects/rbge/Begonia_genomes/Hillebrandia_illumina/* ./reads/
rsync -azvh /home/tmichel/projects/rbge/JiaxChen/Bsub_Bflu_data/skims/*.fq.gz  ./reads/
