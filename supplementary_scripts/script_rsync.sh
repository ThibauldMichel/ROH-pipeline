#!/bin/bash

#SBATCH --job-name="rsync"
#SBATCH --export=ALL




rsync -azvh /home/tmichel/scratch/ROH-pipeline_PNG/* ./
