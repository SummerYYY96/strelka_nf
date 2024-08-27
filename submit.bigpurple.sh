#!/bin/bash

#SBATCH -o slurm-%j.out
#SBATCH -J test_strelka
#SBATCH -p intellispace
#SBATCH --time=5-00:00:00
#SBATCH --ntasks-per-node=1
#SBATCH -c 8
#SBATCH --mem 48G
#SBATCH --export=HOSTNAME

module load singularity/3.9.8
./nextflow run main.nf
