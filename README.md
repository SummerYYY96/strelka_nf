## Motivation  
This exercise is to implement methods in [this publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC9308779/), as "Strelka2 (employing default parameters) with the blood sample as the tumoral input and the tumor sample as control (reverse calling)." The goal is to call clonal hematopoiesis in our targeted gene panels.

## How to run this script  
The script requires a list of samples ID, tumor bam (dedupped, realigned and recalibrated). The example process_samplesheet.py produces a samplesheet given a data directory with all bam files and a csv with sample tumor and normal IDs.  
The script then run nextflow processes by running the following on command line.  
```
make update # install nextflow under current directory
make submit # job submission to slurm
```
 
