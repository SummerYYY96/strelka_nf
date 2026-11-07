# Makefile to run the pipeline
SHELL:=/bin/bash
export NXF_VER:=25.10.9

install:
	curl -fsSL get.nextflow.io | bash

submit:
	sbatch submit.bigpurple.sh
