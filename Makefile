# Makefile to run the pipeline
SHELL:=/bin/bash
export NXF_VER:=19.07.0

install:
	curl -fsSL get.nextflow.io | bash

submit:
	sbatch submit.bigpurple.sh
