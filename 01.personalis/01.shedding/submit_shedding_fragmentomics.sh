#!/bin/bash

#$ -l tmem=16G
#$ -l h_vmem=16G
#$ -l h_rt=06:00:00
#$ -S /bin/bash
#$ -N fragmentomics
#$ -cwd
#$ -o /SAN/colcc/tracerx_personalis_pipeline/hannah/01.personalis/01.shedding/logs/
#$ -e /SAN/colcc/tracerx_personalis_pipeline/hannah/01.personalis/01.shedding/logs/

/SAN/colcc/tracerx_personalis_pipeline/hannah/py_env/bin/python -u /SAN/colcc/tracerx_personalis_pipeline/hannah/01.personalis/01.shedding/shedding_fragmentomics_analysis.py