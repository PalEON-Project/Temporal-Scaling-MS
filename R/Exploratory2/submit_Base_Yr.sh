#!/bin/sh
#$ -wd /projectnb/dietzelab/paleon/Temporal-Scaling-MS/R/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -m e
#$ -M crollinson@gmail.com
#$ -l h_rt=120:00:00
#$ -q 'geo*'
#$ -pe omp 10
#$ -v OMP_NUM_THREADS=10
#$ -N Base_Yr
#cd /projectnb/dietzelab/paleon/Temporal-Scaling-MS/R/
R CMD BATCH 2_process_drivers_all_drivers_Base_Yr.R 
