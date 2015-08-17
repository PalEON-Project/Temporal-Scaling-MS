#!/bin/sh
#$ -wd /projectnb/dietzelab/paleon/Temporal-Scaling-MS/R/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -m e
#$ -M crollinson@gmail.com
#$ -l h_rt=120:00:00
#$ -q 'geo*'
#$ -N GSbyRes
#cd /projectnb/dietzelab/paleon/Temporal-Scaling-MS/R/
R CMD BATCH 2a_process_drivers_all_drivers_byResolution_GS.R 
