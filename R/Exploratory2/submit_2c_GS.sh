#!/bin/sh
#$ -wd /projectnb/dietzelab/paleon/Temporal-Scaling-MS/R/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -m e
#$ -M crollinson@gmail.com
#$ -l h_rt=120:00:00
#$ -q 'geo*'
#$ -N GSbySite
#cd /projectnb/dietzelab/paleon/Temporal-Scaling-MS/R/
R CMD BATCH 2c_process_drivers_all_drivers_bySite_GS.R 
