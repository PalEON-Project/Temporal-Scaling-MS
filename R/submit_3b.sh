#!/bin/sh
#$ -wd /projectnb/dietzelab/paleon/Temporal-Scaling-MS/R/
#$ -j y
#$ -S /bin/bash
#$ -V
#$ -m e
#$ -M crollinson@gmail.com
#$ -l h_rt=168:00:00
#$ -N Interactions_4drivers
#cd /projectnb/dietzelab/paleon/Temporal-Scaling-MS/R/
R CMD BATCH 3b_process_interactions_4Drivers.R