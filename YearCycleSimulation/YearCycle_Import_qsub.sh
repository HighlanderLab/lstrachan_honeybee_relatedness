#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N YearCycle
#$ -cwd       
#$ -l h_vmem=40G
#$ -pe sharedmem 1
#$ -o YearCycle.txt
#$ -l h_rt=20:00:00
#$ -j yes
#$ -P roslin_HighlanderLab


# Initialise the environment modules
. /etc/profile.d/modules.sh
module load roslin/R/4.1.0
 
# Run the program
Rscript YearCycle_SpringerRelationships_Import.R &
./cpumemlog.sh $!
