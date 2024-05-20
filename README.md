# Genetic relatedness of individual honeybees using SIMplyBee

## Introduction
This repository contain R scripts that demonstrates the principles of expected and realised genetic relatedness among individual honeybees using [SIMplyBee](https://cran.r-project.org/web/packages/SIMplyBee/index.html). 

If you have no experience with SIMplyBee, we suggest you go to the software's website [simplybee.info](http://www.simplybee.info) and study the tutorials. This repository uses SIMplyBee to simulate honeybee populations as described in the below manuscript and to calculate and summarise relatedness between individual honeybees.

TODO: Cite pre-print/manuscript

## Respository contents
- ```README.md``` is this file
- ```SimulateFounderPop.R``` uses SIMplyBee to create founder population genomes. Since this process can be time-consuming, we have also provided the output file ```FounderGenomes_ThreePop_16chr.RData``` that can be downloaded and used to by-pass this part of the simulation.
- ```YearCycle_Relatedness.R``` runs a yearly cycle of honeybee population for 10 years and calculates and summarizes genetic relatedness for three different honeybee populations.
- ```Plot_YearCycle_output.R``` processes the output data of the ```YearlyCycle_Relatedness.R``` and create plots to visualise the data.
