# Genetic relatedness of individual honeybees using SIMplyBee

## Introduction
This repository contain R scripts that demonstrates the principles of expected and realised genetic relatedness among individual honeybees using [SIMplyBee](https://cran.r-project.org/web/packages/SIMplyBee/index.html). 

If you have no experience with SIMplyBee, we suggest you go to the software's website [simplybee.info](http://www.simplybee.info) and study the tutorials. This repository uses SIMplyBee to simulate honeybee populations as described in the below manuscript and to calculate and summarise relatedness between individual honeybees.

As of the 20.05.24, we have a paper associated with these scripts in the final stages of writing. In the coming weeks, a link to the pre-print manuscript will be available here. 

## Respository contents
The `YearCycleSimulation` folder contains `Scripts` and `Plotting Data`



`Scripts` contains: 
* ```SimulateFounderPop.R``` uses SIMplyBee to create founder population genomes.
* ```YearCycle_Relatedness.R``` runs a yearly cycle of honeybee population for 10 years and calculates and summarizes genetic relatedness for three different honeybee populations.
* ```Plot_YearCycle_output.R``` processes the output data of the ```YearlyCycle_Relatedness.R``` and create plots to visualise the data.




`Plotting Data` contains: 
* `FounderGenomes_ThreePop_16chr.RData` is the output file of the script ```SimulateFounderPop.R``` that can be downloaded and used to by-pass this part of the simulation, as it can be time-consuming. 
