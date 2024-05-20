# Genetic relatedness of honeybees using SIMplyBee

## Introduction 
This repository contain R scripts that demonstrates the principles of expected and realised genetic relatedness among individual honeybees using [SIMplyBee](https://cran.r-project.org/web/packages/SIMplyBee/index.html). 

If you have no experience with SIMplyBee, we suggest you go to the software's website [simplybee.info](http://www.simplybee.info).

## Respository contents
- ```README.md``` is this file
- ```SimulateFounderPop.R``` demonstrates how SIMplyBee functions are used to create founder population genomes. Since this process can be timely, we have also provided the output file ```FounderGenomes_ThreePop_16chr.RData``` which can be downloaded and used to explore our scripts.
- ```YearCycle_Relatedness.R``` is the script of a 10 year simulation cycle, where genomic relationships are calculated for three different honeybee populations
- ```Plot_YearCycle_output.R``` demonstrates how to process the output data of the ```YearlyCycle_Relatedness.R``` and create plots to visually view the data.
