# Q-rator

Q-rator is an R lang. script that automates creation of files for FlexQTL and other genetic datasets, written for R version 4.0.5+. It requires several input files (which must be converted to comma-separated value files before running) and creates several output files, which we will enumerate below.  

## Required Files

**Map file:** intmap11_20k 
**Marker data:** 20k_8koverlap 
**Phenotype file:** The formatting of your phenotype data must follow these rules. List all individuals in column A with column name “Index”. Data for each individual goes in columns B and onward. Each column of phenotype data must have a descriptive column name in Row 1.

