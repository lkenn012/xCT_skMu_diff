# xCT_skMu_diff
 Code and data for the analysis and plotting of C2C12 transcriptomics data as reported in Kanaan _et al_ (2023) {DOI}

## Code
Contains all Python and R files used in obtaining GEO transcriptomics datasets, analysis of the datasets, and plotting figures used in publication main text and supplementary

## Data
All files produced from the analyses of C2C12 code as well as ShinyGO outputs used in creating enrichment plots

## Environment information

The following version of Python and R along with all accompanying libraries are used in the code provided here. All code was implemented on Windows with an installation time of approximately 15 minutes. 

Python   --3.8.18

### Scientific programming & machine learning packages
- numpy   --1.24.3
- pandas   --2.0.3
- scikit-learn   --1.3.0
- scipy   --1.10.1
- statsmodels   --0.14.0

### Plotting packages
- matplotlib   --3.7.2
- seaborn   --0.12.2

Installation information:
```
conda install numpy
conda install pandas 
conda install scikit-learn   
conda install scipy
conda install statsmodels  
conda install matplotlib  
conda install seaborn  
```

R  --4.2.2
RStudio  --2022.12.0+353

### R libraries
- GEOquery_2.66.0
- readxl_1.4.2
- Biobase_2.58.0
- BiocGenerics_0.44.0

Installation information:
RStudio installation download from https://posit.co/download/rstudio-desktop/
```
BiocManager::install("GEOquery")
install.packages("readxl")
BiocManager::install("Biobase")
BiocManager::install("BiocGenerics")
```
