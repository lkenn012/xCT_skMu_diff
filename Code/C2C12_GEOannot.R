### Code to annotate and extract gene expression data from GEO

## want to download GEO datasets by GSE ID
# BiocManager::install("GEOquery", force = TRUE)  # install GEOquery library

# need libraries
library(GEOquery) # for querying GEO database
library(readxl) # for loading xlsx files

# Use GEOquery to download the expression data
# download data for specified GSE ID
# GSE IDs for publish C2C12 transciptomics data over differentiation: 
#   GSE989, GSE84158, GSE4694, GSE46492, GSE16992, GSE148294, GSE126370, GSE11415, GSE110957, GSE108503
gse_id <- "GSE84158"    # Specify dataset of interest
gse <- getGEO(gse_id)

# Need to get information about time points in data (for each GSM in GSE)
# Check headers for the GSE data, may not be same in each dataset
slotNames(gse[[1]])
gse[[1]]@phenoData
varLabels(gse[[1]]@phenoData)

# This gives us the title which contain the time points of each GSM
# Below are several titles for accessing time point info, need to manually enter the title for a specific dataset

# gsm_titles <- gse[[1]]@phenoData$`time point:ch1`
#gsm_titles <- gse[[1]]@phenoData$`time:ch1`
# gsm_titles <- gse[[1]]@phenoData$`differentiation stage:ch1`
#gsm_titles <- gse[[1]]@phenoData$`description`
# gsm_titles <- gse[[1]]@phenoData$`Age:ch1`
gsm_titles <- gse[[1]]@phenoData$`characteristics_ch1`

# Code to extract the hours from time point data when format is like "{number}hours_..."  and convert to days
gsm_times <- as.numeric(gsub("[^-0-9\\.]", "", gsm_titles))   # extract the numeric values
gsm_times <- gsm_times / 24   # convert hours to days


# poor formatting means we will often have to manually set times based on information on GEO or in gsm_titles
gsm_times
gsm_times <- c(-1,-1,-1,0,0,0,0.25,0.25,0.25,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5)


# Annotate expression data by replacing row names (probe IDs) with gene symbols
gene_dat <- fData(gse[[1]])
gene_dat

# Get gene symbols which will be used later in formatting
# Below are several names for gene symbol columns, will need to select one manually.
g_symbols <- fData(gse[[1]])[,'Gene Symbol']
# g_symbols <- fData(gse[[1]])[,'GENE_SYMBOL']
# g_symbols <- fData(gse[[1]])[,'Symbol']

# Now need to annotate the expression data
exprs_df <- exprs(gse[[1]])

# Certain data have treatments, need to remove these columns first based on title, then set cols to the times
# The following code removes columns for any ethanol treated samples. 
# This can be replaced with other terms, or column numbers can be selected based on 'phenoData' or GEO informations (see line 81)
title_vec <- as.vector(gsm_titles)  # we collected title from phenoData, convert to vector
colnames(exprs_df) <- title_vec   # set as column names

col_to_keep <- grep("EtOH", colnames(exprs_df), invert=TRUE)  # this grabs all columns which DO NOT contain "EtOH"
exprs_df <- exprs_df[, col_to_keep]   # subset df without treatment

# Check that lengths of our variables match
length(g_symbols)
length(gsm_times)
nrow(exprs_df)
ncol(exprs_df)

# set new names
rownames(exprs_df) <- g_symbols
colnames(exprs_df) <- gsm_times

# convert to dataframe object
exprs_df <- as.data.frame(exprs_df)
exprs_df

# Save full expression set
##  NOTE: MAY WANT TO EXCLUDE SOME IRRELEVANT DATA COLUMNS HERE
exprs_df <- exprs_df[, 1:(ncol(exprs_df) - 6)]    # Select relevant subset of data
exprs_df

# Name and save the file in our desired directory
f_name <- paste0("data/", gse_id, "_allExpression_data.csv")
write.csv(exprs_df, file=f_name)
