# Human-gut-microbiome-in-IBD-a-network-based-approach
This repository contains basic R code used to reproduce the network analysis presented in the manuscript:
"The human gut microbiota in IBD, characterizing hubs, the core microbiota and terminal nodes: a network-based approach."

We apply network-based methods to analyze 16S rRNA gene sequencing data of the human gut microbiome in the context of inflammatory bowel disease (IBD), comparing local and global network properties, as well as identifying hub and terminal nodes in IBD patients and healthy controls.

Data Availability:
Due to data availability constraints, the dataset used in this analysis cannot be uploaded here. However, the code is structured to work with any two input .RDS files, here referred to as BLcases and BLcontrols.
Expected Data Format:
Each input file should be a data frame with the following structure:
The first column contains the sample IDs as character strings. All subsequent columns contain genus-level abundance counts as integers. Column names for genera should follow the format: G__GenusName. Row-wise: each row represents a sample.
