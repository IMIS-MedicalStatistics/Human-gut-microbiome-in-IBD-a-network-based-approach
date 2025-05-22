#load libraries, basically all thar are needed for everything
library(plyr)
library(devtools)
library(ggplot2)
library(gtools) 
library(tidyverse) # for plotting and wrangling data
library(SpiecEasi) # Has sparcc and also does clr transforms
library(otuSummary)
#library(dplyr) 
library(reshape2) # has the melt funct  ion, which I use to wrangle data

library(psych) # for calculating regular correlations with p values
library(glue) # to get sth comparable to fstrings in python
library(NetCoMi)
library(scales)
library(janitor)

library(scales)
library(vegan)
library(ggpubr)
library(fossil)


### for more detailed describtion see NetCoMi: https://rdrr.io/github/stefpeschel/NetCoMi/




#construct nw ibd vs non ibd 
#here all asmples with na values (also covariables) get removed, no need to take care of before by hand


##### FILTERING ADAPTED
#filter samples: total reads at least 10000
#filter taxa: number of reads a.l. x% of total reads 
#PLUS numb Samp (observes al in x samples), x being 1% of minimum library size


nwComp <- netConstruct(
    data = BLcases,
    data2 = BLcontrolls,
    # or pearson, sparcc, cclasso,...
    measure = "sparcc",
    normMethod = "none", 
    sparsMethod = "t-test", 
    alpha = 0.05,
    adjust= "fdr",
                   
    verbose = 1,
                   
    # dont put in absolute parameters, bc dataframe size will differ
    # if two groups differ use smaller size
    filtSamp = "totalReads",
    filtSampPar = list(totalReads = 10000),
    filtTax = c("relFreq", "numbSamp"), 
    filtTaxPar = list(relFreq= 0.0001, numbSamp = 0.01*dim(BLcontrolls)[1]),
                  )
saveRDS(nwComp, file="NW_construct_case_cont_BL.RDS", compress=F) always for with that nw from now on

# netAnalyze
# will automatically create GCMs
file_name <- "GCM_cases_controls_LCC.png" # no need to chose the same name
#png(file_name, width = 4200, height = 4200, res = 300) 

# make figures bigger
options(repr.plot.width=15, repr.plot.height=15)

props_nwComp <- netAnalyze(nwComp_read,
            # if true (default) computes cent values only for largetest connected component (LCC)
            centrLCC = FALSE,
            
            #cluster algo
            clustMethod = "cluster_fast_greedy", # choce of clustering algortihm
            hubQuant = 0.90,          # 90th percentile
            hubPar = c("betweeness", "degree", "closeness"), # chooce of centrality measures
            weightDeg = TRUE,
            normDeg = TRUE,
            normBetw = TRUE,
            normClose = TRUE,
            normEigen = TRUE,
            verbose=1,
                          )
dev.off()

#netCompare: computational high demand

comp_nw <- netCompare(props_nwComp,
           permTest = TRUE,
           nPerm = 1000L,  # 1000 permutations
           verbose = FALSE,  
           adjust = "fdr", # adjust for multiple testing
           jaccQuant = 0.90,
                     )

saveRDS(comp_nw, file="NW_comp_case_cont_BL_1000L.RDS", compress=F) # do save, because comutational high demand (days on cluster)

comp_nw <- readRDS("NW_comp_case_cont_BL_1000L.RDS")

sumComp <- summary(comp_nw,
       groupNames = c("cases", "controls"),
       showCentr = c("betweenness", "closeness", "degree"),
       numbNodes = 10L
       )
sumComp

# on netAnalyze

sum <- summary(props_nwComp,
       groupNames = c("cases", "controls"),
       showCentr = c("betweenness", "closeness", "degree"),
       numbNodes = 100L
       )
sum

# Set plot size for Jupyter notebook rendering
options(repr.plot.width=25, repr.plot.height=25)

# Output file name for the plot
file_name <- "revised_cases_controls_NWPic.png"
# Uncomment to save the plot to a file
#png(file_name, width = 4200, height = 4200, res = 300)

# Plot the network
plot(props_nwComp, 
     layout = "spring",
     groupNames = c("cases", "controls"),
     sameLayout = TRUE,
     repulsion = 0.5, # Node repulsion for layout
     shortenLabels = "simple",
     labelLength = 8,
     charToRm = "G__",
     labelFont = 4,
     
     ## Node Settings
     nodeColor = "cluster",  # node color represents cluster they belong to
     colorVec = c("deepskyblue4", "orange", "purple3", "aquamarine", "brown", "red"),
     nodeSize = "normCounts", # Size of nodes based on normalized counts
     highlightHubs = TRUE, # Highlight hub nodes with borders
     hubBorderWidth = 10, # Border width for hubs
     labelScale = FALSE,
     cexLabels = 1.8, # Label font size
     cexNodes = 8, # Node size
     cexTitle = 2, # Title font size
     
     ## Edge Settings
     edgeFilter = "threshold",
     edgeFilterPar = 0.15,
)

legend("topleft", inset = c(0.385, +0.79), 
       legend = c("Low normalized counts, cluster 1", "Medium normalized counts, cluster 1",
                  "High normalized counts, cluster 2"),
       col = "black", pt.cex = c(4, 7, 10), pch = 21, pt.bg = c("purple3", "purple3", "aquamarine"), 
       cex = 1.2, bty = "n")


legend("topleft", inset = c(0.35, +0.0), 
       legend = c("Low positive correlation", "High negative correlation"), 
       lwd = c(1, 5), 
       col = c("green", "red"),  # Green for Low, Red for High
       cex = 1.2, bty = "n")#, 
       #title = "Edge Thickness: Absolute Correlation Value", text.width = 3)


# Close the plot device (if saving to file)
#dev.off()

