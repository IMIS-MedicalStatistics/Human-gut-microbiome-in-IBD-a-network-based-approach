#load libraries
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
library(microbiome)

### START WITH OLD DEFINITION OF CORE (AUNDANCE AND PREVALENCE)

library(phyloseq)
library(microbiome)
library(microViz)

# create phyloseq object

OTU <- otu_table(otu_asv, taxa_are_rows = FALSE) # asv table input
TAX  <- tax_table(tax)                          # taxonomy table
META <-  sample_data(meta)                      #metadata
rownames(TAX) <- rownames(tax)      
colnames(TAX) <- colnames(tax)

physeq <- phyloseq(OTU,TAX,META)
physeq

# filtering, like in NW stuff

#first aggregate on genus level
pgc<- aggregate_taxa(physeq, "Genus")

# then split in groups, here cases and controls, disease status from metadata
pgc_cases <- subset_samples(pgc, disease == "case")
pgc_controls <- subset_samples(pgc, disease == "controls")
pgc_cases
pgc_controls

##samples

# Calculate total reads per sample
total_reads_cases <- sample_sums(pgc_cases)
total_reads_controls <- sample_sums(pgc_controls)

# Set a threshold for minimum total reads
min_total_reads <- 1000  # Set your desired threshold here

# Identify samples with total reads greater than or equal to the threshold
samples_to_keep_cases <- total_reads_cases >= min_total_reads
samples_to_keep_controls <- total_reads_controls >= min_total_reads

# Filter samples based on total reads
pgc_cases_filtered_samples <- prune_samples(samples_to_keep_cases, pgc_cases)
pgc_controls_filtered_samples <- prune_samples(samples_to_keep_controls, pgc_controls)

##taxa

#at least observes in 4 samples (1% min lib size (see nw stuff) from smaller subgroup )
#physeq_filtered, min prevalence and min total abundance

pgc_cases_filtered_samples_taxa <- pgc_cases_filtered_samples %>% tax_filter(
                                                    min_prevalence = 4,
                                                    min_total_abundance = 0.0001
                                                        )
pgc_controls_filtered_samples_taxa <- pgc_controls_filtered_samples %>% tax_filter(
                                                    min_prevalence = 4,
                                                    min_total_abundance = 0.0001
                                                        )
# do clr trafo
pgc_cases_filter <- microbiome::transform(pgc_cases_filtered_samples_taxa, "clr")
pgc_controls_filter <- microbiome::transform(pgc_controls_filtered_samples_taxa, "clr")

pgc_cases_filter
pgc_controls_filter

### standard definition
# detection: threshold presence/absence in samples
# prevalence: threshold in [0,1]


core_genera_cases <- microbiome::core_members(pgc_cases_filter,
                                              detection = 0.1/100, prevalence = 50/100)
core_genera_cases
length(core_genera_cases)


core_genera_controls <- microbiome::core_members(pgc_controls_filter,
                                              detection = 0.1/100, prevalence = 50/100)
core_genera_controls
length(core_genera_controls)

### NOW DO NW STUFF DEFINITION OF CORE: HUBS (CENTRALITY MEASURES HIGHEST VALUES)

# load preprocessed data, HERE ONLY UNRELATED SUBJECTS

BLcases <- readRDS("cases.RDS")
BLcontrolls <- readRDS("controlls.RDS")
dim(BLcontrolls)
dim(BLcases)#cut off the G__ in front
colnames(BLcases) <- sub("^G__", "", colnames(BLcases))
colnames(BLcontrolls) <- sub("^G__", "", colnames(BLcontrolls))

# use same nw trough whole analysis

nwComp <- readRDS("NW_construct_case_cont_BL.RDS") # netwrok from before

# here adjust hub quantile value to match taxa numbers in old def core
# for 0.54 the sum of hubs for cases and controls is 76, just as for old def

props_nwComp <- netAnalyze(nwComp,
            # if true (default) computes cent values only for largetest connected component (LCC)
            centrLCC = FALSE,
            
            #cluster algo
            clustMethod = "cluster_fast_greedy",
            hubQuant = 0.55,          
            hubPar = c("betweeness", "degree", "closeness"),
            weightDeg = TRUE,
            normDeg = TRUE,
            normBetw = TRUE,
            normClose = TRUE,
            normEigen = TRUE,
            verbose=1,
                          )

sum <- summary(props_nwComp,
       groupNames = c("cases", "controls"),
       showCentr = c("betweenness", "closeness", "degree"),
       numbNodes = 2000L
       )
#sum
# new_id in here, delete that!!

sum$hubs <- sum$hubs[sum$hubs$cases != "new_id",]
length(sum$hubs$cases)

sum$hubs$cases <- gsub("^G__", "", sum$hubs$cases)
sum$hubs$controls <- gsub("^G__", "", sum$hubs$controls)

hubs_cases <- sum$hubs$cases
hubs_controls <- sum$hubs$controls
sum$hubs
#delete empty enries in controls, else errors in following stepssum$hubs

# want supplementary table for manuscript, containing all core members (cases and controls def 1 and def 2)
coreMembers <- sum$hubs
rownames(coreMembers) <- NULL
colnames(coreMembers)[colnames(coreMembers) =="cases"] <- "cases definition 2"
colnames(coreMembers)[colnames(coreMembers) =="controls"] <- "controls definition 2"

#next coloumns are shorter, so fill them with ""
# Determine the length of the dataframe and the length of your character vector
num_rows <- nrow(coreMembers)
char_length <- length(core_genera_cases)
# Fill the remaining rows with empty strings
coreMembers$"cases definition 1" <- c(core_genera_cases, rep("", num_rows - char_length))

#bc next def is 40 member, make df longer
coreMembers <- rbind(coreMembers, "", "")
coreMembers$"controls definition 1" <- c(core_genera_controls)

#sort each coloumn alphabetically
coreMembers[coreMembers == ""] <- NA
for (col in colnames(coreMembers)) {
  coreMembers[[col]] <- sort(coreMembers[[col]], na.last=TRUE)
}
coreMembers[is.na(coreMembers)] <- ""

#reorder colojumns
coreMembers <- coreMembers[,c(4,2,3,1)]
coreMembers
#write.csv(coreMembers, file = "coreMembers.csv", row.names = TRUE)



typeof(hubs_controls)
hubs_controls <- hubs_controls[1:(length(hubs_controls) - 3)]
length(hubs_controls)
hubs_controls

### COMPARE


#cases, 30 from 38 different -> a lot!
case_hubs_diff1 <- setdiff(core_genera_cases, hubs_cases)
case_hubs_diff1
length(case_hubs_diff1)
case_hubs_diff2 <- setdiff(hubs_cases, core_genera_cases)
case_hubs_diff2
length(case_hubs_diff2)

#controls, 24 from 38 different -> a lot too!
controls_hubs_diff1 <- setdiff(core_genera_controls, hubs_controls)
controls_hubs_diff1
length(controls_hubs_diff1)
controls_hubs_diff2 <- setdiff(hubs_controls, core_genera_controls)
controls_hubs_diff2
length(controls_hubs_diff2)

## COMPARE FURTHER

#therefore adjust data to containing only respective core taxa

cases_core_oldDef <- BLcases[,core_genera_cases]
cases_core_NWDef1 <- BLcases[, hubs_cases]
#must convert to numeric, why ever...
cases_core_NWDef <- sapply(cases_core_NWDef1, as.numeric)

controls_core_oldDef <- BLcontrolls[,core_genera_controls]
hubs_controls
controls_core_NWDef <- BLcontrolls[, hubs_controls]
controls_core_NWDef <- sapply(controls_core_NWDef, as.numeric)
cases_core_NWDef

# compare total abundance of all core taxa old def vs nw def, seperatly for cases and controls

tot_abundance_casesOldDef <- mean(rowSums(cases_core_oldDef))
tot_abundance_casesOldDef
tot_abundance_casesNWDef <- mean(rowSums(cases_core_NWDef))
tot_abundance_casesNWDef

tot_abundance_controlsOldDef <- mean(rowSums(controls_core_oldDef))
tot_abundance_controlsOldDef
tot_abundance_controlsNWDef <- mean(rowSums(controls_core_NWDef))
tot_abundance_controlsNWDef

# in cases NW core taxa have higher total abundance, whereas in cotnrols old def have higher abundance


# average prevvalence: prevalence of core taxa on average per sample


## cases old

taxa_presence_casesOld <- mean((colSums(cases_core_oldDef != 0))/dim(cases_core_oldDef)[1])
# Calculate the average prevalence

## cases new

taxa_presence_casesNW <- mean((colSums(cases_core_NWDef != 0))/dim(cases_core_NWDef)[1])
# Calculate the average prevalence

## controls old

taxa_presence_controlsOld <- mean((colSums(controls_core_oldDef != 0))/dim(controls_core_oldDef)[1])
# Calculate the average prevalence


## controls old

taxa_presence_controlsNW <- mean((colSums(controls_core_NWDef != 0))/dim(controls_core_NWDef)[1])
# Calculate the average prevalence


## compare sum( centrality values: degree+betw+clos) for each taxa, mean over all taxa

#cases old def

#want them from the nw where all are in tho, dont do a new nw with just those genera
deg <- sum$central$degree1
betw <- sum$central$between1
clos <- sum$central$close1
# have G__ in front, must go


#give it a name
colnames(deg)[1] <- "genera"
colnames(deg)[2] <- "degree"
colnames(deg)[3] <- "degree_control"
colnames(betw)[1] <- "genera"
colnames(betw)[2] <- "betw"
colnames(betw)[3] <- "betw_control"
colnames(clos)[1] <- "genera"
colnames(clos)[2] <- "clos"
colnames(clos)[3] <- "clos_control"
#only need first half, order of cents doesnt matter
betw  <-  betw[1:ceiling(nrow(betw)/2),]
clos <- clos[1:ceiling(nrow(clos)/2),]

deg$genera <- gsub("^G__", "", deg$genera)
betw$genera <- gsub("^G__", "", betw$genera)
clos$genera <- gsub("^G__", "", clos$genera)

CasesCoreOldDef <- deg[deg$genera %in% core_genera_cases,c(1,2)]
CasesCoreOldDef <- CasesCoreOldDef[1:length(core_genera_cases),]
CasesCoreOldDef <- merge(CasesCoreOldDef, betw[,c("genera","betw")], by = "genera", all.x = TRUE)
CasesCoreOldDef <- merge(CasesCoreOldDef, clos[,c("genera","clos")], by = "genera", all.x = TRUE)
rownames(CasesCoreOldDef) <- CasesCoreOldDef$genera
CasesCoreOldDef$genera <- NULL
CasesCoreOldDef <- sapply(CasesCoreOldDef, as.numeric)
CasesCoreOldDef_mean  <- mean(rowSums(CasesCoreOldDef))


## cases new def

CasesCoreNWDef <- deg[deg$genera %in% hubs_cases, c(1,2)]
CasesCoreNWDef <- CasesCoreNWDef[1:length(hubs_cases),]
CasesCoreNWDef <- merge(CasesCoreNWDef, betw[,c("genera","betw")], by = "genera", all.x = TRUE)
CasesCoreNWDef <- merge(CasesCoreNWDef, clos[,c("genera","clos")], by = "genera", all.x = TRUE)
rownames(CasesCoreNWDef) <- CasesCoreNWDef$genera
CasesCoreNWDef$genera <- NULL
CasesCoreNWDef <- sapply(CasesCoreNWDef, as.numeric)
CasesCoreNWDef_mean  <- mean(rowSums(CasesCoreNWDef))


## controls old def

ControlCoreOldDef <- deg[deg$genera %in% core_genera_controls,c(1,3)]
ControlCoreOldDef <- ControlCoreOldDef[1:length(core_genera_controls),]
ControlCoreOldDef <- merge(ControlCoreOldDef, betw[,c("genera","betw_control")], by = "genera", all.x = TRUE)
ControlCoreOldDef <- merge(ControlCoreOldDef, clos[,c("genera","clos_control")], by = "genera", all.x = TRUE)
rownames(ControlCoreOldDef) <- ControlCoreOldDef$genera
ControlCoreOldDef$genera <- NULL
ControlCoreOldDef <- sapply(ControlCoreOldDef, as.numeric)
ControlCoreOldDef_mean  <- mean(rowSums(ControlCoreOldDef))


## controls new def

ControlCoreNWDef <- deg[deg$genera %in% hubs_controls,c(1,3)]
ControlCoreNWDef <- ControlCoreNWDef[1:length(hubs_controls),]
ControlCoreNWDef <- merge(ControlCoreNWDef, betw[,c("genera","betw_control")], by = "genera", all.x = TRUE)
ControlCoreNWDef <- merge(ControlCoreNWDef, clos[,c("genera","clos_control")], by = "genera", all.x = TRUE)
rownames(ControlCoreNWDef) <- ControlCoreNWDef$genera
ControlCoreNWDef$genera <- NULL
ControlCoreNWDef <- sapply(ControlCoreNWDef, as.numeric)
ControlCoreNWDef_mean  <- mean(rowSums(ControlCoreNWDef))


### visualize, grouped bar plot for each measure: total abundance (A), prevalence (P) in %,
# sum centrality measures mean (C)
# first is old def second is NW def

# hard code values in, reuslts from above, here just DUMMY VALUES

A_cases <- c(1,2)
A_controls <- c(3,4)
P_cases <- c(5,6)
P_controls <- c(7,8)
C_cases <- c(9,10)
C_controls <- c(11,12)

### abundance

df <- data.frame(
  Group = rep(c("Cases", "Controls"), each = 2),
  Definition = rep(c("Abundance-based", "Network-based"), times = 2),
  Total_Abundance = c(A_cases, A_controls)
)

# Combine the data into a data frame
df_prevalence <- data.frame(
  Group = rep(c("Cases", "Controls"), each = 2),
  Definition = rep(c("Abundance-based", "Network-based"), times = 2),
  Prevalence_Percentage = c(P_cases, P_controls)
)

# Combine the data into a data frame
df_centrality <- data.frame(
  Group = rep(c("Cases", "Controls"), each = 2),
  Definition = rep(c("Abundance-based", "Network-based"), times = 2),
  Average_Sum_Centrality = c(C_cases, C_controls)
)
df_centrality

# Plotting
#helper for colors
df$colors <- c("Cases Definition 1", "Cases Definition 2", "Controls Definition 1", "Controls Definition 2")
df_prevalence$colors <- c("Cases Definition 1", "Cases Definition 2", "Controls Definition 1", "Controls Definition 2")
df_centrality$colors <- c("Cases Definition 1", "Cases Definition 2", "Controls Definition 1", "Controls Definition 2")

coreMeasureA <- ggplot(df, aes(x = Group, y = Total_Abundance, fill = colors)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "Cumulative Abundance") +  # Remove title
  scale_y_continuous(labels = scales::comma) +  # Format y-axis labels without scientific notation
  scale_fill_manual(values = c("#4e79a7", "#a0cbe8", "#e78ac3", "#d37295")) +  # Assign colors 2,4,1,3
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),  # Adjust the size of x-axis labels
        axis.text.y = element_text(size = 12),  # Adjust the size of y-axis labels
        axis.title = element_text(size = 14),  # Adjust the size of axis labels
        plot.title = element_text(size = 16))+
guides(fill = FALSE)  # Remove legend# Adjust the size of plot title
ggsave("CoreMeasuresA_new.png", coreMeasureA, width = 6, height = 6, units = "in", dpi = 300)
coreMeasureA

coreMeasureP <- ggplot(df_prevalence, aes(x = Group, y = Prevalence_Percentage, fill = colors)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "Average Prevalence") +
  scale_y_continuous(labels = scales::percent) +
scale_fill_manual(values = c("#4e79a7", "#a0cbe8", "#e78ac3", "#d37295")) +  # Assign colors 2,4,1,3
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16))+
guides(fill = FALSE)  # Remove legend
ggsave("CoreMeasuresP_new.png", coreMeasureP, width = 6, height = 6, units = "in", dpi = 300)
coreMeasureP

# Plotting for Average of Sum Centrality Measures
coreMeasureC <- ggplot(df_centrality, aes(x = Group, y = Average_Sum_Centrality, fill = colors)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "Average Sum Centrality") +
scale_fill_manual(values = c("#4e79a7", "#a0cbe8", "#e78ac3", "#d37295")) +  # Assign colors 2,4,1,3
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16))+
guides(fill = FALSE)  # Remove legend
ggsave("CoreMeasuresC_new.png", coreMeasureC, width = 6, height = 6, units = "in", dpi = 300)

# Define the values for each group and definition
Cases_Definition1 <- c(11689795, 0.0651340996168582, 13.2025420588235)
Controls_Definition1 <- c(7586669, 0.115068493150685, 4.39063238095238)
Cases_Definition2 <- c(11386311, 0.0747126436781609, 14.8795353846154)
Controls_Definition2 <- c(6532096, 0.104109589041096, 5.28060315789474)

# Create a data frame for plotting
df <- data.frame(
  Measure = rep(c("Total Abundance", "Prevalence Percentage", "Sum Centrality Measures Mean"), 4),
  Group = rep(c("Cases", "Controls"), each = 2),
  Definition = rep(c("Definition 1", "Definition 2"), times = 2),
  Value = c(Cases_Definition1, Controls_Definition1, Cases_Definition2, Controls_Definition2)
)
df

# visualize different core member, venn diagramm
library(VennDiagram)
library(ggvenn)
# Create a Venn diagram with taxa names displayed
x = list("Definition 1" = core_genera_controls, "Definition 2" = hubs_controls)
venn.plot <- venn.diagram(x,
  category.names = c("Definition 1", "Definition 2"),
  filename = NULL,
)

# Plot the Venn diagram
#grid.draw(venn.plot)

vennControls <- ggvenn(x, show_elements = F, label_sep = "\n", text_size = 6, fill_color = c("#e78ac3", "#d37295")
              , auto_scale = TRUE,stroke_size = 1)+
          ggtitle("Controls")+
          theme(plot.title = element_text(hjust = 0.5, vjust = -10, size = 20))

vennControls

ggsave("venn_diagramControls.png", vennControls, width = 8, height = 6, units = "in", dpi = 300)


# visuaize cor members, cases and controls
x = list("Definition 1" = core_genera_cases, "Definition 2" = hubs_cases)
venn.plot <- venn.diagram(x,
  category.names = c("Definition 1", "Definition 2"),
  filename = NULL,
)

vennCases <- ggvenn(x, show_elements = F, label_sep = "\n", text_size = 6, fill_color = c("#4e79a7", "#a0cbe8")
              , auto_scale = TRUE,stroke_size = 1)+
          ggtitle("Cases")+
          theme(plot.title = element_text(hjust = 0.5, vjust = -10, size = 20))

vennCases

ggsave("venn_diagramCases.png", vennCases, width = 8, height = 6, units = "in", dpi = 300)



