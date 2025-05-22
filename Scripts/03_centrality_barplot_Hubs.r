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

# use nw from before
nwComp <- readRDS("NW_construct_case_cont_BL.RDS")

options(repr.plot.width=15, repr.plot.height=15)

props_nwComp <- netAnalyze(nwComp,
            # if true (default) computes cent values only for largetest connected component (LCC)
            centrLCC = FALSE,
            clustMethod = "cluster_fast_greedy",
            hubQuant = 0.90,          
            hubPar = c("betweeness", "degree", "closeness"),
            weightDeg = TRUE,
            normDeg = TRUE,
            normBetw = TRUE,
            normClose = TRUE,
            normEigen = TRUE,
            verbose=1,
  )


summary <- summary(props_nwComp,
       groupNames = c("cases", "controls"),
       showCentr = c("betweenness", "closeness", "eigenvector", "degree"),
       numbNodes = 25L
       )
summary

# look at hubs

setdiff(summary$hubs$cases, summary$hubs$controls)
setdiff(summary$hubs$controls, summary$hubs$cases)
intersect(summary$hubs$cases, summary$hubs$controls)

#cases are 0 to 25

#degree from hubs cases
casesDeg <- summary$central$degree1
casesDeg$cases <- as.numeric(casesDeg$cases)
casesDeg <- casesDeg[c(1:25),]
#do min mx normalization to be able to compare cent values of different measures!!!
casesDeg$cases <- (casesDeg$cases - min(casesDeg$cases)) / (max(casesDeg$cases) - min(casesDeg$cases))
casesDeg <- casesDeg[casesDeg$` ` %in% summary$hubs$cases,c(1,2)]
colnames(casesDeg)[1] <- "genera"
colnames(casesDeg)[2] <- "degree"
casesDeg

#betweenness 
casesB <- summary$central$between1
casesB$cases <- as.numeric(casesB$cases)
casesB <- casesB[c(1:25),]
casesB$cases <- (casesB$cases - min(casesB$cases)) / (max(casesB$cases) - min(casesB$cases))
casesB <- casesB[casesB$` ` %in% summary$hubs$cases,c(1,2)]
colnames(casesB)[1] <- "genera"
colnames(casesB)[2] <- "betw"
casesB

#closeness
casesC <- summary$central$close1
casesC$cases <- as.numeric(casesC$cases)
casesC <- casesC[c(1:25),]
casesC$cases <- (casesC$cases - min(casesC$cases)) / (max(casesC$cases) - min(casesC$cases))
casesC <- casesC[casesC$` ` %in% summary$hubs$cases,c(1,2)]
colnames(casesC)[1] <- "genera"
colnames(casesC)[2] <- "close"
casesC

#merge
case_prev <- merge(casesDeg, casesB, by= "genera", sort=TRUE)
case <- merge(case_prev, casesC, by= "genera", sort=TRUE)
case$cases <- as.numeric(case$degree+case$betw+case$close)
case

# same for controlss

#controls are 27 to 51

#degree from hubs controls
controlsDeg <- summary$central$degree1
controlsDeg$controls <- as.numeric(controlsDeg$controls)
controlsDeg <- controlsDeg[c(27:51),]
#do min mx normalization to be able to compare cent values of different measures!!!
controlsDeg$controls <- (controlsDeg$controls - min(controlsDeg$controls)) / (max(controlsDeg$controls) - min(controlsDeg$controls))
controlsDeg <- controlsDeg[controlsDeg$` ` %in% summary$hubs$controls,c(1,3)]
colnames(controlsDeg)[1] <- "genera"
colnames(controlsDeg)[2] <- "degree"
controlsDeg

#betweenness 
controlsB <- summary$central$between1
controlsB$controls <- as.numeric(controlsB$controls)
controlsB <- controlsB[c(27:51),]
controlsB$controls <- (controlsB$controls - min(controlsB$controls)) / (max(controlsB$controls) - min(controlsB$controls))
controlsB <- controlsB[controlsB$` ` %in% summary$hubs$controls,c(1,3)]
colnames(controlsB)[1] <- "genera"
colnames(controlsB)[2] <- "betw"
controlsB

#closeness
controlsC <- summary$central$close1
controlsC$controls <- as.numeric(controlsC$controls)
controlsC <- controlsC[c(27:51),]
controlsC$controls <- (controlsC$controls - min(controlsC$controls)) / (max(controlsC$controls) - min(controlsC$controls))
controlsC <- controlsC[controlsC$` ` %in% summary$hubs$controls,c(1,3)]
colnames(controlsC)[1] <- "genera"
colnames(controlsC)[2] <- "close"
controlsC

controls_prev <- merge(controlsDeg, controlsB, by= "genera", sort=TRUE)
controls <- merge(controls_prev, controlsC, by= "genera", sort=TRUE)
controls$controls <- as.numeric(controls$degree+controls$betw+controls$close)
controls

#copy and paste cd and uc stuff
cd <- read.csv("cdCents.csv")
cd
uc <- read.csv("ucCents.csv")
uc


#also need infromation on cent values for non hubs

# want merges dataframes cases with contols, fill up centrality values for no hubs in different color


case_control <- merge(controls[,c("genera", "controls")], case[,c("genera", "cases")], by= "genera", all=TRUE)
case_control
case_control_cd <- merge(case_control, cd[,c("genera", "cd")], by= "genera", all=TRUE)
case_control_cd_uc <- merge(case_control_cd, uc[,c("genera", "uc")], by= "genera", all=TRUE)
rownames(case_control_cd_uc) <- case_control_cd_uc$genera
case_control_cd_uc$genera <- NULL

case_control_cd_uc <- round(case_control_cd_uc, digits=2)
#case_control_cd_uc

#fill NA values

rownames(case_control_cd_uc) <- gsub("^G__", "", rownames(case_control_cd_uc))
case_control_cd_uc
#save for manuscipt
#write.csv(hubs_cents, file = "hubsCents.csv", row.names = TRUE)

## make bar plots
###  controls ###

library(ggplot2)
library(scales)

### COMPARABLE SCALE OF VALUES FOR 3 CENTS NEEDED, ESP DEGREE


#make plot bigger
options(repr.plot.width=15, repr.plot.height=15)

#define  colors 
Tableau20 <- c("#edc948", "#f4a546", "#bda293")

# Reshape data into long format
df_long <- tidyr::pivot_longer(controls, cols = c("degree", "betw", "close"), names_to = "centrality", values_to = "value")
df_long$taxon <- substring(df_long$genera, 4)
# Create bar plot
ggplot(df_long, aes(x = reorder(taxon,-value), y = value, fill = centrality)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Genera") +
  ylab("Centrality") +
  ggtitle("Controls") +
  theme_minimal() +
  scale_fill_manual(values = Tableau20) +
  theme(legend.text = element_text(size=22), 
        legend.title = element_text(size=22),
        axis.text.x = element_text(size = 25, angle=75, face = "bold"),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 25, vjust = 1), # increase y-axis label size
        axis.title.x = element_text(size = 25, vjust=15), # increase y-axis tick mark size
        plot.title = element_text(size = 35))+ # increase title size
        guides(fill="none")
#ggsave("controlCentsBar.png", width=15, height=15, dpi=300)

## make bar plots
###  controls ###

library(ggplot2)
library(scales)

### COMPARABLE SCALE OF VALUES FOR 3 CENTS NEEDED, ESP DEGREE
#make plot bigger
options(repr.plot.width=15, repr.plot.height=15)

#define  colors 
Tableau20 <- c("#edc948", "#f4a546", "#bda293")
# Reshape data into long format
df_long <- tidyr::pivot_longer(case, cols = c("degree", "betw", "close"), names_to = "centrality", values_to = "value")
df_long$taxon <- substring(df_long$genera, 4)
# Create bar plot
ggplot(df_long, aes(x = reorder(taxon,-value), y = value, fill = centrality)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Genera") +
  ylab("Centrality") +
  ggtitle("Cases") +
  theme_minimal() +
  scale_fill_manual(values = Tableau20) +
  theme(legend.text = element_text(size=22), 
        legend.title = element_text(size=22),
        axis.text.x = element_text(size = 25, angle=75, face = "bold"),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 25, vjust = 1), # increase y-axis label size
        axis.title.x = element_text(size = 25, vjust=15), # increase y-axis tick mark size
        plot.title = element_text(size = 35))+ # increase title size
        guides(fill="none")
#ggsave("caseCentsBar.png", width=15, height=15, dpi=300)


