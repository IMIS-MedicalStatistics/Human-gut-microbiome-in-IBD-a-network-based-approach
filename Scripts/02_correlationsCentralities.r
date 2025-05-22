## load libraries, don't need all here
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

# for paper always for with same network 

nwComp_read <- readRDS("NW_construct_case_cont_BL.RDS")

# netAnalyze
# will automatically create GCMs
file_name <- "GCM_cases_controls_LCC.png" # or whatever name one would prefer
#png(file_name, width = 4200, height = 4200, res = 300) 

# make figures bigger
options(repr.plot.width=15, repr.plot.height=15)

props_nwComp <- netAnalyze(nwComp_read,
            # if true (default) computes cent values only for largetest connected component (LCC)
            centrLCC = FALSE,
            
            #cluster algo
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
dev.off()

comp_nw <- readRDS("NW_comp_case_cont_BL_1000L.RDS") # saved from NetwrokConstruction script

sumComp <- summary(comp_nw,
       groupNames = c("cases", "controls"),
       showCentr = c("betweenness", "closeness", "degree"),
       numbNodes = 10L
       )
sumComp


## controls (group 2)

str(sumComp$properties$deg2)
data <- data.frame(close = sumComp$properties$close2,
                   betw = sumComp$properties$betw2,
                  deg = sumComp$properties$deg2)

# Compute correlations
cor_cb <- cor(data$close, data$betw, method = "spearman", use = "complete.obs")
cor_db <- cor(data$deg, data$betw, method = "spearman", use = "complete.obs")
cor_dc <- cor(data$deg, data$close, method = "spearman", use = "complete.obs")

# Create labels with subscript syntax
cor_labels <- c(
  paste("ρ[clos_betw] == ", round(cor_cb, 2)),
  paste("ρ[deg_betw] == ", round(cor_db, 2)),
  paste("ρ[deg_clos] == ", round(cor_dc, 2))
)

# Label positions (stacked vertically at top-left)
label_data <- data.frame(
  x = min(log(data$close), na.rm = TRUE),
  y = seq(from = max(data$betw, na.rm = TRUE), by = -0.05 * diff(range(data$betw, na.rm = TRUE)), length.out = 3),
  label = cor_labels
)

# Create a scatter plot with ggplot2
a <- ggplot(data, aes(x = log(close), y = betw, color = deg)) +
  geom_point(shape = 16, size = 5) +  # Set fill color based on degree variable
  scale_color_viridis_c(option = "plasma") +  # Use viridis color scale
  labs(title = "Controls",
       x = "log-scaled Closeness Centrality",
       y = "Betweenness Centrality",
       fill = "Degree") +
  geom_text(data = label_data, aes(x = x, y = y, label = label),inherit.aes = FALSE,
            hjust = -0.3, vjust = 1.0, size = 6, fontface = "italic", parse = TRUE) +
  theme_minimal() +
  theme(plot.title = element_text(size = 28),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))
plot(a)
#ggsave("ClosBetwControls.png2.png", width=10, height=10, dpi=300)


### cases (group 1)

str(sumComp$properties$deg1)
data <- data.frame(close = sumComp$properties$close1,
                   betw = sumComp$properties$betw1,
                  deg = sumComp$properties$deg1)

# Compute correlations
cor_cb <- cor(data$close, data$betw, method = "spearman", use = "complete.obs")
cor_db <- cor(data$deg, data$betw, method = "spearman", use = "complete.obs")
cor_dc <- cor(data$deg, data$close, method = "spearman", use = "complete.obs")

# Create labels with subscript syntax
cor_labels <- c(
  paste("ρ[clos_betw] == ", round(cor_cb, 2)),
  paste("ρ[deg_betw] == ", round(cor_db, 2)),
  paste("ρ[deg_clos] == ", round(cor_dc, 2))
)
# Label positions (stacked vertically at top-left)
label_data <- data.frame(
  x = min(log(data$close), na.rm = TRUE),
  y = seq(from = max(data$betw, na.rm = TRUE), by = -0.05 * diff(range(data$betw, na.rm = TRUE)), length.out = 3),
  label = cor_labels
)

# Create a scatter plot with ggplot2
b <- ggplot(data, aes(x = log(close), y = betw, color = deg)) +
  geom_point(shape = 16, size = 5) +  # Set fill color based on degree variable
  scale_color_viridis_c(option = "plasma") +  # Use viridis color scale
  labs(title = "Cases",
       x = "log-scaled Closeness Centrality",
       y = "Betweenness Centrality",
       fill = "Degree") +
 geom_text(data = label_data, aes(x = x, y = y, label = label),inherit.aes = FALSE,
            hjust = -0.3, vjust = 1.0, size = 6, fontface = "italic", parse = TRUE) +
  theme_minimal() +
  theme(plot.title = element_text(size = 28),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))
#ggsave("ClosBetwCases2.png", width=10, height=10, dpi=300)


