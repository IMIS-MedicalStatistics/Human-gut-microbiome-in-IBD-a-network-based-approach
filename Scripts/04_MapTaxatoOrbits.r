## load libraries
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



## map taxa back to orbits from gcms

### based on output from netANalyze when having to input data (cases and controls), here called props_nwComp

# map taxa back to orbits


### cases 1


top=25 # looking at the top 25 taxa

# Convert the matrix to a data frame, take row new_id out [-1,]
orbit_df <- as.data.frame(props_nwComp$graphlet$ocount1[-1,])

# List to store top 10 genera for each orbit
top_by_orbit <- list()

# Specify the orbits of interest, end  node orbits
selected_orbits <- c("O1","O4", "O6", "O9")

# Iterate over each orbit
for (orbit in selected_orbits) {
  # Subset the data frame to include only the current orbit
  selected_df <- subset(orbit_df, select = orbit)
  
  # Calculate the sum of occurrences for each genus in the current orbit, rowSums doesnt really do anything here, just for safety test
  genus_occurrences <- rowSums(selected_df)
  
  # Sort genera based on occurrences
  sorted_genera <- names(sort(genus_occurrences, decreasing = TRUE))
  
  # Select the top genera for the current orbit
  top_genera <- sorted_genera[1:top]
  
  # Store the top genera for the current orbit in the list
  top_by_orbit[[orbit]] <- top_genera
}

# Print the top  genera for each orbit
for (orbit in selected_orbits) {
  cat("Top",top, "genera for orbit", orbit, ":", paste(top_by_orbit[[orbit]], collapse = ", "), "\n")
}


# Create an empty data frame to store the top genera for each orbit
top_df <- data.frame(matrix(ncol = length(selected_orbits), nrow = top))
colnames(top_df) <- selected_orbits

# Fill the data frame with the top genera for each orbit
for (i in 1:length(selected_orbits)) {
  orbit <- selected_orbits[i]
  top_df[, i] <- top_by_orbit[[orbit]]
}

# Print the data frame
top_df

# Check for common genera across all four orbits
common_genera_cases <- Reduce(intersect, top_by_orbit)

# Print the common genera
cat("Common genera across all four orbits:", paste(common_genera_cases, collapse = ", "))



## similar for controls , group 2


top=25

# Convert the matrix to a data frame, take row new_id out
orbit_df <- as.data.frame(props_nwComp$graphlet$ocount2[-1,])

# List to store top genera for each orbit
top_by_orbit <- list()

# Specify the orbits of interest
selected_orbits <- c("O1","O4", "O6", "O9")

# Iterate over each orbit
for (orbit in selected_orbits) {
  # Subset the data frame to include only the current orbit
  selected_df <- subset(orbit_df, select = orbit)
  
  # Calculate the sum of occurrences for each genus in the current orbit
  genus_occurrences <- rowSums(selected_df)
  
  # Sort genera based on occurrences
  sorted_genera <- names(sort(genus_occurrences, decreasing = TRUE))
  
  # Select the top genera for the current orbit
  top_genera <- sorted_genera[1:top]
  
  # Store the top genera for the current orbit in the list
  top_by_orbit[[orbit]] <- top_genera
}

# Print the top genera for each orbit
for (orbit in selected_orbits) {
  cat("Top", top,"genera for orbit", orbit, ":", paste(top_by_orbit[[orbit]], collapse = ", "), "\n")
}


# Create an empty data frame to store the top genera for each orbit
top_df <- data.frame(matrix(ncol = length(selected_orbits), nrow = top))
colnames(top_df) <- selected_orbits

# Fill the data frame with the top genera for each orbit
for (i in 1:length(selected_orbits)) {
  orbit <- selected_orbits[i]
  top_df[, i] <- top_by_orbit[[orbit]]
}

# Print the data frame
top_df

# Check for common genera across all four orbits
common_genera_controls <- Reduce(intersect, top_by_orbit)

# Print the common genera
cat("Common genera across all four orbits:", paste(common_genera_controls, collapse = ", "))



# compare top 30, therefrom select intersection in all 4 orbits, and now difference between cases and control
diff_genera1 <- setdiff(common_genera_cases, common_genera_controls)
diff_genera2 <- setdiff(common_genera_controls, common_genera_cases)
diff_genera <- c(diff_genera1, diff_genera2)
diff_genera1 # in cases but not controls
diff_genera2 # in controls but not cases

length(diff_genera)
length(common_genera_cases)
length(common_genera_controls)

## 20 from 26 different

same_genera1 <- intersect(common_genera_cases, common_genera_controls)
same_genera2 <- intersect(common_genera_controls, common_genera_cases)
same_genera <- c(same_genera1, same_genera2)
same_genera # 6 total
length(same_genera)

# Determine the length of the longer variable

# Remove "G__" prefix from each character in diff_genera1
diff_genera1 <- gsub("^G__", "", diff_genera1)

# Remove "G__" prefix from each character in diff_genera2
diff_genera2 <- gsub("^G__", "", diff_genera2)

max_length <- max(length(diff_genera1), length(diff_genera2))

# Create a vector of the same length as diff_genera2 filled with NA values
diff_genera1_filled <- c(diff_genera1, rep(" ", max_length - length(diff_genera1)))

# Create the data frame
diff_genera_df <- data.frame("Only in cases" = diff_genera1_filled, "Only in controls" = diff_genera2)

# Display the data frame
diff_genera_df