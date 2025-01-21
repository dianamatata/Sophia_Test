
# Load necessary libraries and data  -----------------------
library(ggplot2)
library(dplyr)
library(tidyr)


file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/mobi_data_omim_splice_clinvarentries.txt"
mobi_data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# plot MPA -----------------------

print(colnames(mobi_data))
mobi_data$"MPA.score."  <- as.numeric(mobi_data$"MPA.score.")
unique_mpa_scores <- unique(mobi_data$`MPA.score.`)
print(unique_mpa_scores)

# Round the MPA scores to the nearest integer between 0 and 10, and handle NAs
mobi_data$bin <- ifelse(is.na(mobi_data$`MPA.score.`), NA, round(mobi_data$`MPA.score.`))

# Ensure bins are between 0 and 10 (including handling NAs)
mobi_data <- mobi_data %>%
  mutate(bin = ifelse(bin < 0 | bin > 10, NA, bin))

# Convert MPA scores to bins
mobi_data$bin <- ifelse(is.na(mobi_data$`MPA.score.`), "NA", as.character(round(mobi_data$`MPA.score.`)))

# Convert bin to a factor and make sure "NA" is one of the levels
mobi_data$bin <- factor(mobi_data$bin, levels = c(as.character(0:10), "NA"))

# Calculate bin counts
bin_counts <- mobi_data %>%
  group_by(bin) %>%
  summarise(count = n(), .groups = "drop")

# Assign colors based on bin value
bin_counts$color <- ifelse(bin_counts$bin %in% c("8", "9", "10"), "red", 
                           ifelse(is.na(bin_counts$bin), "grey", "blue"))

# Create the plot with color conditionality
plot <- ggplot(bin_counts, aes(x = bin, y = count, fill = color)) +
  geom_bar(stat = "identity", color = 'black', alpha = 0.7) +
  geom_text(aes(label = count), vjust = -0.5, size = 5) +  # Display count on top of each bar
  labs(title = 'Histogram of MPA Score', x = 'MPA Score', y = 'Frequency') +
  theme_minimal() +
  scale_x_discrete(breaks = c(as.character(0:10), "NA")) +  # Ensure x-axis is discrete with string labels
  scale_fill_identity() +  # Use the predefined fill colors
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Optional: Rotate x-axis labels


ggsave("/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/plots/MPA_score_histogram.png", plot = plot, width = 8, height = 6)


# Plot simplified_clinvar --------------

# Get unique values from simplified_clinvar
unique_clinvar_scores <- unique(mobi_data$simplified_clinvar)

# Replace NA and empty values with "Not_provided"
mobi_data$simplified_clinvar[is.na(mobi_data$simplified_clinvar) | mobi_data$simplified_clinvar == ""] <- "Not_provided"

# Create the count table based on the corrected simplified_clinvar
clinvar_count <- table(mobi_data$simplified_clinvar)

# Print the count table
print(clinvar_count)

# Convert the count table into a dataframe for easier plotting
clinvar_df <- as.data.frame(clinvar_count)

# Rename the columns for easier plotting
colnames(clinvar_df) <- c("simplified_clinvar", "count")

# Define the desired order of levels
clinvar_order <- c("Benign", "Benign/Likely_benign", "Likely_benign", 
                   "Conflicting_classifications_of_pathogenicity", 
                   "Uncertain_significance", "Likely_pathogenic", 
                   "Pathogenic/Likely_pathogenic", "Pathogenic", 
                   "drug_response", "Not_provided")

# Relevel the 'simplified_clinvar' column based on the defined order
clinvar_df$simplified_clinvar <- factor(clinvar_df$simplified_clinvar, levels = clinvar_order)

plot2 <- ggplot(clinvar_df, aes(x = simplified_clinvar, y = count, fill = simplified_clinvar)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Not_provided" = "grey", 
                               "Likely_pathogenic" = "red", 
                               "Pathogenic/Likely_pathogenic" = "red", 
                               "Pathogenic" = "red", 
                               "drug_response" = "red", 
                               "Benign" = "blue", 
                               "Benign/Likely_benign" = "blue", 
                               "Likely_benign" = "blue", 
                               "Conflicting_classifications_of_pathogenicity" = "orange", 
                               "Uncertain_significance" = "orange")) +
  labs(x = "Simplified Clinvar", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  geom_text(aes(label = count), vjust = -0.3, size = 3)  # Add labels on top of bars
ggsave("/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/plots/Clinvar_histogram_all.png", plot = plot2, width = 8, height = 6)


# Filter out the "Not_provided" category for better visibility
clinvar_df_filtered <- clinvar_df[clinvar_df$simplified_clinvar != "Not_provided", ]

plot3 <- ggplot(clinvar_df_filtered, aes(x = simplified_clinvar, y = count, fill = simplified_clinvar)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Likely_pathogenic" = "red", 
                               "Pathogenic/Likely_pathogenic" = "red", 
                               "Pathogenic" = "red", 
                               "drug_response" = "red", 
                               "Benign" = "blue", 
                               "Benign/Likely_benign" = "blue", 
                               "Likely_benign" = "blue", 
                               "Conflicting_classifications_of_pathogenicity" = "orange", 
                               "Uncertain_significance" = "orange")) +
  labs(x = "Simplified Clinvar", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  geom_text(aes(label = count), vjust = -0.3, size = 3)  # Add labels on top of bars

ggsave("/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/plots/Clinvar_histogram.png", plot = plot3, width = 8, height = 6)

# Simplified groups --------

# Group the categories as per your request

# Define a new column for the groups
mobi_data$simplified_clinvar_group <- case_when(
  mobi_data$simplified_clinvar %in% c("Benign", "Benign/Likely_benign", "Likely_benign") ~ "Benign",
  mobi_data$simplified_clinvar %in% c("Conflicting_classifications_of_pathogenicity", "Uncertain_significance") ~ "Uncertain_significance",
  mobi_data$simplified_clinvar %in% c("Likely_pathogenic", "Pathogenic/Likely_pathogenic", "Pathogenic", "drug_response") ~ "Pathogenic",
  TRUE ~ "Not_provided"
)

clinvar_group_count <- table(mobi_data$simplified_clinvar_group)
print(clinvar_group_count)
clinvar_group_df <- as.data.frame(clinvar_group_count)
colnames(clinvar_group_df) <- c("simplified_clinvar_group", "count")
clinvar_group_order <- c("Benign", "Uncertain_significance", "Pathogenic", "Not_provided")
clinvar_group_df$simplified_clinvar_group <- factor(clinvar_group_df$simplified_clinvar_group, levels = clinvar_group_order)

clinvar_group_df_filtered <- clinvar_group_df[clinvar_group_df$simplified_clinvar_group != "Not_provided", ]


plot4 <- ggplot(clinvar_group_df_filtered, aes(x = simplified_clinvar_group, y = count, fill = simplified_clinvar_group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Not_provided" = "grey", 
                               "Pathogenic" = "red", 
                               "Benign" = "blue", 
                               "Uncertain_significance" = "orange")) +
  labs(x = "Simplified Clinvar Group", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  geom_text(aes(label = count), vjust = -0.3, size = 3)  # Add labels on top of bars


ggsave("/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/plots/Clinvar_histogram_simplified.png", plot = plot4, width = 8, height = 6)


# plot "gnomAD.v4.Genome."  ------------------
mobi_data$gnomAD.v4.Genome.
mobi_data$gnomAD_value <- ifelse(mobi_data$gnomAD.genome. == "No match in gnomAD genome", 0, mobi_data$gnomAD.v4.Genome.)
# gnomAD.genome.

# ggplot(mobi_data, aes(x = gnomAD_value)) +
#   geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
#   labs(title = "Distribution of gnomAD.genome.", 
#        x = "gnomAD.genome.", 
#        y = "Frequency") +
#   theme_minimal()

distrib_gnomad_freq <- ggplot(mobi_data, aes(x = gnomAD.v4.Genome.)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
  labs(title = "Distribution of gnomAD.v4.Genome", 
       x = "gnomAD.v4.Genome", 
       y = "Frequency") +
  geom_vline(xintercept = 0.0001, color = "red", linetype = "dashed", size = 1) +  # Add vertical line at x = 0.0001
  theme_minimal()

ggsave("/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/plots/gnomadfreq.png", plot = distrib_gnomad_freq, width = 8, height = 6)


# If gnomAD.v4.Genome is numeric and you want a log-transformed histogram
distrib_gnomad_freq2 <- ggplot(mobi_data, aes(x = gnomAD.v4.Genome.)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  scale_x_log10() +  # Apply log10 transformation to the x-axis
  labs(title = "Log-Distribution of gnomAD.v4.Genome", 
       x = "gnomAD.v4.Genome (log scale)", 
       y = "Frequency") +
  geom_vline(xintercept = 0.0001, color = "red", linetype = "dashed", size = 1) +  # Add vertical line at x = 0.0001
  theme_minimal()

ggsave("/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/plots/gnomadfreq_log10.png", plot = distrib_gnomad_freq2, width = 8, height = 6)



# Plot pvs1 --------
# Filter data where PVS1 is different from 0
mobi_data$LOF_binned <- cut(mobi_data$LOF_pathogenic_count,
                            breaks = c(-Inf, 0, 1, 50, 100, Inf),
                            labels = c("inf", "0", "2-50", ">50",">100"),
                            right = FALSE)

# Filter data where PVS1 is not 0
filtered_data <- mobi_data[mobi_data$PVS1 != 0, ]

# Create the plot
variant_count <- dim(filtered_data)[1]

LOF <- ggplot(filtered_data, aes(x = LOF_binned)) +
  geom_bar(fill = "skyblue", color = "black") +
  labs(title = paste("Distribution of LOF_pathogenic_count in", variant_count, "variants that are LOF"), 
       x = "LOF_pathogenic_count", 
       y = "Frequency") +
  theme_minimal()

ggsave("/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/plots/LOF.png", plot = LOF, width = 8, height = 6)



# variant type -------
table(mobi_data$variant_type)
filtered_data <- mobi_data[mobi_data$variant_type != "", ]


variant_typeplot <- ggplot(filtered_data, aes(x = variant_type)) +
  geom_bar(fill = "skyblue", color = "black") +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5) +  # Add count on top of bars
  labs(title = "Distribution of Variant Types",
       x = "Variant Type",
       y = "Count") +
  theme_minimal()

ggsave("/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/plots/variant_type.png", plot = variant_typeplot, width = 8, height = 6)



