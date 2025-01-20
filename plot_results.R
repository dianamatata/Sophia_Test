
# Load necessary libraries and data  -----------------------
library(ggplot2)
library(dplyr)

file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/mobi_data_temp.txt"
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
bin_counts$color <- ifelse(bin_counts$bin %in% c("8", "9", "10"), "red", "blue")

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

