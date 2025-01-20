
# Load necessary libraries
library(ggplot2)
library(dplyr)

# Example DataFrame (replace this with mobi_data)
# mobi_data <- data.frame('MPA score: ' = c(1.2, 3.5, 7.8, 3.3, 4.5, 2.1, 8.9, 1.3, 3.7, 2.9))
file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/mobi_data_omim_splice.txt"

mobi_data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
mobi_data$`MPA score: ` <- as.numeric(mobi_data$`MPA score: `)

# Create bins (0-10 with intervals, e.g., 10 bins)
bin_width <- 1
mobi_data$bin <- cut(mobi_data$`MPA score: `, breaks = seq(0, 10, by = bin_width), right = FALSE)

# Calculate bin counts for text labels
bin_counts <- mobi_data %>%
  group_by(bin) %>%
  summarise(count = n()) %>%
  mutate(bin_center = as.numeric(gsub("[\\(\\)]", "", bin)) + bin_width / 2)

# Plot the histogram
ggplot(mobi_data, aes(x = `MPA score: `)) +
  geom_histogram(binwidth = bin_width, fill = 'blue', color = 'black', alpha = 0.7) +
  geom_text(
    data = bin_counts,
    aes(x = bin_center, y = count, label = count),
    va = 'bottom',
    size = 5
  ) +
  labs(title = 'Histogram of MPA Score', x = 'MPA Score', y = 'Frequency') +
  theme_minimal()
