# Load required R package
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(dplyr)

#Setting the home directory
setwd("D:/KCL2024/Courses/7BBG1002_Cloud_computing/Project")

# 1. Read the file
file_path <- "data/processed/counts_folder/counts_alignment_SRR1166447.txt"  # Replace with your file path
data <- read.table(file_path, header = TRUE, sep = "\t")

# View the structure of the data
print(head(data))

# 2. Extract the required columns: gene length (Length) and read counts (last column)
gene_length <- data$Length
read_counts <- data[, ncol(data)]

# 3. Calculate the total read counts (sequencing depth)
library_size <- sum(read_counts)

# 4. Calculate RPKM
calc_rpkm <- function(counts, gene_length, library_size) {
  rpkm <- (counts / gene_length) / (library_size / 1e6)
  return(rpkm)
}
rpkm <- calc_rpkm(read_counts, gene_length, library_size)

# 5. Add RPKM values to the data frame
data$rpkm <- rpkm

# 6. Perform "By totals" normalization
normalize_by_totals <- function(rpkm) {
  normalized <- rpkm / sum(rpkm) * 1e6  # Scale to millions
  return(normalized)
}
data$normalized_rpkm <- normalize_by_totals(rpkm)

# 7. View the results
print(head(data))

# 8. Save the results to a CSV file
write.csv(data, "data/processed/Normalized_data/normalized_SRR1166447.csv", row.names = FALSE)
