# Bacterial DNA Sequence Analysis
print("My Bacterial DNA Analysis Project")

# Load the Biostrings package
library(Biostrings)

# Read bacterial DNA sequences from the FASTA file
sequences <- readDNAStringSet("data/bacterial_16S_sequences.fasta")

# Display the sequences
sequences

# Check the length of each DNA sequence
sequence_lengths <- width(sequences)

# Display the lengths
sequence_lengths

# Count the four DNA bases
base_counts <- letterFrequency(sequences, letters = c("A", "T", "G", "C"))

# Display the base counts
base_counts

# Calculate GC content percentage
gc_content <- (base_counts[, "G"] + base_counts[, "C"]) / sequence_lengths * 100

# Display GC content
gc_content

# Create a results table
results <- data.frame(
  Species = c("Escherichia coli",
              "Pseudomonas aeruginosa",
              "Klebsiella pneumoniae",
              "Bacillus subtilis",
              "Staphylococcus aureus",
              "Salmonella enterica"),
  Length = sequence_lengths,
  A = base_counts[, "A"],
  T = base_counts[, "T"],
  G = base_counts[, "G"],
  C = base_counts[, "C"],
  GC_Content = gc_content
)

# Display the results table
results

# Save the results table as a CSV file
write.csv(results, "results/nucleotide_analysis.csv", row.names = FALSE)

# Load ggplot2 package
library(ggplot2)

# Create a GC content comparison graph
ggplot(results, aes(x = Species, y = GC_Content)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(
    title = "GC Content of Six Bacterial Species",
    x = "Bacterial Species",
    y = "GC Content (%)"
  ) +
  theme_minimal()

# Save the GC content graph
ggsave("results/gc_content_comparison.png",
       width = 10,
       height = 6,
       dpi = 300)

# Convert nucleotide counts into a data frame
nucleotide_data <- data.frame(
  Species = rep(results$Species, 4),
  Base = rep(c("A", "T", "G", "C"), each = 6),
  Count = c(base_counts[, "A"],
            base_counts[, "T"],
            base_counts[, "G"],
            base_counts[, "C"])
)

# Display nucleotide data
nucleotide_data

# Create nucleotide composition graph
ggplot(nucleotide_data, aes(x = Species, y = Count, fill = Base)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Nucleotide Composition of Six Bacterial Species",
    x = "Bacterial Species",
    y = "Number of Nucleotides"
  ) +
  theme_minimal()

# Save nucleotide composition graph
ggsave("results/nucleotide_composition.png",
       width = 12,
       height = 7,
       dpi = 300)

# Load DECIPHER package
library(DECIPHER)

# Align the six bacterial DNA sequences
aligned_sequences <- AlignSeqs(sequences)

# Display aligned sequences
aligned_sequences

# Save aligned sequences in FASTA format
writeXStringSet(aligned_sequences,
                "results/aligned_sequences.fasta",
                format = "fasta")

# Calculate pairwise distances between bacterial sequences
distance_matrix <- DistanceMatrix(aligned_sequences)

# Convert distances into percentage similarity
similarity_matrix <- (1 - distance_matrix) * 100

# Display the similarity matrix
round(similarity_matrix, 2)


# Give the sequences shorter names for the graph
rownames(similarity_matrix) <- c(
  "E. coli",
  "P. aeruginosa",
  "K. pneumoniae",
  "B. subtilis",
  "S. aureus",
  "S. enterica"
)

colnames(similarity_matrix) <- c(
  "E. coli",
  "P. aeruginosa",
  "K. pneumoniae",
  "B. subtilis",
  "S. aureus",
  "S. enterica"
)

# Convert the matrix into a format for ggplot2
similarity_data <- as.data.frame(as.table(similarity_matrix))

# Create the heatmap
ggplot(similarity_data, aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f%%", Freq)), size = 3) +
  scale_fill_gradient(low = "white", high = "steelblue",
                      limits = c(0, 100)) +
  labs(
    title = "16S rRNA Sequence Similarity Heatmap",
    x = "Bacterial Species",
    y = "Bacterial Species",
    fill = "Similarity (%)"
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the heatmap
ggsave("results/sequence_similarity_heatmap.png",
       width = 10, height = 8, dpi = 300)

# Save the pairwise distance matrix
write.csv(
  distance_matrix,
  "results/sequence_distance_matrix.csv"
)

# Save the similarity matrix
write.csv(
  similarity_matrix,
  "results/sequence_similarity_matrix.csv"
)

