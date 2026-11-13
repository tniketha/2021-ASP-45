# Run from the repository root with Rscript code/01_bacterial_sequence_analysis.R.
library(Biostrings)
library(DECIPHER)
library(ggplot2)

input_file <- "data/bacterial_16S_sequences.fasta"
output_dir <- "results"

# Match labels by versioned accession and use a consistent order for the plots.
metadata <- data.frame(
  Accession = c("NR_024570.1", "NR_117678.1", "NR_117686.1",
                "NR_102783.2", "NR_037007.2", "NR_074910.1"),
  Species = c("Escherichia coli", "Pseudomonas aeruginosa", "Klebsiella pneumoniae",
              "Bacillus subtilis", "Staphylococcus aureus", "Salmonella enterica"),
  Label = c("E. coli", "P. aeruginosa", "K. pneumoniae",
            "B. subtilis", "S. aureus", "S. enterica")
)

# Read as text first so invalid symbols cannot be silently discarded on import.
records <- readBStringSet(input_file, format = "fasta")
headers <- names(records)
accessions <- sub("[[:space:]].*$", "", headers)
if (anyDuplicated(accessions) || !setequal(accessions, metadata$Accession)) {
  stop("The FASTA must contain exactly one record for each expected accession.")
}
records <- records[match(metadata$Accession, accessions)]
dna <- toupper(as.character(records))
if (any(!grepl("^[ACGTRYSWKMBDHVN]+$", dna))) {
  stop("Each unaligned FASTA record must contain only IUPAC DNA bases.")
}
sequences <- DNAStringSet(dna)
names(sequences) <- names(records)
coverage <- ifelse(grepl("partial sequence$", names(sequences)), "Partial", "Complete")

bases <- c("A", "T", "G", "C")
base_counts <- letterFrequency(sequences, letters = bases)
sequence_lengths <- width(sequences)
called_bases <- rowSums(base_counts)
ambiguous_bases <- sequence_lengths - called_bases
if (any(called_bases == 0L)) {
  stop("Each sequence must contain at least one unambiguous base.")
}

# Keep the original total-length GC denominator and report ambiguous calls.
results <- data.frame(
  Accession = metadata$Accession,
  Species = metadata$Species,
  Coverage = coverage,
  Length = sequence_lengths,
  base_counts,
  Ambiguous_Bases = ambiguous_bases,
  GC_Content = 100 * (base_counts[, "G"] + base_counts[, "C"]) / sequence_lengths
)

dir.create(output_dir, showWarnings = FALSE)
write.csv(results, file.path(output_dir, "nucleotide_analysis.csv"), row.names = FALSE)

# Use proportions for the composition plot to compare records of unequal length.
plot_labels <- factor(metadata$Label, levels = metadata$Label)
composition <- data.frame(
  Species = rep(plot_labels, times = length(bases)),
  Base = factor(rep(bases, each = nrow(metadata)), levels = bases),
  Percent = as.vector(100 * base_counts / sequence_lengths)
)

composition_plot <- ggplot(composition, aes(Species, Percent, fill = Base)) +
  geom_col(position = "dodge") +
  scale_fill_brewer(palette = "Dark2") +
  labs(title = "Nucleotide composition of bacterial 16S rRNA records",
       x = NULL, y = "Percentage of sequence positions (%)") +
  theme_minimal(base_size = 12)

gc_plot <- ggplot(results, aes(plot_labels, GC_Content)) +
  geom_col(fill = "#297D85", width = 0.65) +
  geom_text(aes(label = sprintf("%.2f", GC_Content)), vjust = -0.5, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "GC content of bacterial 16S rRNA records",
       x = NULL, y = "GC content (%)") +
  theme_minimal(base_size = 12)

ggsave(file.path(output_dir, "nucleotide_composition.png"), composition_plot,
       width = 10, height = 5, dpi = 300)
ggsave(file.path(output_dir, "gc_content_comparison.png"), gc_plot,
       width = 10, height = 5, dpi = 300)

# Align the records and compare overlapping positions without penalizing gaps.
aligned_sequences <- AlignSeqs(sequences, processors = 1, verbose = FALSE)
writeXStringSet(aligned_sequences, file.path(output_dir, "aligned_sequences.fasta"))

distance_matrix <- DistanceMatrix(
  aligned_sequences, method = "overlap", type = "matrix",
  includeTerminalGaps = FALSE, penalizeGapLetterMatches = FALSE,
  minCoverage = 0, correction = NA, processors = 1, verbose = FALSE
)
if (any(!is.finite(distance_matrix))) {
  stop("Some pairs have no comparable alignment positions.")
}

aligned_accessions <- sub("[[:space:]].*$", "", names(aligned_sequences))
labels <- metadata$Label[match(aligned_accessions, metadata$Accession)]
dimnames(distance_matrix) <- list(labels, labels)
similarity_matrix <- (1 - distance_matrix) * 100

write.csv(distance_matrix, file.path(output_dir, "sequence_distance_matrix.csv"))
write.csv(similarity_matrix, file.path(output_dir, "sequence_similarity_matrix.csv"))

similarity_data <- as.data.frame(as.table(similarity_matrix),
                                 responseName = "Similarity")
similarity_plot <- ggplot(similarity_data, aes(Var2, Var1, fill = Similarity)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", Similarity)), size = 3.5) +
  scale_fill_gradient(low = "#F2F6F9", high = "#71A9CC",
                      name = "Similarity (%)") +
  scale_y_discrete(limits = rev(labels)) +
  labs(title = "Pairwise 16S rRNA sequence similarity", x = NULL, y = NULL) +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1))

ggsave(file.path(output_dir, "sequence_similarity_heatmap.png"), similarity_plot,
       width = 8, height = 6, dpi = 300)

# Use a three-group cut to summarize this exploratory average-linkage tree.
cluster_result <- hclust(as.dist(distance_matrix), method = "average")
group_count <- 3L
cluster_groups <- cutree(cluster_result, k = group_count)
cluster_results <- data.frame(
  Accession = metadata$Accession,
  Species = metadata$Species,
  Cluster = unname(cluster_groups[metadata$Label])
)
write.csv(cluster_results, file.path(output_dir, "bacterial_cluster_groups.csv"),
          row.names = FALSE)

png(file.path(output_dir, "bacterial_clustering.png"),
    width = 2400, height = 1500, res = 300)
par(mar = c(5, 5, 4, 2))
plot(cluster_result, hang = -1,
     main = "Average-linkage clustering of 16S rRNA records",
     xlab = "", sub = "Boxes indicate the three-group cut",
     ylab = "Pairwise dissimilarity at merge")
par(xpd = TRUE)
rect.hclust(cluster_result, k = group_count, border = "#297D85")
invisible(dev.off())

# Recalculate distances using the same unambiguous alignment columns for every pair.
alignment <- as.matrix(aligned_sequences)
shared_columns <- colSums(matrix(alignment %in% bases, nrow = nrow(alignment))) ==
  nrow(alignment)
if (!any(shared_columns)) {
  stop("The alignment contains no columns with a called base in every record.")
}
shared_alignment <- alignment[, shared_columns, drop = FALSE]
pairs <- combn(seq_along(labels), 2L)
first <- pairs[1L, ]
second <- pairs[2L, ]
shared_distances <- rowMeans(shared_alignment[first, , drop = FALSE] !=
                               shared_alignment[second, , drop = FALSE])
shared_matrix <- matrix(0, length(labels), length(labels),
                        dimnames = list(labels, labels))
shared_matrix[cbind(first, second)] <- shared_distances
shared_matrix[cbind(second, first)] <- shared_distances
shared_clusters <- hclust(as.dist(shared_matrix), method = "average")
shared_groups <- cutree(shared_clusters, k = group_count)

comparison <- data.frame(
  Record_1 = labels[first],
  Record_2 = labels[second],
  Shared_Columns = ncol(shared_alignment),
  Overlap_Similarity_Percent = similarity_matrix[cbind(first, second)],
  Shared_Similarity_Percent = 100 * (1 - shared_distances),
  Same_Overlap_Group = unname(cluster_groups[labels[first]] ==
                               cluster_groups[labels[second]]),
  Same_Shared_Group = unname(shared_groups[labels[first]] ==
                              shared_groups[labels[second]])
)
write.csv(comparison, file.path(output_dir, "shared_position_comparison.csv"),
          row.names = FALSE)
