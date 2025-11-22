# read the file
file_path_1 <- commandArgs()[6]
lines_1 <- readLines(file_path_1)
# Extract sequence
sequences_1 <- sub("^\\s*\\d+\\s*", "", lines_1)
sequences_1 <- sequences_1[grepl("^[ACGT]+$", sequences_1)]

# Sequence grouping into pairs
seq_pairs_1 <- split(sequences_1, ceiling(seq_along(sequences_1)/2))

# initialize
mutations_1 <- numeric(length(seq_pairs_1))

# Calculate the number of mutation sites in each pair of sequences
for (i in seq_along(seq_pairs_1)) {
  seq1 <- unlist(strsplit(seq_pairs_1[[i]][1], ""))
  seq2 <- unlist(strsplit(seq_pairs_1[[i]][2], ""))
  mutations_1[i] <- sum(seq1 != seq2)
}
mutation_rate_corrected <- (-3/4)*log(1-((4*mutations_1)/(3*1000)))

# Calculate the Ks distribution-mean and Ks distribution-variance
mean <- mean(mutation_rate_corrected)
variance <- var(mutation_rate_corrected)
mean <- sprintf("%.10f",mean)
variance <- sprintf("%.10f",variance)
cat(mean,"\t",variance,sep="")

