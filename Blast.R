# Initial stuff
library(Biostrings)
library(msaR)
library(caTools)
library(seqinr)
library(dplyr)

# Function to load Allergen data
prepare_data_alg <- function(stringSet){
  data <- readAAStringSet(stringSet)
  train_alg_fasta <- msaR::as.fasta(data)
  seqinr::write.fasta(train_alg_fasta, names(train_alg_fasta), "train_alg.fasta")
}

# Function to load Non allergen data
prepare_data_nlg <- function(stringSet){
  data <- readAAStringSet(stringSet)
  test_nlg_fasta <- msaR::as.fasta(data)
  seqinr::write.fasta(test_nlg_fasta, names(test_nlg_fasta), "test_nlg.fasta")
}


############################################################################
# Part 1: Prepare data for BLAST
setwd("~/Trabajo/Allergen Predictor TFM")
prepare_data_alg("allergen_data.fasta")
prepare_data_nlg("Top10_nalg.fasta")

############################################################################
# Part 2: Run BLAST
system("makeblastdb -in train_alg.fasta -dbtype 'prot'")
system("blastp -db train_alg.fasta -query test_nlg.fasta -out blast_nlg -outfmt '6 qseqid sseqid pident evalue'")


############################################################################
# Part 3: Select best homology % for each query
nlg <- data.frame(read.table("blast_nlg"))
results <- nlg %>%
  group_by(V1) %>%
  filter(V4 == min(V4)) %>%
  filter(V3 == max(V3))

colnames(results) <- c("Query ProtName", "Search ProtName", "% Homology", "e-Value")

write.table(results, "blast_results.txt", row.names = F, quote = F, col.names = T, sep = '\t')
