###############################################################
############                                       ############
############          Prepare non allergen         ############
############                dataset                ############
############         Design: Pablo Aliaga          ############
############                                       ############
###############################################################

# Initial stuff
library(Biostrings)
library(seqinr)
library(optparse)
library(docstring)

# Optparse options
option_list <- list(
  make_option(c("-i","--input"), action = "store", dest = "input",
              help = "FASTA file with filtered Uniprot search."),
  make_option(c("-n", "--number"), action = "store", dest = "number", default = 21900,
              help = "Number of proteins to be selected"),
  make_option(c("-o", "--output"), action = "store", dest = "output", default = "non-allergen",
              help = "Give a name to the FASTA file with non allergens.
              [default: %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Function to be loaded
makeProtSeqDF <- function(stringset) {
  #' Make a DataFrame from FASTA
  #' 
  #' @param stringset The FASTA file to be converted
  #' @return DataFrame with protein names and sequences
  first_DF <- readAAStringSet(stringset)
  seq_name <- names(first_DF)
  sequence <- paste(first_DF)
  protSeqDF <- data.frame(seq_name, sequence, stringsAsFactors = FALSE)
  return(protSeqDF)
}

###############################################################
# Load Uniprot data
uniprot_total <- makeProtSeqDF(opt$input)

# Set random seed and take a sample
set.seed(150)
uniprot_filtered <- uniprot_total[sample(nrow(uniprot_total), opt$number), ]

# Save new fasta
write.fasta(as.list(uniprot_filtered[,2]), uniprot_filtered[,1], file.out=paste(opt$output,
                                                                                ".fasta", sep=""))





