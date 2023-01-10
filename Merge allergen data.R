###############################################################
############                                       ############
############          Merge allergen data          ############
############                                       ############
############         Design: Pablo Aliaga          ############
############                                       ############
###############################################################

# Initial stuff
library(readxl)
library(optparse)

# Optparse options
option_list <- list(
  make_option("--input1", action = "store", dest = "input1",
              help = "Input Allergen Online protein excel table."),
  make_option("--input2", action = "store", dest = "input2",
              help = "Input SDAP protein excel table."),
  make_option("--input3", action = "store", dest = "input3",
              help = "Input allergen.org protein '.txt' list."),
  make_option(c("-o", "--output"), action = "store", dest = "output", default = "allergens",
              help = "Name for the '.txt' file containing the total list of allergens.
              [Default: %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))

########################################################
# Load data from Allergen Online excel
allergen_online <- read_excel(opt$input1)
allergen_online <- unlist(allergen_online[,2])

# Load data from SDAP excel
sdap <- read_excel(opt$input2)
sdap <- unlist(sdap[,2])

# Load data from allergen.org txt file
allergen_org <- read.table(opt$input3)
allergen_org <- unname(unlist(allergen_org))

# Merge allergens
alr_up_total <- c(allergen_org, sdap, allergen_online)
alr_filtered <- unique(alr_up_total)

# Write list with Uniprot IDs
write.table(alr_filtered, paste(opt$output, ".txt", sep=""), quote = F, row.names = F, col.names = F)

cat("\n###########################################################\n\n")
cat("Use UniProt Retrieve/ID mapping tool link: 'https://www.uniprot.org/id-mapping/' 
to convert UniprotKB IDs to FASTA file.\n\n")
cat("###########################################################\n\n")
