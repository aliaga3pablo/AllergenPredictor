###############################################################
############                                       ############
############          Download data from           ############
############            allergen.org               ############
############         Design: Pablo Aliaga          ############
############                                       ############
###############################################################

# Initial stuff
library(rvest)
library(magrittr)
library(Biostrings)
library(seqRFLP)
library(xml2)
library(stringi)
library(stringr)
library(optparse)

# Optparse options
option_list <- list(
  make_option(c("-i","--input"), action = "store", dest = "input",
              help = "Input html ['http://example.com'] for the allergen.org table."),
  make_option(c("-o", "--output"), action = "store", dest = "output", default = "allergen_org",
              help = "Give a name to the txt file with proteins. [Example: allergen_org].
              [default: %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))

###############################################################
# Read html table link
html <- read_html(opt$input)

# Get the names for the different allergens from the table
population <- html %>% html_nodes(xpath='//*[@id="example"]') %>%
  html_table %>% extract2(1)
alrgNames <- unique(population[population[,2] != "",2])
  
# Get the links from the allergens from the table and them full links as input 
# for the next part
links <- paste0("http://www.allergen.org/", html %>% html_nodes(xpath = "//td/a") %>%
                  html_attr("href"))

cat("Searching Database...")
# Follow these links to get the Uniprot IDs
alr_up_list <- sapply(links,
                      function(u) u %>% read_html %>%
                        html_nodes(xpath = '//*[@id="isotable"]') %>% html_table() %>%
                        extract2(1) %>% extract("UniProt")
)

# Remove 1 repeated allergen in alr_up_list
alr_up_list <- alr_up_list[names(alr_up_list) %in%
                             "http://www.allergen.org/viewallergen.php?aid=963.UniProt" == FALSE]
  
# Convert to vector and rename
alr_up <- setNames(unlist(alr_up_list, use.names=F), 
                   rep(unlist(alrgNames, use.names = F), 
                       times=lengths(alr_up_list, use.names = F)))

# Clean up Ids
alr_up <- alr_up[! is.na(alr_up)]
alr_up <- alr_up[! alr_up == ""]

write.table(alr_up, paste(opt$output, ".txt", sep = ""), quote = F, row.names = F, col.names = F)

cat("\nDone!\n")


