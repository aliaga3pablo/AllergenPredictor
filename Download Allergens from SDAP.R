###############################################################
############                                       ############
############          Download data from           ############
############               SDAP web                ############
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
              help = "Input html ['http://example.com'] for the SDAP web table."),
  make_option(c("-o", "--output"), action = "store", dest = "output", default = "sdap",
              help = "Give a name to the txt file with proteins. [Example: sdap].
              [default: %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))

###############################################################
# Download data from SDAP
# Read html table link
html <- read_html(opt$input)

# Get the names for the different allergens from the table
population <- html %>% html_nodes(xpath='//*[@width="90%"]') %>%
  html_table %>% extract2(1)
alrgNames <- unique(population[population[,2] != "",2])

# Remove "header" row
alrgNames <- alrgNames[-1,]

# Get the links from the allergens from the table and them full links as input 
# for the next part
links <- paste0("https://fermi.utmb.edu/cgi-bin/SDAP/", 
                 html %>% html_nodes(xpath = '//td/a') %>% html_attr("href"))

# Remove first items (not useful) and those for allergen.org. Remove final white space too
links <- links[-c(1:8)]
links <- links[! links == "https://fermi.utmb.edu/cgi-bin/SDAP/http://www.allergen.org/"]

links <- sapply(links, function(u) gsub(" ", "", u))
links <- unname(links)


# Follow these links to get the IDs, first it's necessary to change encoding to UTF-8
alr_up_list <- c()

cat("Searching Database...")
for(i in 1:length(links)){
  bytes <- readLines(links[i])
  utf8 <- iconv(bytes, from = "windows-1252", to = "UTF-8")
  aux_html <- read_html(charToRaw(paste(utf8, collapse = "\n")), encoding = "UTF8")
  aux_table <- aux_html %>% html_nodes(xpath = '//*[@border="2"]') %>% html_table()
  aux_table2 <- extract2(aux_table, 2) %>% setNames(., .[1,]) %>% extract(-1,) 
  if ("Link to Source" %in% names(aux_table2)){
    alr_up_list <- c(alr_up_list, aux_table2[1,2])
  } else {
    aux_table2 <- extract2(aux_table, 3) %>% setNames(., .[1,]) %>% extract(-1,)
    if ("Link to Source" %in% names(aux_table2)){
      alr_up_list <- c(alr_up_list, aux_table2[1,2])
    } else {
      aux_table2 <- extract2(aux_table, 4) %>% setNames(., .[1,]) %>% extract(-1,)
      alr_up_list <- c(alr_up_list, aux_table2[1,2])
    }
  }
}

# Clean up Ids
alr_up <- unname(unlist(alr_up_list))
alr_up <- alr_up[! is.na(alr_up)]
alr_up <- alr_up[! alr_up == ""]

write.table(alr_up, paste(opt$output, ".txt", sep = ""),quote = F, row.names = F, col.names = F)

cat("\nDone!")

cat("\n\n###########################################################\n\n")
cat("Use UniProt Retrieve/ID mapping tool link: 'https://www.uniprot.org/id-mapping/' 
to convert RefSeq_Protein, GI_number, UniportKB_AC-ID, EMBL-GenBank-DDBJ and
EMBL-GenBank-DDBJ_CDS to UniprotKB accesion ID.\n\n")
cat("Merge all this data and form an unique table.")
cat("\n\n###########################################################\n\n")

