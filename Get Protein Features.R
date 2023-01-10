###############################################################
############                                       ############
############          Get protein features         ############
############            from FASTA file            ############
############         Design: Pablo Aliaga          ############
############                                       ############
###############################################################

# Initial stuff
library(Biostrings)
library(Peptides)
library(ggplot2)
library(optparse)
library(docstring)

# Optparse options
option_list <- list(
  make_option("--file1", action = "store", dest = "file1",
              help = "Input stringset ['name.fa'] for the positive condition."),
  make_option("--file0", action = "store", dest = "file0",
              help = "Input stringset ['name.fa'] for the negative condition."),
  make_option(c("-d", "--dframe"), action = "store", dest = "dframe",
              help = "Give a name to the data.frame with features. [Example: DF]"),
  make_option(c("-s", "--stdf"), action = "store_true", dest = "stdf",
              default = FALSE, help = "Make a stadistic data.frame.
              [default: %default]"),
  make_option(c("-v", "--violp"), action = "store_true", dest = "violinplot",
              default = FALSE, help = "Make violinplot to compare features.
              [default: %default]"),
  make_option(c("-p", "--dplot"), action = "store_true", dest = "denplot",
              default = FALSE, help = "Make density to plots to compare features.
              [default: %default")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Functions to be loaded
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

getAAComposition <- function(protSeqs) {
  #' Get aminoacid composition
  #' 
  #' Helper function to be called from within getProtFeatures to obtain
  #' the aminoacid composition
  #' 
  #' @param protSeqs List of sequences
  #' @return Aminoacid DataFrame
  AADataFrame <- data.frame(Peptides::aaComp(protSeqs))
  AADataFrame <- AADataFrame[, c(FALSE, TRUE)]
  AADataFrame <- data.frame(t(AADataFrame))
  rownames(AADataFrame) <- NULL
  colnames(AADataFrame) <- c("aaTiny_pc", "aaSmall_pc", "aaAliphatic_pc", "aaAromatic_pc",
                             "aaNonPolar_pc", "aaPolar_pc", "aaCharged_pc", "aaBasic_pc",
                             "aaAcidic_pc")
  return(AADataFrame)
}

getProtFeatures <- function(protSeqDF) {
  #' Get protein features
  #' 
  #' Obtain physicochemical features from protein DataFrame
  #' 
  #' @param protSeqDF Protein DataFrame
  #' @return DataFrame with physicochemical features and aminoacid composition
  protSeqs <- protSeqDF[, "sequence"]
  molecular_weight <- Peptides::mw(protSeqs, monoisotopic = FALSE)
  pepLength <- Peptides::lengthpep(protSeqs)
  isoelectric_point <- Peptides::pI(protSeqs, pKscale = "EMBOSS")
  instability <- Peptides::instaIndex(protSeqs)
  aliphaticIndex <- Peptides::aIndex(protSeqs)
  bomanIndex <- Peptides::boman(protSeqs)
  AADataFrame <- getAAComposition(protSeqs)
  protNames <- protSeqDF[, 1]
  protFeatsDF <- data.frame(protNames, molecular_weight, pepLength, 
                            isoelectric_point, instability, aliphaticIndex, bomanIndex, 
                            AADataFrame)
  return(protFeatsDF)
}

finalDF_seqs <- function(seqs, varName, condBin) {
  #' Final DataFrame assembly
  #' 
  #' Add Condition and binCondition columns to annotation DataFrame
  #' 
  #' @param seqs Annotation DataFrame
  #' @param varName Condition name [example: Allergen]
  #' @param condBin Condition binary name [example: 1]
  condition_seqs <- rep(varName, length(seqs[, 1]))
  binary_condition <- rep(condBin, length(seqs[, 1]))
  data.frame("condition" = condition_seqs, "binCondition" = binary_condition, seqs)
}

makeBoxplot <- function(DF, feature, yName){
  #' Make a boxplot for a feature
  #' 
  #' Helper function to be called from within getAllBoxplots. This function make
  #' a boxplot comparing conditions in DF.
  #' 
  #' @params DF DataFrame that contains the final dataset
  #' @params feature Feature to be represented
  #' @params yName Label name for Y axis
  ggplot(DF, aes_string(x = "condition", y = feature)) + geom_boxplot(
    color="blue", fill="blue", alpha=0.2, notch = TRUE, notchwidth = 0.8, outlier.colour = "red",
    outlier.fill = "red"
  ) + ylab(yName) + theme_minimal()
}

getAllBoxplots <- function(DF){
  #' Get all boxplots
  #' 
  #' Make a boxplot for each feature in the DataFrame and save them 
  #' in a PDF file
  #' 
  #' @param DF DataFrame that contains the final dataset
  col_names_to_plot <- names(DF)[-c(1:3)]
  boxplot_list <- list()
  for (col_name in col_names_to_plot) {
    boxplot_list[[col_name]] <- makeBoxplot(DF, feature = DF[[col_name]], 
                                            yName = col_name)
  }
  pdf(paste(opt$dframe, "_boxplots.pdf", sep = ""))
  for (col_name in col_names_to_plot) {
    print(boxplot_list[[col_name]])
  }
}

makeViolinplot <- function(DF, feature, yName){
  #' Make a violinplot for a feature
  #' 
  #' Helper function to be called from within getAllViolinplots. This function make
  #' a violinplot comparing conditions in DF.
  #' 
  #' @params DF DataFrame that contains the final dataset
  #' @params feature Feature to be represented
  #' @params yName Label name for Y axis
  ggplot(DF, aes_string(x = "condition", y = feature, fill = "condition")) + geom_violin(
    width = 1.4, alpha = 0.6) + scale_fill_manual(values=c("#E69F00", "#56B4E9")) + 
    geom_boxplot(color = "grey", width = 0.1, alpha = 0.2) + 
    ylab(yName) + theme_minimal() 
}

getAllViolinplots <- function(DF){
  #' Get all violinplots
  #' 
  #' Make a violinplot for each feature in the DataFrame and save them 
  #' in a PDF file
  #' 
  #' @param DF DataFrame that contains the final dataset
  col_names_to_plot <- names(DF)[-c(1:3)]
  violinplot_list <- list()
  for (col_name in col_names_to_plot) {
    if (col_name == "molecular_weight" | col_name == "pepLength") {
      violinplot_list[[col_name]] <- makeViolinplot(DF,
                                                       feature = log2(DF[[col_name]]), 
                                                       yName = col_name)
      
    } else { violinplot_list[[col_name]] <- makeViolinplot(DF, 
                                                              feature = DF[[col_name]],
                                                              yName = col_name)
    } 
  }
  pdf(paste(opt$dframe, "_violinplots.pdf", sep = ""))
  for (col_name in col_names_to_plot) {
    print(violinplot_list[[col_name]])
  }
}

makeAllTest <- function(DF) {
  #' Make statistical test
  #' 
  #' Make statistical test for each feature in the final DataFrame
  #' 
  #' @param DF DataFrame that contains the final dataset
  #' @return List with the test results
  col_names_to_test <- names(DF)[-c(1:3)]
  mwTest_list <- list()
  for (col_name in col_names_to_test) {
    mwTest_list[[col_name]] <- wilcox.test(DF[[col_name]] ~ condition, data = DF)
  }
  test_list <- list("mwTest_list" = mwTest_list)
  return(test_list)
}

presentStatistic <- function(DF, statList) {
  #' Make DataFrame with statistic results
  #' 
  #' Calculate and presents the median, its logarithm and the p-value for the Mann-Whitney test
  #' in a DataFrame
  #' 
  #' @param DF DataFrame that contains the final dataset
  #' @param statList List with the statistical test results
  #' @return DataFrame with statistic results   

  # Define feature matrix - for each feature (row) we will put the results of
  # the statistical tests
  featureM <- matrix(nrow=dim(DF)[2]-3, ncol = 4)
  colnames(featureM) <- c("Allergen median",
                          "Non allergen median", "Log2 medians",
                          "P value median")
  i <- 1
  for (protType in unique(DF$condition)) {
    featuresDF <- DF[DF$condition == protType, -c(1:3)]
    j <- 1
    medianV <- vector()
    for (feat in names(featuresDF)) {
      medianV[j] <- median(featuresDF[, feat])
      j <- j + 1
    }
    featureM[, i] <- medianV
    i <- i + 1
  }
  pValueMedian <- vector(length=dim(featuresDF)[2])
  names(pValueMedian) <- names(featuresDF)
  for (feat in names(featuresDF)) {
    pValueMedian[[feat]] <- c(statList[["mwTest_list"]][[feat]][[3]])
  }
  featureM[, "P value median"] <- pValueMedian
  log2Medians <- log2(featureM[, "Allergen median"]/featureM[, "Non allergen median"])
  featureM[, "Log2 medians"] <- log2Medians
  row.names(featureM) <- names(featuresDF)
  as.data.frame(featureM)
}

makeDensityPlot <- function(DF, feature, xName) {
  #' Make a density plot for a feature
  #' 
  #' Helper function to be called from within getAllDensityPlots. This function make
  #' a density plot comparing conditions in DF.
  #' 
  #' @params DF DataFrame that contains the final dataset
  #' @params feature Feature to be represented
  #' @params xName Label name for X axis
  ggplot(DF, aes_string(x = feature,
                        color = "condition", fill = "condition")) + geom_density(alpha=0.6) +
    xlab(xName) + theme_minimal() + scale_fill_manual(values=c("#E69F00", "#56B4E9"))
}

getAllDensityPlots <- function(DF) {
  #' Get all density plots
  #' 
  #' Make a density plot for each feature in the DataFrame and save them 
  #' in a PDF file
  #' 
  #' @param DF DataFrame that contains the final dataset
  col_names_to_plot <- names(DF[-c(1:3)])
  density_plot_list <- list()
  for (col_name in col_names_to_plot) {
    if (col_name == "molecular_weight" | col_name == "pepLength") {
      density_plot_list[[col_name]] <- makeDensityPlot(DF,
                                                       feature = log2(DF[[col_name]]), 
                                                       xName = col_name)
      
    } else { density_plot_list[[col_name]] <- makeDensityPlot(DF, 
                                                              feature = DF[[col_name]],
                                                              xName = col_name)
    } 
  }
  pdf(paste(opt$dframe, '_density_plots.pdf', sep = ""))
  for (col_name in col_names_to_plot) {
    print(density_plot_list[[col_name]])
  }
}

###############################################################
# Part 1: Load raw-data
algDF <- makeProtSeqDF(stringset = file.path(opt$file1))
nalgDF <- makeProtSeqDF(stringset = file.path(opt$file0))
###############################################################
# Part 2: Obtain feature annotation for each sequence using Peptides functions
algDF_features <- getProtFeatures(protSeqDF = algDF)
nalgDF_features <- getProtFeatures(protSeqDF = nalgDF)

###############################################################
# Part 3: Produce final data.frame
algsDF <- finalDF_seqs(seqs = algDF_features, varName = "Allergen",
                       condBin = 1)
nalgsDF <- finalDF_seqs(seqs = nalgDF_features, varName = "Non allergen",
                        condBin = 0)
DF <- rbind(algsDF, nalgsDF)

write.table(DF, paste(opt$dframe, ".txt", sep = ""), row.names = F,
            quote = F, sep = '\t')
  
###############################################################
# Part 5: Produce violin plots and density plots for each feature and statistical tests
# Make violin plots for features
if (opt$violinplot == TRUE) {
  getAllViolinplots(DF = DF)
}
# Represent densitiy for each feature and condition
if (opt$denplot == TRUE) {
  densityPlot_list <- getAllDensityPlots(DF = DF)
}
# Make Mann Whitney test for features
stat <- makeAllTest(DF = DF)
################################################################################
# Part 6: Make a statistic dataframe
if (opt$stdf == TRUE) { 
  statDF <- presentStatistic(DF = DF, statList = stat)
  write.table(statDF, paste("stat", opt$dframe, ".tsv", sep = ""), row.names = T,
              quote = F,col.names = NA, sep = '\t')
}


