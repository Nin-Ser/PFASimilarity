#############################################################
# Script: Similarity calculation for PFAS vs standards
# Author: Ninon Serre
# Date: 2026
# R version: 4.5.2
# Purpose:
#   - Calculate similarity scores between detected PFAS and standards by combining ECFP6 and MACCS fingerprints and fluorine content derived from SMILES representations.
#
# Usage:
#   - Place "standards.csv" and "pfas.csv" in working directory
#   - CSV format: two columns (Name, SMILES)
#   - Set working directory if needed (e.g., setwd("C:/path/to/folder"))
#   - Run script: output "similarity.csv" generated
#############################################################

#### --- Load required libraries --- ####
# Install if missing
if (!require("rcdk")) install.packages("rcdk", repos='http://cran.r-project.org') 
if (!require("fingerprint")) install.packages("fingerprint", repos='http://cran.r-project.org')
if (!require("stringr")) install.packages("stringr", repos='http://cran.r-project.org')

# Load libraries with versions used for reproducibility
library(rcdk) # version 3.8.2
library(fingerprint) # version 3.5.7
library(stringr) # version 1.6.0

#############################################################
# Functions
#############################################################

# Count F atoms in a SMILES string
count_F <- function(smiles) length(str_extract_all(smiles, "F")[[1]])

# Parse SMILES safely (returns NULL if invalid)
safe_parse <- function(smiles) {
  mol <- tryCatch(parse.smiles(smiles)[[1]], error = function(e) NULL)
  return(mol)
}

# Compute fingerprints (returns list with ECFP6 and MACCS)
compute_fps <- function(mol) {
  if (is.null(mol)) return(NULL)
  list(
    ecfp6 = get.fingerprint(mol, type = "circular", circular.type = "ECFP6"),
    maccs = get.fingerprint(mol, type = "maccs")
  )
}

# Compute weighted similarity score between two molecules
similarity_score <- function(fp1, len1, fp2, len2) {
  if (is.null(fp1) | is.null(fp2)) return(NA)
  
  tanimoto_ecfp6 <- distance(fp1$ecfp6, fp2$ecfp6, method = "tanimoto")
  tanimoto_maccs <- distance(fp1$maccs, fp2$maccs, method = "tanimoto")
  
  # Chain length similarity
  if (len1 == 0 & len2 == 0) {
    length_score <- 1
  } else {
    length_score <- 1 - abs(len1 - len2) / max(len1, len2)
  }
  
  # Weighted final score
  return((0.25 * tanimoto_ecfp6) + (0.50 * tanimoto_maccs) + (0.25 * length_score))
}

#############################################################
# Main script
#############################################################

# Read input CSVs
standards <- read.csv("standards.csv", stringsAsFactors = FALSE)
colnames(standards) <- c("Name", "SMILES")

pfas <- read.csv("pfas.csv", stringsAsFactors = FALSE)
colnames(pfas) <- c("Name", "SMILES")

# Precompute molecule objects and fingerprints for standards
standards$Mol <- lapply(standards$SMILES, safe_parse)
standards$FP  <- lapply(standards$Mol, compute_fps)
standards$Len <- sapply(standards$SMILES, count_F)

# Precompute molecule objects and fingerprints for PFAS
pfas$Mol <- lapply(pfas$SMILES, safe_parse)
pfas$FP  <- lapply(pfas$Mol, compute_fps)
pfas$Len <- sapply(pfas$SMILES, count_F)

# Track invalid PFAS
invalid_pfas_smiles <- which(is.na(pfas$SMILES) | pfas$SMILES == "")

invalid_pfas_idx <- unique(c(
  invalid_pfas_smiles,
  which(sapply(pfas$FP, is.null))
))

if(length(invalid_pfas_idx) > 0) {
  warning(length(invalid_pfas_idx), " PFAS skipped in the pfas.csv file due to invalid SMILES: ", 
          paste(pfas$Name[invalid_pfas_idx], collapse = ", "))
}

# Track invalid standards
invalid_standards_smiles <- which(is.na(standards$SMILES) | standards$SMILES == "")

invalid_standards_idx <- unique(c(
  invalid_standards_smiles,
  which(sapply(standards$FP, is.null))
))

if(length(invalid_standards_idx) > 0) {
  warning(length(invalid_standards_idx), " standards skipped in the standards.csv file due to invalid SMILES: ", 
          paste(standards$Name[invalid_standards_idx], collapse = ", "))
}

# Compute best match for each valid PFAS
results_list <- vector("list", nrow(pfas))

for (i in seq_len(nrow(pfas))) {
  
  # Skip invalid PFAS
  if (i %in% invalid_pfas_idx) next
  
  fp_pfas  <- pfas$FP[[i]]
  len_pfas <- pfas$Len[i]
  
  scores <- mapply(function(fp_std, len_std) {
    similarity_score(fp_pfas, len_pfas, fp_std, len_std)
  }, standards$FP, standards$Len)
  
  best_idx <- which.max(scores)
  
  results_list[[i]] <- data.frame(
    PFAS_Name     = pfas$Name[i],
    PFAS_SMILES   = pfas$SMILES[i],
    Best_Standard = standards$Name[best_idx],
    Similarity_Score = round(scores[best_idx], 3),
    stringsAsFactors = FALSE
  )
}

# Combine and save results
results <- do.call(rbind, results_list)
write.csv(results, "similarity.csv", row.names = FALSE)

message("Similarity calculation finished. Results saved in 'similarity.csv'")

