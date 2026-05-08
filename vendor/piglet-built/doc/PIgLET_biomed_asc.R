## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
suppressMessages({
  library(data.table)
  library(tigger)
  library(piglet)
})

## ----eval=FALSE---------------------------------------------------------------
# url <- "https://bitbucket.org/yaarilab/piglet/raw/70b7d4491e25e7197e2a94bd890ce5b6e3b506a8/data-raw/HVGERM_OGRDB.fasta"
# tmp_dest_file <- file.path(tempdir(), "HVGERM_OGRDB.fasta")
# download.file(url, tmp_dest_file, mode = "wb")
# 
# ref_ogrdb <- readIgFasta(tmp_dest_file)
# 
# ref_ogrdb_frw1 <- piglet::artificialFRW1Germline(ref_ogrdb)
# 

## ----eval=FALSE---------------------------------------------------------------
# asc_frw1 <- inferAlleleClusters(ref_ogrdb_frw1)
# allele_table_frw1 <- setDT(asc_frw1@alleleClusterTable)[, .(imgt_allele, new_allele)]
# setnames(allele_table_frw1, c("allele", "asc_allele"))

## ----eval=FALSE---------------------------------------------------------------
# allele_table_piglet <- fread("https://bitbucket.org/yaarilab/piglet/raw/70b7d4491e25e7197e2a94bd890ce5b6e3b506a8/data-raw/allele_threshold_table.tsv")
# 
# allele_table_frw1$threshold <- 1e-04
# allele_table_frw1$threshold <- apply(allele_table_frw1, 1, function(x){
#   gene <- unlist(strsplit(x[["allele"]],"[*]"))
#   alleles <- unlist(strsplit(gene[2],"_"))
#   gene <- gene[1]
#   alleles <- paste0(gene,"*",alleles)
#   thresh <- allele_table_piglet[allele %in% alleles, sum(threshold)]
#   thresh
# })
# 
# allele_table_frw1 <- rbind(
#   allele_table_frw1[,.(allele,asc_allele,threshold)],
#   allele_table_piglet[!grepl("V",allele),]
# )
# allele_table_frw1[,tag:=substr(allele, 4, 4)]
# 

## ----eval=FALSE---------------------------------------------------------------
# url <- "https://bitbucket.org/yaarilab/piglet/raw/70b7d4491e25e7197e2a94bd890ce5b6e3b506a8/data-raw/HVGERM_ogrdb_asc_partial.fasta"
# tmp_dest_file <- file.path(tempdir(), "HVGERM_ogrdb_asc_partial.fasta")
# download.file(url, tmp_dest_file, mode = "wb")
# 
# asc_germline <- readIgFasta(tmp_dest_file)
# 
# allele_table <- fread("https://bitbucket.org/yaarilab/piglet/raw/70b7d4491e25e7197e2a94bd890ce5b6e3b506a8/data-raw/allele_threshold_table_ogrdb_partial.tsv")

## ----eval=FALSE---------------------------------------------------------------
# data <- tigger::AIRRDb
# data$v_call_or <- data$v_call
# 
# allele_table_split <- allele_table[, {
#   parts <- strsplit(allele, "\\*")[[1]]
#   gene <- parts[1]
#   alleles <- strsplit(parts[2], "_")[[1]]
#   expanded <- paste0(gene, "*", alleles)
#   .(allele = expanded, asc_allele, threshold, tag)
# }, by = .I]
# 
# allele_table_split[, I := NULL]
# 
# asc_data <- assignAlleleClusters(data, allele_table_split, v_call = "v_call", from_col = "allele", to_col = "asc_allele")
# 

## ----eval=FALSE---------------------------------------------------------------
# 
# # using the asc annotations
# 
# asc_genotype <- inferGenotypeAllele(
#   asc_data,
#   allele_threshold_table = allele_table,
#   call = "v_call", # change to the column call you want to genotype
#   asc_annotation = TRUE, # if you use iuis names then set to FALSE
#   single_assignment = TRUE, # if you want to use the single assignment algorithm
#   find_unmutated = FALSE # change to TRUE to filter mutated reads
#   # germline_db = asc_germline # Uncomment if you want to filter mutated reads
# )
# 
# # using the biomed annotations, make sure to convert the v_call to the collapsed biomed annotations
# allele_table_biomed <- allele_table
# allele_table_biomed[, asc_allele := allele]
# 
# allele_table_split <- allele_table_biomed[, {
#   parts <- strsplit(allele, "\\*")[[1]]
#   gene <- parts[1]
#   alleles <- strsplit(parts[2], "_")[[1]]
#   expanded <- paste0(gene, "*", alleles)
#   .(allele = expanded, asc_allele, threshold, tag)
# }, by = .I]
# 
# allele_table_split[, I := NULL]
# 
# biomed_data <- assignAlleleClusters(data, allele_table_split, v_call = "v_call", from_col = "allele", to_col = "asc_allele")
# 
# asc_genotype_biomed <- inferGenotypeAllele(
#   biomed_data,
#   allele_threshold_table = allele_table,
#   call = "v_call", # change to the column call you want to genotype
#   asc_annotation = FALSE, # if you use iuis names then set to FALSE
#   single_assignment = TRUE
# )
# 

