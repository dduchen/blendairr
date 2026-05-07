### prepare the datasets for the package.

## germline
HVGERM <- tigger::readIgFasta("data-raw/HVGERM.fasta")

usethis::use_data(HVGERM, overwrite = TRUE)

## functionality
hv_functionality <- data.table::fread("data-raw/ighv_functionality.tsv")

usethis::use_data(hv_functionality, overwrite = TRUE)

## allele cluster table
allele_threshold_table <- data.table::fread("data-raw/allele_threshold_table.tsv")
allele_threshold_table[,tag:=substr(allele, 4, 4)]
usethis::use_data(allele_threshold_table, overwrite = TRUE)
