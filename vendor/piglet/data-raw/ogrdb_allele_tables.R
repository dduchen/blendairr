## load the ogrdb reference set. match the VH alleles, and replace the names.
library(data.table)
library(tigger)
library(piglet)

ref_ogrdb <- readIgFasta("data-raw/HVGERM_OGRDB.fasta")
ref_piglet <- readIgFasta("data-raw/HVGERM.fasta")

ref_ogrdb_df <- data.table(allele = names(ref_ogrdb), seq = unname(ref_ogrdb))
ref_piglet_df <- data.table(allele = names(ref_piglet), seq = unname(ref_piglet))

ref_df <- merge(ref_ogrdb_df, ref_piglet_df, by = "seq", suffixes = c("_ogrdb","_piglet"), all.x = T)

ref_df_filter <- ref_df[!duplicated(ref_df$allele_ogrdb),]

ref <- setNames(ref_df_filter$seq, ref_df_filter$allele_ogrdb)

asc <- inferAlleleClusters(ref)

allele_table <- setDT(asc@alleleClusterTable)
allele_table <- allele_table[,.(iuis_allele, new_allele)]
names(allele_table) <- c("allele","asc_allele")

allele_table_piglet <- fread("data-raw/allele_threshold_table.tsv")

allele_table <- merge(allele_table, allele_table_piglet, by = "allele", suffixes = c("","_piglet"), all.x = T)
allele_table[is.na(threshold), threshold := 1e-04]

allele_table_vdjbase_piglet <- rbind(
  allele_table[,.(allele,asc_allele,threshold)],
  allele_table_piglet[!grepl("V",allele),]
)

allele_table_vdjbase_piglet[,tag:=substr(allele, 4, 4)]

fwrite(allele_table_vdjbase_piglet, "data-raw/allele_threshold_table_ogrdb.tsv", sep = "\t", quote = F, row.names = F)
writeFasta(asc@alleleClusterSet, "data-raw/HVGERM_ogrdb_asc.fasta")

ref_ogrdb_partial <- readIgFasta("data-raw/HVGERM_OGRDB_PARTIAL.fasta")

asc_partial <- inferAlleleClusters(ref_ogrdb_partial)

allele_table_partial <- setDT(asc_partial@alleleClusterTable)
allele_table_partial <- allele_table_partial[,.(iuis_allele, new_allele)]
names(allele_table_partial) <- c("allele","asc_allele")
allele_table_partial$threshold <- 1e-04
allele_table_partial$threshold <- apply(allele_table_partial, 1, function(x){
  gene <- unlist(strsplit(x[["allele"]],"[*]"))
  alleles <- unlist(strsplit(gene[2],"_"))
  gene <- gene[1]
  alleles <- paste0(gene,"*",alleles)
  thresh <- allele_table_vdjbase_piglet[allele %in% alleles, sum(threshold)]
  thresh
})

allele_table_vdjbase_piglet_partial <- rbind(
  allele_table_partial[,.(allele,asc_allele,threshold)],
  allele_table_piglet[!grepl("V",allele),]
)
allele_table_vdjbase_piglet_partial[,tag:=substr(allele, 4, 4)]

fwrite(allele_table_vdjbase_piglet_partial, "data-raw/allele_threshold_table_ogrdb_partial.tsv", sep = "\t", quote = F, row.names = F)
writeFasta(asc_partial@alleleClusterSet, "data-raw/HVGERM_ogrdb_asc_partial.fasta")
