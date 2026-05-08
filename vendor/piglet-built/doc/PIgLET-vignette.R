## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  crop = knitr::hook_pdfcrop
)
suppressMessages({
  library(data.table)
  library(tigger)
  library(piglet)
  library(htmltools)
})

## ----echo=FALSE---------------------------------------------------------------
htmltools::img(src = knitr::image_uri("piglet_logo.svg"), 
               alt = 'logo', 
               style = 'position:absolute; top:0; right:0; padding:10px;border: none !important;')

## ----plot-amplicon, fig.width=11, fig.height=2, echo=FALSE, fig.cap="**V library amplicon length.** Each row is a different V coverage, S1 for full length, S2 for BIOMED-2 primers, and S3 for adaptive coverage. The colors indicates the V regions according to IMGT numbering, where dark gray represents the IMGT gaps."----

FWR1 <- "caggtgcagctggtgcagtctggggct...gaggtgaagaagcctggggcctcagtgaaggtctcctgcaaggcttct"
CDR1 <- "ggatacaccttc............accggctactat"
FWR2 <- "atgcactgggtgcgacaggcccctggacaagggcttgagtggatgggacgg"
CDR2 <- "atcaaccctaac......agtggtggcaca"
FWR3 <- "aactatgcacagaagtttcag...ggcagggtcaccagtaccagggacacgtccatcagcacagcctacatggagctgagcaggctgagatctgacgacacggtcgtgtattactgt"
CDR3 <- "gcgagaga"

regions <- c(FWR1 = FWR1, CDR1 = CDR1, FWR2 = FWR2, CDR2 = CDR2, FWR3 = FWR3, CDR3 = CDR3)
full_seq <- paste(regions, collapse = "")

# Find primer location in ungapped sequence
primer <- "GGCCTCAGTGAAGGTCTCCTGCAAG"
seq_ungapped <- gsub("[.]", "", toupper(full_seq))
loc <- unlist(aregexec(primer, seq_ungapped, max.distance = 4))
cutoff_ungapped <- loc[1] + nchar(primer)

# Map positions to regions and handle masking
seq_chars <- unlist(strsplit(toupper(full_seq), ""))
base_vals <- ifelse(seq_chars == ".", "gap", rep(names(regions), nchar(regions)))

# Calculate start index in gapped sequence (index of the cutoff_ungapped-th non-dot char)
cutoff_s2 <- which(seq_chars != ".")[cutoff_ungapped]

# Prepare Data
s1 <- base_vals
s2 <- s1; s2[seq_len(cutoff_s2 - 1)] <- NA
s3 <- s1; s3[seq_len(224)] <- NA

df <- data.frame(S1 = s1, S2 = s2, S3 = s3)
df$Position <- seq_len(nrow(df))
df_long <- tidyr::pivot_longer(df, cols = c("S1", "S2", "S3"), names_to = "Protocol", values_to = "Region")

# Factor Ordering
df_long$Protocol <- factor(df_long$Protocol, levels = c("S3", "S2", "S1"))
df_long$Region <- factor(df_long$Region, levels = c(names(regions), "gap"))

# Plot
fam_col <- c("brown4", "darkblue", "darkorchid4", "darkgreen", "firebrick", "darkorange3")
color_map <- structure(c(fam_col, "#00000099"), names = levels(df_long$Region))

ggplot2::ggplot(df_long, ggplot2::aes(x = Position, y = Protocol, fill = Region)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_manual(values = color_map, na.value = "gray90") +
  ggplot2::theme_minimal() +
  ggplot2::labs(x = "Position", y = "Coverage Protocol") +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 12, face = "bold"),
                 legend.position = "bottom",
                 panel.grid.major = ggplot2::element_blank())


## ----germline-----------------------------------------------------------------
library(piglet)
data(HVGERM)

## ----functionality------------------------------------------------------------
data(hv_functionality)

## ----lst-germlineset----------------------------------------------------------
germline <- HVGERM
## keep only functional alleles
germline <- germline[hv_functionality$allele[hv_functionality$functional=="F"]]
## keep only alleles that start from the first position of the V sequence
germline <- germline[!grepl("^[.]", germline)]
## keep only alleles that are at minimum 318 nucleotide long
germline <- germline[nchar(germline) >= 318]
## keep only localized alleles (remove NL)
germline <- germline[!grepl("NL", names(germline))]

## ----lst-germlineset-code, ref.label='lst-germlineset', anchor="block", eval=FALSE----
# germline <- HVGERM
# ## keep only functional alleles
# germline <- germline[hv_functionality$allele[hv_functionality$functional=="F"]]
# ## keep only alleles that start from the first position of the V sequence
# germline <- germline[!grepl("^[.]", germline)]
# ## keep only alleles that are at minimum 318 nucleotide long
# germline <- germline[nchar(germline) >= 318]
# ## keep only localized alleles (remove NL)
# germline <- germline[!grepl("NL", names(germline))]

## ----lst-asc, fig.cap=""------------------------------------------------------
asc <- inferAlleleClusters(
  germline_set = germline, 
  trim_3prime_side = 318, 
  mask_5prime_side = 0, 
  family_threshold = 75, 
  allele_cluster_threshold = 95)

## ----asc-plot, fig.height=18, fig.width=15, fig.cap="**Allele similarity clusters.** The out most circle is the allele names, the second layer are the ASC groups, each group is labeled and colored. The third circle is the clustering dendrogram, the branches are colored by the ASC families. The blue and orange dashed lines are the 95% and 75% similarity ASC threshold."----
plot(asc)

## ----frw1---------------------------------------------------------------------
germline_frw1 <- artificialFRW1Germline(germline, mask_primer = T)

## ----eval=FALSE---------------------------------------------------------------
# zenodo_doi <- "10.5281/zenodo.7401189"
# asc_archive <-
#   recentAlleleClusters(doi = zenodo_doi, get_file = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# allele_cluster_table <- extractASCTable(archive_file = asc_archive)

## ----echo=FALSE---------------------------------------------------------------
allele_cluster_table <- read.delim('asc_alleles_table.tsv',sep='\t')

## ----echo=FALSE---------------------------------------------------------------
asc_tables <- dplyr::left_join(
  asc$alleleClusterTable %>%
    dplyr::select(new_allele, imgt_allele) %>% dplyr::group_by(new_allele) %>%
    dplyr::summarise(imgt_allele = paste0(sort(
      unique(imgt_allele)
    ), collapse = "/")),
  allele_cluster_table %>% dplyr::select(new_allele, imgt_allele) %>% dplyr::group_by(new_allele) %>%
    dplyr::summarise(imgt_allele = paste0(sort(
      unique(imgt_allele)
    ), collapse = "/")),
  by = "new_allele",
  suffix = c(".piglet", ".zenodo")
)

## -----------------------------------------------------------------------------
# loading TIgGER AIRR-seq b cell data
data <- tigger::AIRRDb

## -----------------------------------------------------------------------------
allele_cluster_table <-
  allele_cluster_table %>% dplyr::group_by(new_allele, func_group, thresh) %>%
  dplyr::summarise(imgt_allele = paste0(sort(unique(imgt_allele)), collapse = "/"),
                   .groups = "keep")

## -----------------------------------------------------------------------------
# storing original v_call values
data$v_call_or <- data$v_call
# assigning the ASC alleles
asc_data <- assignAlleleClusters(data, allele_cluster_table)

head(asc_data[, c("v_call", "v_call_or")])

## -----------------------------------------------------------------------------
# reforming the germline set
asc_germline <- germlineASC(allele_cluster_table, germline = HVGERM)

## -----------------------------------------------------------------------------
# inferring the genotype
asc_genotype <- inferGenotypeAllele_asc(
  asc_data,
  alleleClusterTable = allele_cluster_table,
  germline_db = asc_germline,
  find_unmutated = T
)

head(asc_genotype)

## -----------------------------------------------------------------------------
# get the genotype alleles
alleles <- unlist(strsplit(asc_genotype$genotyped_imgt_alleles, ","))
# get the genes
genes <- gsub("[*][0-9]+", "", alleles)
# extract the alleles
alleles <- sapply(strsplit(alleles, "[*]"), "[[", 2)
# make sure to extract only alleles
alleles <- gsub("([0-9]+).*$", "\\1", alleles)
# create the genotype
genotype <- data.frame(alleles = alleles, gene = genes)
# plot the genotype
tigger::plotGenotype(genotype = genotype)

