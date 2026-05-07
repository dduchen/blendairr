pkgname <- "piglet"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "piglet-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('piglet')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("allele_diff")
### * allele_diff

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: allele_diff
### Title: Alleles nucleotide position difference
### Aliases: allele_diff

### ** Examples

{
reference_allele = "AAGG"
sample_allele = "ATGA"

# setting position_threshold = 0 will return all differences
diff <- allele_diff(reference_allele, sample_allele)
# "A2T", "G4A"
print(diff)

# setting position_threshold = 3 will return the differences from position three onward
diff <- allele_diff(reference_allele, sample_allele, position_threshold = 3)
# "G4A"
print(diff)

# setting snps = FALSE will return the differences as indices
diff <- allele_diff(reference_allele, sample_allele, snps = FALSE)
# 2, 4
print(diff)

}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("allele_diff", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("allele_diff_indices")
### * allele_diff_indices

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: allele_diff_indices
### Title: Calculate differences between characters in columns of germs and
###   return their indices as an int vector.
### Aliases: allele_diff_indices

### ** Examples

germs = c("ATCG", "ATCC") 
X = 3 
result = allele_diff_indices(germs, X)
# 1, 2, 3



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("allele_diff_indices", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("allele_diff_indices_parallel2")
### * allele_diff_indices_parallel2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: allele_diff_indices_parallel2
### Title: Calculate SNPs or their count for each germline-input sequence
###   pair with optional parallel execution.
### Aliases: allele_diff_indices_parallel2

### ** Examples

# Example usage
germs <- c("ATCG", "ATCC")
inputs <- c("ATTG", "ATTA")
X <- 0

# Return indices of SNPs
result_indices <- allele_diff_indices_parallel2(germs, inputs, X, 
parallel = TRUE, return_count = FALSE)
print(result_indices)  # list(c(4), c(3, 4))

# Return counts of SNPs
result_counts <- allele_diff_indices_parallel2(germs, inputs, X, 
parallel = FALSE, return_count = TRUE)
print(result_counts)  # c(1, 2)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("allele_diff_indices_parallel2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("allele_diff_strings")
### * allele_diff_strings

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: allele_diff_strings
### Title: Calculate differences between characters in columns of germs and
###   return them as a string vector.
### Aliases: allele_diff_strings

### ** Examples

germs = c("ATCG", "ATCC") 
X = 3 
result = allele_diff_strings(germs, X) 
# "A2T", "T3C", "C2G"



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("allele_diff_strings", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("assignAlleleClusters")
### * assignAlleleClusters

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: assignAlleleClusters
### Title: Assign allele similarity clusters
### Aliases: assignAlleleClusters

### ** Examples



# preferably obtain the latest ASC cluster table
# asc_archive <- recentAlleleClusters(doi="10.5281/zenodo.7429773", get_file = TRUE)

# allele_cluster_table <- extractASCTable(archive_file = asc_archive)

# example allele similarity cluster table
data(allele_cluster_table)

# loading TIgGER AIRR-seq b cell data
data <- tigger::AIRRDb

asc_data <- assignAlleleClusters(data, allele_cluster_table)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("assignAlleleClusters", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("detect_communities_leiden")
### * detect_communities_leiden

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: detect_communities_leiden
### Title: Leiden community detection
### Aliases: detect_communities_leiden

### ** Examples

data(HVGERM)
d <- igDistance(HVGERM[1:10], method = "hamming")
g <- distance_to_graph(d)
comm <- detect_communities_leiden(g, resolution = 0.5)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("detect_communities_leiden", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("distance_to_graph")
### * distance_to_graph

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: distance_to_graph
### Title: Convert distance matrix to weighted graph
### Aliases: distance_to_graph

### ** Examples

data(HVGERM)
d <- igDistance(HVGERM[1:10], method = "hamming")
g <- distance_to_graph(d)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("distance_to_graph", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("extractASCTable")
### * extractASCTable

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: extractASCTable
### Title: Extracts the allele cluster table from the archive file.
### Aliases: extractASCTable

### ** Examples


## No test: 
asc_archive <- recentAlleleClusters(doi="10.5281/zenodo.7429773", get_file = TRUE)

allele_cluster_table <- extractASCTable(archive_file = asc_archive)
## End(No test)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("extractASCTable", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("germlineASC")
### * germlineASC

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: germlineASC
### Title: Converts IGHV germline set to ASC germline set.
### Aliases: germlineASC

### ** Examples


# preferably obtain the latest ASC cluster table
# asc_archive <- recentAlleleClusters(doi="10.5281/zenodo.7429773", get_file = TRUE)

# allele_cluster_table <- extractASCTable(archive_file = asc_archive)

data(HVGERM)

# example allele similarity cluster table
data(allele_cluster_table)

asc_germline <- germlineASC(allele_cluster_table, germline = HVGERM)






base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("germlineASC", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("igDistance")
### * igDistance

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: igDistance
### Title: Germline set alleles distance
### Aliases: igDistance

### ** Examples

data(HVGERM)
# Using DECIPHER method (default, for V segments)
d1 <- igDistance(HVGERM[1:10], method = "decipher")

# Using Hamming distance
d2 <- igDistance(HVGERM[1:10], method = "hamming")

# Using Levenshtein distance (good for D/J segments)
d3 <- igDistance(HVGERM[1:10], method = "lv")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("igDistance", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("inferAlleleClusters")
### * inferAlleleClusters

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: inferAlleleClusters
### Title: Allele similarity cluster
### Aliases: inferAlleleClusters

### ** Examples

# load the initial germline set
## No test: 
data(HVGERM)

germline <- HVGERM[!grepl("^[.]", HVGERM)]

# Hierarchical clustering (default)
asc <- inferAlleleClusters(germline)

# Leiden community detection
asc_leiden <- inferAlleleClusters(germline[1:50],
                                  clustering_method = "leiden",
                                  target_clusters = 10)

## plotting the clusters
plot(asc)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("inferAlleleClusters", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("inferGenotypeAllele")
### * inferGenotypeAllele

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: inferGenotypeAllele
### Title: Allele based genotype inference
### Aliases: inferGenotypeAllele

### ** Examples



# loading TIgGER AIRR-seq b cell data
data <- tigger::AIRRDb

# allele threshold table
data(allele_threshold_table)

data(HVGERM)

# inferring the genotype
genotype <- inferGenotypeAllele(
data = data,
allele_threshold_table = allele_threshold_table,
germline_db = HVGERM, find_unmutated=TRUE)

# filter alleles with z_score >= 0 

head(genotype[genotype$z_score >= 0,])




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("inferGenotypeAllele", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("inferGenotypeAllele_asc")
### * inferGenotypeAllele_asc

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: inferGenotypeAllele_asc
### Title: Allele similarity cluster based genotype inference Testing
###   function
### Aliases: inferGenotypeAllele_asc

### ** Examples



# loading TIgGER AIRR-seq b cell data
data <- tigger::AIRRDb

# preferably obtain the latest ASC cluster table
# asc_archive <- recentAlleleClusters(doi="10.5281/zenodo.7429773", get_file = TRUE)

# allele_cluster_table <- extractASCTable(archive_file = asc_archive)

# example allele similarity cluster table
data(allele_cluster_table)

data(HVGERM)

# reforming the germline set
asc_germline <- germlineASC(allele_cluster_table, germline = HVGERM)

# assigning the ASC alleles
asc_data <- assignAlleleClusters(data, allele_cluster_table)

# inferring the genotype
asc_genotype <- inferGenotypeAllele_asc(
data = asc_data,
alleleClusterTable = allele_cluster_table,
germline_db = asc_germline, find_unmutated=TRUE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("inferGenotypeAllele_asc", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("insert_gaps2_vec")
### * insert_gaps2_vec

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: insert_gaps2_vec
### Title: Insert gaps into an ungapped sequence based on a gapped
###   reference sequence.
### Aliases: insert_gaps2_vec

### ** Examples

# Example usage
gapped <- c("caggtc..aact", "caggtc---aact")
ungapped <- c("caggtcaact", "caggtcaact")

# Sequential execution
result <- insert_gaps2_vec(gapped, ungapped, parallel = FALSE)
print(result)  # "caggtc..aact", "caggtc---aact"

# Parallel execution
result_parallel <- insert_gaps2_vec(gapped, ungapped, parallel = TRUE)
print(result_parallel)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("insert_gaps2_vec", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotCommunityNetwork")
### * plotCommunityNetwork

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotCommunityNetwork
### Title: Plot community network
### Aliases: plotCommunityNetwork

### ** Examples

## No test: 
data(HVGERM)
asc <- inferAlleleClusters(HVGERM[1:30],
                           clustering_method = "leiden",
                           target_clusters = 5)
plotCommunityNetwork(asc)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotCommunityNetwork", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotSilhouetteOptimization")
### * plotSilhouetteOptimization

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotSilhouetteOptimization
### Title: Plot silhouette optimization results
### Aliases: plotSilhouetteOptimization

### ** Examples

## No test: 
data(HVGERM)
d <- igDistance(HVGERM[1:30], method = "hamming")
g <- distance_to_graph(d)
opt <- optimize_resolution(g, d, target_clusters = 5)
plotSilhouetteOptimization(opt)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotSilhouetteOptimization", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotTruncatedTree")
### * plotTruncatedTree

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotTruncatedTree
### Title: Plot truncated tree visualization
### Aliases: plotTruncatedTree

### ** Examples

## No test: 
data(HVGERM)
asc <- inferAlleleClusters(HVGERM[1:50])

# Basic truncated tree
if (requireNamespace("ggtree", quietly = TRUE)) {
  plotTruncatedTree(asc)
}
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotTruncatedTree", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("recentAlleleClusters")
### * recentAlleleClusters

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: recentAlleleClusters
### Title: Retrieving allele similarity clusters Zenodo archive
### Aliases: recentAlleleClusters

### ** Examples


## No test: 
recentAlleleClusters(doi="10.5281/zenodo.7401189")
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("recentAlleleClusters", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("zenodoArchive")
### * zenodoArchive

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: zenodoArchive
### Title: zenodoArchive
### Aliases: zenodoArchive

### ** Examples

## No test: 
  zenodo_archive <- zenodoArchive$new(
     doi = "10.5281/zenodo.7401189"
  )

  # view available version ins the archive
  archive_versions <- zenodo_archive$get_versions()

  # Getting the available files in the latest zenodo archive version
  files <- zenodo_archive$get_version_files()

  # downloading the first file from the latest archive version
  zenodo_archive$download_zenodo_files()
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("zenodoArchive", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
