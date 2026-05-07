**assignAlleleClusters** - *Assign allele similarity clusters*

Description
--------------------

`assignAlleleClusters` uses the allele clusters annotation to change the preliminary allele
assignments to the new annotations before inferring a genotype.


Usage
--------------------
```
assignAlleleClusters(data, alleleClusterTable, v_call = "v_call")
```

Arguments
-------------------

data
:   data.frame in AIRR format, containing V allele calls from a single subject and the sample IMGT-gapped V(D)J sequences under seq.

alleleClusterTable
:   A data.frame of the allele clusters new annotations relative to the original reference set. See details.

v_call
:   name of the V allele call column. Default is `v_call`




Value
-------------------

A modified input `data.frame` with the new assigned



Examples
-------------------

```R
### Not run:
asc_archive <- recentAlleleClusters(doi="10.5281/zenodo.7401239", get_file = TRUE)

```

*Files will be downloaded to tmp directory: /tmp/Rtmpn83kLa**doi supplied is not an 'all versions doi'
for viewing all of the archive records change the doi to:10.5281/zenodo.7401189*
```R
# 
# allele_cluster_table <- extractASCTable(archive_file = asc_archive)
# 
# # loading TIgGER AIRR-seq b cell data
# data <- tigger::AIRRDb
# 
# asc_data <- assignAlleleClusters(data, allele_cluster_table)
```








