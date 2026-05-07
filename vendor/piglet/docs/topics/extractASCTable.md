**extractASCTable** - *Extracts the allele cluster table from the archive file.*

Description
--------------------

Extracts the allele cluster table from the archive file.


Usage
--------------------
```
extractASCTable(archive_file = NULL)
```

Arguments
-------------------

archive_file
:   A path to the asc archive file. Default is null. (see details)




Value
-------------------

Returns the allele cluster table.

The table columns:
`new_allele` - the ASC given allele name
`func_group` - the ASC cluster number
`imgt_allele` - the original IUIS/IMGT allele name
`thresh` - the allele threshold for ASC-based genotype inference
`amplicon_length` - is the original length of the reference set.


Details
-------------------

For downloading the latest archive file with the updated allele cluster table, use the function `recentAlleleClusters`.



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
```








