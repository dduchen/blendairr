**germlineASC** - *Converts IGHV germline set to ASC germline set.*

Description
--------------------

Converts IGHV germline set to ASC germline set.


Usage
--------------------
```
germlineASC(allele_cluster_table, germline)
```

Arguments
-------------------

allele_cluster_table
:   The allele cluster table.

germline
:   An IGHV germline set with matching names to the "imgt_allele" column in the allele_cluster_table.




Value
-------------------

Returns the IGHV germline set with the ASC allele names.



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
# data(HVGERM)
# 
# asc_germline <- germlineASC(allele_cluster_table, germline = HVGERM)
```








