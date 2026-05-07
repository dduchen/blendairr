**recentAlleleClusters** - *Retrieving allele similarity clusters Zenodo archive*

Description
--------------------

A wrapper function for `zenodoArchive`, download the most recent allele similarity clusters and thresholds from the zenodo archive.
The clusters and thresholds are based on [https://yaarilab.github.io/IGHV_reference_book/](https://yaarilab.github.io/IGHV_reference_book/)
At the moment only available for human IGHV reference set.


Usage
--------------------
```
recentAlleleClusters(
doi = "10.5281/zenodo.7401189",
path,
get_file = F,
quite = F
)
```

Arguments
-------------------

doi
:   The doi for the archive to download. Default is the IGHV set.

path
:   The output folder for saving the archive files. Default is to a temporary directory.

get_file
:   Logical (FALSE by default). Do you want to return the path for the file downloaded.

quite
:   Logical (FALSE by default). Do you want to suppress informative messages




Value
-------------------

If get_file is TRUE, the function returns the path to the archive file



Examples
-------------------

```R
### Not run:
recentAlleleClusters(doi="10.5281/zenodo.7401189")
```

*Files will be downloaded to tmp directory: /tmp/Rtmpn83kLa*






