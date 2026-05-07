**zenodoArchive** - *zenodoArchive*

Description
--------------------

zenodoArchive

zenodoArchive






Format
-------------------

`[R6Class](http://www.rdocumentation.org/packages/R6/topics/R6Class)` object.


Value
-------------------

Object of `[R6Class](http://www.rdocumentation.org/packages/R6/topics/R6Class)` for modelling an zenodoArchive for ASC cluster files


Public fields
-------------------



`doi`
:   zenodoArchive doi, NULL is not supplied

`all_versions`
:   zenodoArchive if to return all versions, `true` when not specified

`sort`
:   zenodoArchive how to sort the records, `mostrecent` when not specified

`page`
:   zenodoArchive which page to pull in query, `1` when not specified

`size`
:   zenodoArchive how many records per page,  `20` when not specified

`zenodoVersions`
:   zenodoArchive doi available version, a storing variable.

`zenodoQuery`
:   zenodoArchive doi version query, a storing variable.

`download_file`
:   zenodoArchive doi downloads files, a storing variable.

`download_url`
:   zenodoArchive doi downloads urls, a storing variable.




Methods
-------------------

Public methods

+  [`zenodoArchive$new()`](#method-zenodoArchive-new)
+  [`zenodoArchive$clean_doi()`](#method-zenodoArchive-clean_doi)
+  [`zenodoArchive$zenodo_query()`](#method-zenodoArchive-zenodo_query)
+  [`zenodoArchive$get_versions()`](#method-zenodoArchive-get_versions)
+  [`zenodoArchive$get_version_files()`](#method-zenodoArchive-get_version_files)
+  [`zenodoArchive$download_zenodo_files()`](#method-zenodoArchive-download_zenodo_files)
+  [`zenodoArchive$clone()`](#method-zenodoArchive-clone)





Method `new()`
initializes the zenodoArchive
Usage
zenodoArchive$new(
doi,
page = 1,
size = 20,
all_versions = "true",
sort = "mostrecent"
)


Arguments


`doi`
:   A zenodo doi. To retrieve all records supply a concept doi (a generic doi common to all versions).

`page`
:   Which page to query. Default is 1

`size`
:   How many records per page. Default is 20

`all_versions`
:   If to return all concept doi versions. If `true` returns all, if `false` returns the latest. Default is `ture`

`sort`
:   Which sorting to apply on the records. Default is `mostrecent`. Possible sortings "bestmatch", "mostrecent", "-mostrecent" (ascending), "version", "-version" (ascending).







Method `clean_doi()`
cleans the doi record for query
Usage
zenodoArchive$clean_doi(doi = self$doi)


Arguments


`doi`
:   The zenodo archive doi



Returns
the clean doi





Method `zenodo_query()`
Query the zenodo archive according to the initial parameters.
Usage
zenodoArchive$zenodo_query(...)


Arguments


`...`
:   Excepts the self created by `initialize`



Returns
a list with the query values.





Method `get_versions()`
Extract all concept doi available versions.
Usage
zenodoArchive$get_versions(...)


Arguments


`...`
:   Excepts the self created by `initialize`



Returns
a data.frame of the available versions.





Method `get_version_files()`
get the chosen doi archive version available files
Usage
zenodoArchive$get_version_files(version = "latest")


Arguments


`version`
:   which archive version files to get. Default to latest. To see all available version use `get_versions`



Returns
a list of the available files in the archive version.





Method `download_zenodo_files()`
get the chosen doi archive version available files
Usage
zenodoArchive$download_zenodo_files(
file = NULL,
path = tempdir(),
version = "latest",
all_files = F,
get_file_path = F,
quite = F
)


Arguments


`file`
:   If supplied, downloads the specific file from the archive.

`path`
:   The output folder for saving the archive files. Default is to a temporary directory.

`version`
:   which archive version files to get. Default to latest. To see all available version use `get_versions`

`all_files`
:   Logical (FALSE by default). Do you want to download all files in the archive.

`get_file_path`
:   Logical (FALSE by default). Do you want to return the path for the file downloaded.

`quite`
:   Logical (FALSE by default). Do you want to suppress informative messages



Returns
If `get_file_path` is TRUE, the function returns the path to the archive file





Method `clone()`
The objects of this class are cloneable with this method.
Usage
zenodoArchive$clone(deep = FALSE)


Arguments


`deep`
:   Whether to make a deep clone.







Examples
-------------------

```R
### Not run:
zenodo_archive <- zenodoArchive$new(
# doi = "10.5281/zenodo.7401189"
# )
# 
# # view available version ins the archive
# archive_versions <- zenodo_archive$get_versions()
# 
# # Getting the available files in the latest zenodo archive version
# files <- zenodo_archive$get_version_files()
# 
# # downloading the first file from the latest archive version
# zenodo_archive$download_zenodo_files()
```

**Error**: <text>:14:0: unexpected end of input
12: # # downloading the first file from the latest archive version
13: # zenodo_archive$download_zenodo_files()
   ^






