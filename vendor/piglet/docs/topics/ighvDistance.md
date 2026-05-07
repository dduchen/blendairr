**ighvDistance** - *Germline set alleles distance*

Description
--------------------

Calculates the distance between pairs of alleles based on their aligned germline sequences.
The function assume the germline set sequence are at an even length.
If not the function will pad the sequences with to the longest sequence length with Ns.


Usage
--------------------
```
ighvDistance(germline_set)
```

Arguments
-------------------

germline_set
:   A character list of the IMGT aligned IGHV allele sequences. See details for curating options.




Value
-------------------

A `matrix` of the computed distances between the alleles pairs.


Details
-------------------

The aligned IMGT IGHV allele germline set can be download from the IMGT site [https://www.imgt.org/](https://www.imgt.org/) under the section genedb.









