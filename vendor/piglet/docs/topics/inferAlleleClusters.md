**inferAlleleClusters** - *Allele similarity cluster*

Description
--------------------

A wrapper function to infer the allele clusters. See details for cluster inference


Usage
--------------------
```
inferAlleleClusters(
germline_set,
trim_3prime_side = 318,
mask_5prime_side = 0,
family_threshold = 75,
allele_cluster_threshold = 95
)
```

Arguments
-------------------

germline_set
:   Either a character vector of strings representing Ig sequence alleles, or a path to to the germline set file (must be gapped by IMGT scheme for optimal results).

trim_3prime_side
:   To which nucleotide position to trim the sequences. Default is 318; NULL will take the entire sequence length.

mask_5prime_side
:   Mimic short sequence libraries, gets the length of nucleotides to mask from the 5' side, the staring position. Default is 0.

family_threshold
:   The similarity threshold for the family level. Default is 75.

allele_cluster_threshold
:   The similarity threshold for the allele cluster level. Default is 95.




Value
-------------------

An object of type [GermlineCluster](GermlineCluster-class.md) that includes the following slots:


Details
-------------------

The distance between pairs of the alleles germline set sequences is calculated, then the alleles are clustered based on two similarity thresholds.
One for the family cluster and the other for the allele cluster. Then the new allele cluster names are generated and the germline set sequences are renamed and duplicated alleles are removed.

The allele cluster names are by the following scheme:
IGHVF1-G1*01 - IGH = chain, V = region, F1 = family cluster numbering,
G1 - allele cluster numbering, and 01 = allele numbering (given by clustering order, no connection to the expression)

To plot the allele clusters dendrogram use the `plot` function on the [GermlineCluster](GermlineCluster-class.md) object


Slots
-------------------



`germlineSet`
:   
+  A character vector with the modified germline set (3' trimming and 5' masking).


`alleleClusterSet`
:   
+  A character vector of renamed input germline set to the ASC name scheme (Without 3' and 5' modifications).


`alleleClusterTable`
:   
+  A data.frame of the allele similarity cluster with the new names and the default thresholds.


`threshold`
:   
+  A list of the input family and allele cluster similarity thresholds.


`hclustAlleleCluster`
:   
+  An hclust object of the germline set hierarchical clustering,




Examples
-------------------

```R
### Not run:
# load the initial germline set
# 
# data(HVGERM)
# 
# germline <- HVGERM
# 
# asc <- inferAlleleClusters(germline)
# 
# ## plotting the clusters
# 
# plot(asc)
```



See also
-------------------

By using the plot function on the returned object, a colorful visualization of the allele clusters dendrogram and threshold is received






