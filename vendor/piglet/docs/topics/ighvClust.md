**ighvClust** - *Allele similarity clustering*

Description
--------------------

Cluster the distance matrix from `ighvDistance` to create the allele clusters based on two thresholds:
75% similarity which represents the family clustering and 95% similarity between alleles which represents the allele clusters


Usage
--------------------
```
ighvClust(
germline_distance,
family_threshold = 75,
allele_cluster_threshold = 95
)
```

Arguments
-------------------

germline_distance
:   A germline set distance matrix created by `ighvDistance`.

family_threshold
:   The similarity threshold for the family level. Default is 75.

allele_cluster_threshold
:   The similarity threshold for the allele cluster level. Default is 95.




Value
-------------------

A names list that includes the `data.frame` of the alleles clusters, the thresholds parameters and the
hierarchical clustering of the germline set.









