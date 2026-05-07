**GermlineCluster-class** - *Output of inferAlleleClusters*

Description
--------------------

`GermlineCluster` contains output from [inferAlleleClusters](inferAlleleClusters.md) function.
It includes the allele cluster table, the germline set hierarchical clustering,
and the threshold parameters.


Usage
--------------------
```
"plot"(x, y = NULL)
```

Arguments
-------------------

x
:   GermlineCluster object

y
:   not in use. missing.




Methods (by generic)
-------------------


+  `plot(x = GermlineCluster, y = missing)`: Plot the dendrogram for the allele clusters.



Slots
-------------------



`germlineSet`
:   the original germline set provided

`alleleClusterSet`
:   the renamed germline set with the allele clusters

`alleleClusterTable`
:   the allele cluster table

`hclustAlleleCluster`
:   the hierarchical clustering object of the germline set.

`threshold`
:   the threshold used for the family and the allele clusters.




See also
-------------------

 \link{inferAlleleClusters}







