**piglet** - *The Program for Ig clusters (PIgLET) package*

Description
--------------------

PIgLET is a suite of computational tools that improves genotype inference and downstream AIRR-seq data analysis.
The package as two main tools. The first is Allele Clusters, this tool is designed to reduce the ambiguity within the IGHV alleles. The ambiguity
is caused by duplicated or similar alleles which are shared among different genes.
The second tool is an allele based genotype, that determined the presence of an allele based on
a threshold derived from a naive population.






Allele Similarity Cluster
-------------------


This section provides the functions that support the main tool of creating the allele similarity cluster form
an IGHV germline set.


+ [inferAlleleClusters](inferAlleleClusters.md):      The main function of the section to create the allele clusters based on a germline set.
+ [ighvDistance](ighvDistance.md):             Calculate the distance between IGHV aligned germline sequences.
+ [ighvClust](ighvClust.md):                Hierarchical clustering of the distance matrix from `ighvDistance`.
+ [generateReferenceSet](generateReferenceSet.md):     Generate the allele clusters reference set.
+ [plotAlleleCluster](plotAlleleCluster.md):        Plots the Hierarchical clustering.
+ [artificialFRW1Germline](artificialFRW1Germline.md):   Artificially create an IGHV reference set with framework1 (FWR1) primers.



Allele based genotype
-------------------


This section provides the functions to infer the IGHV genotype using
the allele based method and the allele clusters thresholds


+ [inferGenotypeAllele](inferGenotypeAllele.md):      Infer the IGHV genotype using the allele based method.
+ [assignAlleleClusters](assignAlleleClusters.md):     Renames the v allele calls based on the new allele clusters.
+ [germlineASC](germlineASC.md):              Converts IGHV germline set to ASC germline set.
+ [recentAlleleClusters](recentAlleleClusters.md):     Download the most recent version of the allele clusters table archive from zenodo.
+ [extractASCTable](extractASCTable.md):          Extracts the allele cluster table from the zenodo archive file.
+ [zenodoArchive](zenodoArchive.md):            An R6 object to query the zenodo api.



References
-------------------


1. ##










