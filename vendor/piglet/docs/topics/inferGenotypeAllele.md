**inferGenotypeAllele** - *Allele similarity cluster based genotype inference*

Description
--------------------

`inferGenotypeAllele` infer an individual's genotype based on the allele-base method.
The method utilize the allele specific threshold to determine the presence of an allele in the genotype.
More specifically, the absolute frequency of each allele is calculated and checked against the threshold.


Usage
--------------------
```
inferGenotypeAllele(
data,
alleleClusterTable,
v_call = "v_call",
single_assignment = FALSE,
germline_db = NA,
find_unmutated = FALSE,
seq = "sequence_alignment"
)
```

Arguments
-------------------

data
:   data.frame in AIRR format, containing V allele calls from a single subject and the sample IMGT-gapped V(D)J sequences under seq.

alleleClusterTable
:   A data.frame of the allele similarity clusters thresholds.

v_call
:   name of the V allele call column. Default is `v_call`

single_assignment
:   if TRUE, the method only considers sequence with single assignment for the genotype inference.

germline_db
:   named vector of sequences containing the germline sequences named in V allele calls and the alleleClusterTable. Only required if find_unmutated is TRUE.

find_unmutated
:   if TRUE, use germline_db to find which samples are unmutated. Not needed if V allele calls only represent unmutated samples.

seq
:   name of the column in data with the aligned, IMGT-numbered, V(D)J nucleotide sequence. Default is sequence_alignment.




Value
-------------------

A a data.frame with the inferred V genotype. The table contains the following columns:<table><tr><td>
gene </td>
<td> alleles </td>
<td> imgt_alleles </td>
<td> counts </td>
<td> absolute_fraction </td>
<td> absolute_threshold </td>
<td> genotyped_alleles </td>
</tr><tr><td> genotype_imgt_alleles </td>
<td>
allele cluster </td>
<td> the present alleles </td>
<td> the imgt nomenclature </td>
<td> the number of reads </td>
<td> the absolute fraction </td>
<td> the population driven allele </td>
</tr><tr><td> the alleles which </td>
<td> the imgt nomenclature </td>
<td>
</td>
<td> in the repertoire </td>
<td> of the alleles </td>
<td> for each alleles </td>
<td> of the alleles </td>
</tr><tr><td> thresholds for genotype presence </td>
<td> entered the genotype </td>
<td> of the alleles </td>
<td>
gene </td>
<td> alleles </td>
<td> imgt_alleles </td>
<td> counts </td>
</tr></table>


Details
-------------------

In naive repertoires, allele calls where more than one assignment is assigned is rare. Hence, in case the data represents the naive repertoire of a subject
it is recommended to use the `find_unmutated=TRUE` option, to remove mutated sequences. For non-naive population, the allele calls in cases of multiple assignment
are treated as belonging to all groups.



Examples
-------------------

```R
### Not run:
# loading TIgGER AIRR-seq b cell data
# data <- tigger::AIRRDb
# 
# # getting the archive
# asc_archive <- recentAlleleClusters(doi="10.5281/zenodo.7401239", get_file = TRUE)
# 
# # extracting the allele cluster table
# allele_cluster_table <- extractASCTable(archive_file = asc_archive)
# 
# data(HVGERM)
# 
# # reforming the germline set
# asc_germline <- germlineASC(allele_cluster_table, germline = HVGERM)
# 
# # assigning the ASC alleles
# asc_data <- assignAlleleClusters(data, allele_cluster_table)
# 
# # inferring the genotype
# asc_genotype <- inferGenotypeAllele(asc_data,
# alleleClusterTable = allele_cluster_table,
# germline_db = asc_germline, find_unmutated=T)
```



See also
-------------------

[inferAlleleClusters](inferAlleleClusters.md) will infer the allele clusters based on a supplied V reference set and set the default allele threshold of 1e-04.
See [recentAlleleClusters](recentAlleleClusters.md) to obtain the latest version of the IGHV allele clusters and the naive population based allele threshold.






