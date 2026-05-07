**artificialFRW1Germline** - *FWR1 artificial dataset generator*

Description
--------------------

A function to artificially create an IGHV reference set with framework1 (FWR1) primers (see Details).


Usage
--------------------
```
artificialFRW1Germline(
germline_set,
mask_primer = T,
trimm_primer = F,
quite = F
)
```

Arguments
-------------------

germline_set
:   A germline set distance matrix created by `ighvDistance`.

mask_primer
:   Logical (TRUE by default). If to mask with Ns the region of the primer from the germline sequence

trimm_primer
:   Logical (FALSE by default). If to trim the region of the primer from the germline sequence. If TRUE then, mask_primer is ignored.

quite
:   Logical (FALSE by default). Do you want to suppress informative messages




Value
-------------------

A `list` with the input germline set allele and the trimmed/masked sequences.


Details
-------------------

The FRW1 primers used in this function were taken from the BIOMED-2 protocol. For more information on the protocol and primer design go to:
van Dongen, J., Langerak, A., Brüggemann, M. et al. Design and standardization of PCR primers and protocols for detection of clonal immunoglobulin and
T-cell receptor gene recombinations in suspect lymphoproliferations: Report of the BIOMED-2 Concerted Action BMH4-CT98-3936.
Leukemia 17, 2257–2317 (2003). https://doi.org/10.1038/sj.leu.2403202Van Dongen, J. J. M., et al. "Design and standardization of PCR primers and protocols for detection of clonal immunoglobulin and T-cell
receptor gene recombinations in suspect lymphoproliferations: report of the BIOMED-2 Concerted Action BMH4-CT98-3936."
Leukemia 17.12 (2003): 2257-2317.









