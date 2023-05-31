To study the phylogenetic relationships of LRR-RK sequences, we first extracted their PK domains and aligned them with the PK domains of reference LRR-RKs using *MAFFT*. The resulting alignment was then used to construct phylogenetic trees using *FastTree*. This was done for all four sets of sequences deposited in the supplementary material: LRR-RKs, Subgroup XII, PERU, and FLS2.

To extract the PERU clade, we performed multiple iterations of alignment and tree construction. In each iteration, we identified a well-supported major branch and extracted it using *Dendroscope*. We then realigned the extracted sequences and constructed a new tree. This process was repeated until we obtained a tree with good resolution for comparing PERU with its closely related sequences (average pairwise alignment ~77%; PK domain only).

Finally, we extracted the full-length sequences for the final clade to obtain a more concise alignment and a better overview of the phylogenetic relationships of both the PERU and FLS2 clades.

```bash

# alignment
mafft --anysymbol sequence.fasta > alignment.afa

# tree construction
FastTree alignment.afa > tree.newick

```

The full-length sequences of PERU and FLS2 were checked for the presence of a signal peptide using SignalP 6 and the LRRNT domain. FLS2 sequences were also filtered to remove sequences shorter than 1150 amino acids.