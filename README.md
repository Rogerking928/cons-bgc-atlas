# Biosynthetic gene clusters of twelve coagulase-negative staphylococci

Analysis code for a complete-genome atlas of the non-ribosomal peptide
synthetase and RiPP biosynthetic gene clusters of the coagulase-negative
staphylococci: 554 complete RefSeq genomes of twelve species that colonise and
infect people, retrieved 20 September 2026. The repository also holds the
in-silico PCR typing scheme, the ESM-2 precursor-peptide detector, and the
checks reported in the paper, including the comparison with antiSMASH and
AMRFinderPlus.

## Running

```
pip install -r requirements.txt
cd code
python3 s1_fetch_reports.py      # assembly reports for the twelve species
python3 s2_download.py           # protein.faa + genomic.gff
python3 s2b_download_fna.py      # genomic.fna (antiSMASH reads nucleotides)
python3 s3_fetch_hmms.py         # the 45 Pfam profiles
python3 s4_hmmscan.py            # hmmsearch --cut_ga
python3 s5_parse.py              # -> out/cds.tsv.gz, out/domains.tsv.gz
python3 s6_call_clusters.py      # cluster calling
python3 s9_sketch_tree.py        # whole-proteome MinHash, lineage strata
python3 s10_metadata.py          # BioSample metadata
python3 s11_families.py          # families   s12_matrix.py  genome x family
python3 s53_amrfinder.py         # AMRFinderPlus over the panel
python3 s53b_amr_from_amrfinder.py
python3 s15_colocalise.py        # cluster-to-determinant distances
python3 s13_stats.py             # associations, stratified estimates
python3 s43_distance_null.py     # the permutation null
python3 s14_primers.py           # primer design + in-silico PCR
python3 s17_precursor_dataset.py # ... s52: the ESM-2 detector and its checks
python3 s54_antismash.py         # antiSMASH over the panel
python3 s55_antismash_compare.py # atlas against antiSMASH, scope-aligned
python3 s56_mibig_kcb.py         # knownclusterblast against MIBiG
python3 s58_mibig_join.py
python3 s57_pseudogenes.py       # pseudogenised anchor genes
python3 tables.py                # the tables
python3 f1_atlas.py ... f7_antismash.py   # the figures
python3 make_doc.py              # Tables_and_Figures.docx
python3 make_xlsx.py             # Source_Data.xlsx
python3 make_supp_data.py        # Supplementary_Data.xlsx
```

Scripts are numbered in the order they run; each writes to `out/` and later
steps read those files. Permutation and bootstrap steps use fixed seeds.
Comments in the code are partly in Chinese.

## Data

No genome is redistributed here. Every assembly is a public NCBI RefSeq record
and the fetch scripts download them from the NCBI Datasets API; the panel is
defined by the inclusion rule in `s1_fetch_reports.py` (complete assemblies of
the twelve species, excluding those flagged atypical). The Pfam profiles come
from InterPro, the resistance reference database from AMRFinderPlus, the
reference peptides from UniProt, and the sequence-type schemes from PubMLST and
BIGSdb-Pasteur. antiSMASH and AMRFinderPlus are run as external programs and
are not vendored; the versions used are antiSMASH 8.0.4 and AMRFinderPlus 4.2.7
with database 2026-08-07.1.

## Software

Python 3.12; package versions in `requirements.txt`.

The analysis was written and run by the author in Python 3, with Claude
(Anthropic) and ChatGPT (OpenAI) used as assistants under the author's
direction; the author reviewed and tested all code and takes full
responsibility for the results it produces.

## Licence

MIT.
