# Goal: import conventional LC-MS proteomics IDs for MALDI-MS1 HiTMaP matching

## Current state (2026-09-28)

Implemented and exported `import_lcms_proteomics_candidates()` in
`R/z_lcms_proteomics_import.R`.

It accepts a data frame or file and auto-detects these common ID exports:

| Producer | Preferred input | Identity fields used |
| --- | --- | --- |
| DIA-NN | Main `report.parquet` (or `report.tsv`) | `Stripped.Sequence`, `Modified.Sequence`, `Protein.Group`, `Protein.Ids` |
| FragPipe | `psm.tsv`/peptide table with proteins and `Assigned Modifications` | `Peptide` or `Peptide Sequence`, `Protein`, `Mapped Proteins` |
| MaxQuant | `evidence.txt` | `Sequence`, `Modifications`, `Leading razor protein`/`Proteins` |
| Proteome Discoverer | exported PSM or peptide-group CSV/TSV | `Annotated Sequence`, `Modifications`, `Master Protein Accessions` |
| TPP | `interact.pep.xml` or ProteinProphet `interact.prot.xml` | pepXML `search_hit`, or protXML protein-group membership, peptide sequence, and `modification_info` |

Sources checked when defining those layouts:

- DIA-NN's [main-output reference](https://github.com/vdemichev/DiaNN) documents the current main report and its sequence/protein columns.
- FragPipe's [result-file guide](https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_outputs.html) documents `combined_peptide.tsv`, `combined_modified_peptide.tsv`, and the peptide/protein fields.
- The PRIDE [submission-format guidance](https://github.com/PRIDE-Archive/pride-resources/blob/main/guidelines/pride-submission-formats-guidelines.md) identifies MaxQuant `evidence.txt` as the feature-level ID/quantification table.
- Thermo's [Proteome Discoverer 3.3 guide](https://docs.thermofisher.com/r/Proteome-Discoverer-3.3-User-Guide/2253872651v2en-US2254007691) documents peptide modifications, protein accessions, and peptide positions.
- The [pepXML schema documentation](https://tools.proteomecenter.org/formats/pepXML/pepXML_xmlspy_docs.html) documents `search_hit`, peptide sequence, alternative proteins, and modification positions.

## Import contract

- Candidate identity is **stripped peptide sequence + PTM composition**. PTM site positions are intentionally discarded, because current HiTMaP IMS is MALDI MS1 and cannot localise PTMs.
- For TPP, use `interact.pep.xml` when peptide-level evidence is the desired source. `interact.prot.xml` is also accepted when ProteinProphet protein groups are the desired source; its `ProteinProphet_group_<group_number>` label is retained in `Protein`, while all group member accessions are retained in `Protein_members` for FASTA mapping.
- Multiple rows/localisations with the same identity collapse into one candidate. `source_rows` records how many source rows were collapsed.
- Protein groups/members are retained as semicolon-separated `Protein` and `Protein_members`; they are not exploded into duplicate IMS formulas.
- `formula` and `neutral_mass` have **no adduct atoms**. `mz` is separately calculated as neutral mass + proton for the expected singly charged MALDI precursor. `adduct = "M+H"` is an ion-coordinate label, not part of `formula`.
- The importer calculates the base peptide formula directly, then resolves each supported PTM's elemental composition and monoisotopic mass from HiTMaP's packaged UniMod data. It does **not** invoke the legacy candidate-generation workflow.
- The importer rejects unsupported PTMs/mass shifts instead of silently creating an incorrect formula.

Supported formula-bearing PTMs: Carbamidomethyl, Oxidation, Acetyl, Phospho, Deamidated, Methyl, Dimethyl, Trimethyl, Formyl, Biotin, and N-terminal pyro-Glu from Q/E. It recognises their common names, relevant UniMod IDs, and common numeric mass annotations; formula composition and mass validation come from the packaged UniMod record rather than a second hard-coded composition table.

## Verification completed

- Parsed the new R source.
- Ran synthetic DIA-NN data with two protein groups for the same oxidised peptide: it collapsed to one `Peptide + Oxidation` row and retained the protein membership.
- Ran synthetic FragPipe numeric modifications, MaxQuant named modifications, and Proteome Discoverer annotated sequences: all produced the same Carbamidomethyl + Oxidation formula.
- Ran a synthetic TPP absolute-residue-mass annotation: it converted the TPP oxidation mass to the neutral formula.
- Added and validated both `interact.pep.xml` and `interact.prot.xml` modes. The protXML parser carries ProteinProphet group membership, primary peptides, and indistinguishable peptide variants. A TPP fixture was imported, mapped to its FASTA, and passed through `protein_scoring()` successfully.
- Downloaded a current public DIA-NN parquet report (100,396 rows) to `C:\Users\gguo484\Documents\expdata\lcms_import_validation`. The importer collapsed it to 23,843 candidates and correctly handled both bare unmodified `Modified.Sequence` values and UniMod 1/4/35 annotations.
- Mapped 51 real DIA-NN candidates from two proteins to the existing human FASTA in `expdata`; every mapped peptide had valid start/end coordinates.
- Exercised `protein_scoring()` using the FASTA-mapped real candidates and synthetic IMS match scores. It returned two protein-score rows. This exposed and fixed the invalid vector default for `scoretype`; `median` and `sum` now have implementations.

## Documentation and code review (2026-09-28)

The producers document the following authoritative fields:

- DIA-NN documents `Modified.Sequence`, `Stripped.Sequence`, `Protein.Group`, and `Decoy` in the main report. Its normal UniMod-form modified sequence is sufficient to identify a PTM, but an `--original-mods` export is not necessarily UniMod-normalised.
- FragPipe documents `combined_modified_peptide.tsv` / `combined_ion.tsv` as the modified-peptide tables. `Assigned Modifications` contains a mass, residue, and peptide location, and specifically represents variable modifications; `Modified Sequence`, `Protein`, and `Mapped Proteins` are also documented. Fixed modifications therefore require the run's `fragger.params` (or another explicit mapping) to be formula-complete.
- MaxQuant documents `Sequence`, `Modified sequence`, `Modifications`, `Proteins`, `Leading razor protein`, and `Reverse` in its output tables. The output also flags `Potential contaminant`.
- Proteome Discoverer documents the Peptide Groups export fields used here: `Annotated Sequence`, `Modifications`, `Master Protein Accessions`, `Protein Accessions`, and `Is Decoy`.
- pepXML documents `search_summary/aminoacid_modification` and `terminal_modification` with `massdiff`, residue/terminus, and description, as well as per-hit `modification_info`. protXML documents ProteinProphet protein groups and peptide modification information.

The current code is a useful compatibility importer, but it is **not yet a strict, provenance-backed importer**. Before calling it complete, address these review findings:

### Identification confidence policy (updated after user clarification)

The intended DIA-NN, FragPipe tabular, MaxQuant, and Proteome Discoverer exports have already passed their producing workflows' FDR filtering. Do not apply another default FDR filter to these tables; retain the upstream threshold and filtering provenance when available. DIA-NN's generic **main report** is not universally prefiltered at 1% (its default precursor FDR can be 5%), so users should check their export settings before calling imported candidates 1% FDR. `PEP` and `1 - PeptideProphet probability` are not q-values.

For TPP and FragPipe XML, `fdr_threshold = 0.01` now filters *before* candidate collapse, using the XML's own error-rate calibration. pepXML uses PeptideProphet/iProphet `roc_data_point` for `charge="all"` and the corresponding per-hit probability; protXML uses ProteinProphet `protein_summary_data_filter` and group probability. The least restrictive documented probability cutoff with estimated error <= the requested threshold is selected. If the XML lacks a usable curve, supply a run-calibrated `probability_cutoff` (such as FragPipe's `filter.log` cutoff) or explicitly set `fdr_threshold = NULL` to import without confidence filtering. Raw FragPipe `*.pepXML` search output lacks validation probabilities and cannot meet the XML default. PepXML and protXML filters control PSM and protein-group confidence, respectively; they do not establish 1% FDR for the final collapsed peptide/PTM identities.

Namespace-qualified pepXML/protXML synthetic fixtures in `C:\Users\gguo484\Documents\expdata\lcms_import_validation` verify the 1% and 5% cutoffs. Obtain real XML plus run logs for an end-to-end regression test.

## Real PRIDE fixtures checked (2026-09-28)

Downloaded original public results into `C:\Users\gguo484\Documents\expdata\lcms_import_validation` (not committed to the repository):

| Producer | PRIDE project and local file | Result |
| --- | --- | --- |
| MaxQuant | [PXD003133](https://www.ebi.ac.uk/pride/archive/projects/PXD003133), `PXD003133_txt_100min.zip`; extracted `txt_100min/evidence.txt`, `mqpar.xml`, `parameters.txt` | Full `evidence.txt`: 90,880 source rows, 90,723 retained after reverse-hit exclusion, 29,155 collapsed candidates. |
| FragPipe | [PXD062399](https://www.ebi.ac.uk/pride/archive/projects/PXD062399), `PXD062399_DDA.rar`; extracted `DDA/mc38_neg_1/psm.tsv`, `interact-mc38_neg1.pep.xml`, `filter.log`, plus `DDA/fragger.params` | Full `psm.tsv`: 3,740 rows, 2,664 collapsed candidates. A real C-containing row had Carbamidomethyl and calculated mass matching source within 0.0001 Da. XML default 1% path cannot infer a cutoff because this file lacks a usable ROC curve; `filter.log` reports a 0.728627 PSM probability cutoff at 0.97% FDR. Supplying that cutoff still exposes an unsupported absolute-mass shift, so XML is not yet fully importable. |
| Proteome Discoverer | [PXD064255](https://www.ebi.ac.uk/pride/archive/projects/PXD064255), `PXD064255_PeptideGroups.txt` | Full 3,943-row table imports to 3,943 candidates, including 150 Biotin-bearing rows after resolving canonical Biotin from UniMod ID 3. A Biotin+Oxidation candidate's calculated MH+ agrees with the source `Theo. MH+ [Da]` within 0.000004 Da. Bracketed `[R].PEPTIDE.[D]` notation and `2x` PTM counts required fixes. |
| TPP | [PXD004754](https://www.ebi.ac.uk/pride/archive/projects/PXD004754), `PXD004754_replicate5-interact.pep.xml` and `.prot.xml`; [PXD001080](https://www.ebi.ac.uk/pride/archive/projects/PXD001080), `PXD001080_interact-mock1.pep.xml` and `.prot.xml` | PXD004754 has stable-isotope K/R labels (~6 Da) outside the current formula engine. PXD001080's Comet configuration uses +0.998 Da on Q/N, likely intending deamidation, although canonical deamidation is +0.984016 Da. With an explicit UniMod:7 mapping and 0.02 Da tolerance, its full pepXML and protXML now import; see the follow-up results below. |

The real-file pass also found and fixed a package-data lookup error: `unimod.list.rda` is installed at `system.file("data", ..., package = "HiTMaP")`, but `utils::data()` could not load it. Real downloads are validation fixtures, not automatically safe to commit because of size and provenance.

Complete successful import outputs are saved beside the source files as `PXD003133_MaxQuant_candidates.csv`, `PXD062399_FragPipe_candidates.csv`, and `PXD064255_PD_candidates.csv`. Each was read back successfully; all candidate formulas and neutral masses were present and finite, and every `mz - neutral_mass` equalled the expected proton mass within 0.00001 Da. The real PD run also exposed a `2x` PTM-count parsing error, now corrected; a real MaxQuant two-oxidation row matches its source mass within 0.000004 Da.

1. Replace nearest-mass matching in `modification_counts()` with a source-provided, exact UniMod ID or an explicit user mapping. A mass shift alone is not a unique PTM identity.
2. Make UniMod resolution dynamic for any record with an elemental composition, rather than rejecting every ID outside the small initial supported set.
3. For FragPipe, require or accept `fragger.params` when fixed modifications are needed. Do not treat variable-only `Assigned Modifications` as the complete formula state.
4. For TPP pepXML, read the `search_summary` modification declarations and resolve `modification_info` against residue/terminus plus declared mass difference. Do not use nearest mass. For protXML, accept the accompanying pepXML/search configuration when modification definitions are not present.
5. Use the downloaded real TPP/FragPipe XML for regression testing of the XML confidence filter and naming variants (`interact.prot.xml` and `interact-prot.xml`). The current synthetic fixtures cover a default namespace. Address the observed unsupported isotope label and noncanonical +0.998 Da shifts by reading declared search modifications, not by nearest-mass guessing. Also correct C-terminal absolute-mass parsing.
6. For Proteome Discoverer, recognise the documented `Is Decoy` column. For MaxQuant, exclude the documented `Potential contaminant` rows by default or expose an explicit option. For TPP/FragPipe XML, recognise the run's decoy protein prefix before returning candidates; `drop_decoys` currently only handles an explicit table flag.
7. Turn the downloaded real fixtures and sidecars (`fragger.params`, `mqpar.xml`/modification configuration, PD workflow export, pepXML search summary) into reproducible regression tests, using small derived samples with attribution where appropriate.

## Next AI session: concrete work items

1. Real original exports from all five producers are now downloaded under `expdata` as listed above. Add repeatable regression tests and resolve the documented PTM/sidecar blockers; TPP still has only a synthetic end-to-end success case.
2. Use `prepare_lcms_candidates_for_hitmap(candidates, database)` before legacy `protein_scoring()`. It explicitly expands semicolon-delimited groups against the LC-MS FASTA and adds `start`, `end`, `pro_end`, numeric protein records, and descriptions. Do not pass the collapsed importer output directly to scoring.
3. Confirm the desired ion convention with a MALDI IMS dataset: current output stores a neutral formula plus [M+H]+ `mz`, whereas existing legacy candidate generation often stores an adducted formula. Update PMF isotope scoring only if it demonstrably expects the latter.
4. Extend the PTM composition map only after validating UniMod compositions and an example export (especially TMT/iTRAQ labels, glycosylation, cross-links, and user-defined modifications).
5. Consider retiring or clearly deprecating the older unexported `Protein_feature_list_table_import()` in `R/Canidate_processing_dev.R`; it has a separate, fragile hard-coded importer and requires Java. Do not remove it without checking existing callers.

## Follow-up implementation (2026-09-28)

- Replaced numeric nearest-mass matching with `modification_map` entries keyed by residue or terminus and signed mass shift. Declared XML search modifications are checked against the chosen UniMod record's monoisotopic mass. The PXD001080 Comet `variable_mod2 = 0.998000 QN` setting appears to intend the common Q→E/N→D deamidation (+0.984016 Da); the source mass is 0.013984 Da high. The importer requires an explicit UniMod:7 mapping and `modification_mass_tolerance = 0.02`, warns about the mismatch, and computes canonical formulas. This is an interpretation of the legacy search setting, not a claim that +0.998 Da is the deamidation elemental mass.
- Any non-isotopic UniMod ID with a parseable elemental composition can now be used directly. The importer rejects isotope-bearing compositions because the current HiTMaP formula and isotope-scoring functions do not represent labeled atoms correctly; PXD004754 therefore remains blocked by its SILAC labels. It is safer to reject these rows than output an incorrect natural-abundance formula.
- FragPipe tables now use `fragger.params` for fixed modifications, including rows where `Assigned Modifications` omits them. Explicit fixed annotations in `psm.tsv` are not counted twice. The PXD062399 C+M row `TICSHVQNMIK` reproduces its source neutral mass. The real PXD062399 pepXML imports 2,718 candidates when supplied the `filter.log` cutoff and exact C/M/N-terminal UniMod mappings.
- MaxQuant `evidence.txt` now requires `mqpar.xml` to include its fixed modifications; `Potential contaminant` is excluded by default. Proteome Discoverer `Is Decoy` is recognized. XML imports use the declared decoy prefix, and both `.prot.xml` and `-prot.xml` filenames are detected. protXML requires its paired pepXML for modification declarations. C-terminal XML modification masses are resolved through the declared terminal mass difference.
- Added `tests/testthat/test-lcms-import.R` with small source-derived PXD003133, PXD062399, and PXD064255 rows, namespaced XML cases, and paired real PXD001080 unmodified and doubly deamidated records. It covers formula/mass agreement, fixed-modification accounting, decoys, dynamic UniMod IDs, confidence filtering, and XML filename variants. `testthat::test_local('.', filter = 'lcms-import')` passes. The full PXD001080 pepXML produced 153 candidates (13 deamidated), and its paired protXML produced 361 candidates (77 deamidated) with the explicit mapping and tolerance.
- Checked the existing `MouseBrain_Trypsin_FT/Summary folder/candidatelist.csv` MALDI candidates against `SCORE_PMF()`: the legacy `formula` is the [M+H]+ atom composition, while `pepmz` is neutral. `prepare_lcms_candidates_for_hitmap()` now preserves imported `neutral_formula` and supplies an ionic `formula` for legacy isotope scoring; the imported candidate table itself remains neutral. A regression test checks that the prepared formula's ionic mass matches `mz`.

Still open: isotope-aware formula and isotope-pattern scoring for PXD004754 SILAC labels, and the decision about the old Java-based importer. The older `Protein_feature_list_table_import()` has a caller in `README.Rmd`, so it was left in place.
