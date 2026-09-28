The two `PXD001080_minimal` XML files are reduced from the original
`interact-mock1.pep.xml` and `interact-mock1.prot.xml` exports in
[PRIDE PXD001080](https://www.ebi.ac.uk/pride/archive/projects/PXD001080).
They retain the original search modification declarations, decoy prefix,
error-curve points, ProteinProphet filter point, and paired records for an
unmodified peptide (`VLGVYEILPK`) and a doubly deamidated peptide
(`ATSTTQYPAFNAIVYNEKK`). Unrelated records and XML metadata were removed to
keep the regression fixture small. The source Comet setting uses +0.998 Da
for Q/N, whereas UniMod deamidation is +0.984016 Da. The test explicitly maps
that setting to UniMod:7 with a 0.02 Da tolerance and checks the corrected
formula mass; the importer warns about the discrepancy.
