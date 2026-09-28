# Rows below are reduced from the original PRIDE exports PXD062399 (FragPipe),
# PXD003133 (MaxQuant), and PXD064255 (Proteome Discoverer). XML is a small
# search-summary fixture using PXD062399's declared masses and decoy prefix.
# https://www.ebi.ac.uk/pride/archive/projects/PXD062399
# https://www.ebi.ac.uk/pride/archive/projects/PXD003133
# https://www.ebi.ac.uk/pride/archive/projects/PXD064255

test_that("fixed and variable modifications produce the reported precursor mass", {
  params <- tempfile(fileext = ".params")
  writeLines(c("add_C_cysteine = 57.02146", "decoy_prefix = rev_"), params)
  frag <- data.frame(Peptide = "TICSHVQNMIK", Protein = "sp|P51410|RL9_MOUSE",
    `Assigned Modifications` = "3C(57.0215),9M(15.9949)", check.names = FALSE)
  mapped <- c("C:+57.0215" = "4", "M:+15.9949" = "35")
  result <- import_lcms_proteomics_candidates(frag, format = "fragpipe",
    fragger_params = params, modification_map = mapped)
  expect_equal(result$Modification, "Carbamidomethyl;Oxidation")
  expect_equal(result$neutral_mass, 1345.648, tolerance = 0.001)
  expect_error(import_lcms_proteomics_candidates(frag, format = "fragpipe",
    fragger_params = params), "exact UniMod ID")

  mqpar <- tempfile(fileext = ".xml")
  writeLines("<MaxQuantParams><fixedModifications><string>Carbamidomethyl (C)</string></fixedModifications></MaxQuantParams>", mqpar)
  mq <- data.frame(Sequence = c("AAAGEFADDPCSSVK", "AAAGEFADDPCSSVK"),
    `Modified sequence` = c("_AAAGEFADDPCSSVK_", "_AAAGEFADDPCSSVK_"),
    Modifications = c("Unmodified", "Unmodified"), Proteins = c("P1", "P1"),
    Reverse = c("", ""), `Potential contaminant` = c("", "+"), check.names = FALSE)
  result <- import_lcms_proteomics_candidates(mq, format = "maxquant", mqpar = mqpar)
  expect_equal(result$Modification, "Carbamidomethyl")
  expect_equal(result$source_rows, 1L)
})

test_that("Proteome Discoverer decoys and Biotin annotations are handled", {
  pd <- data.frame(`Annotated Sequence` = c("[R].SPSSNDKHVMSR.[L]", "[R].SPSSNDKHVMSR.[L]"),
    Modifications = c("1xBiotin [K7]; 1xOxidation [M10]", "1xBiotin [K7]"),
    `Master Protein Accessions` = c("Q9H8V3", "REV_Q9H8V3"),
    `Is Decoy` = c("false", "true"), check.names = FALSE)
  result <- import_lcms_proteomics_candidates(pd, format = "proteome_discoverer")
  expect_equal(nrow(result), 1L)
  expect_equal(result$Modification, "Biotin;Oxidation")
  expect_equal(result$mz, 1586.705, tolerance = 0.001)
})

test_that("any formula-bearing UniMod ID can be imported", {
  diann <- data.frame(Stripped.Sequence = "PEPTIDE",
    Modified.Sequence = "PEP(UniMod:53)TIDE", Protein.Group = "P1")
  result <- import_lcms_proteomics_candidates(diann)
  expect_equal(result$Modification, "53")
  expect_true(is.finite(result$neutral_mass))
  diann$Modified.Sequence <- "PEP(UniMod:188)TIDE"
  expect_error(import_lcms_proteomics_candidates(diann), "isotope-aware formula scoring")
  diann$Modified.Sequence <- "PEP(TMT)TIDE"
  expect_error(import_lcms_proteomics_candidates(diann), "Cannot calculate a formula")
})

test_that("pepXML declarations identify modifications before confidence filtering", {
  pep <- tempfile(fileext = ".pep.xml")
  writeLines(c(
    '<msms_pipeline_analysis xmlns="http://regis-web.systemsbiology.net/pepXML">',
    '<analysis_summary analysis="peptideprophet"><peptideprophet_summary><roc_error_data charge="all"><roc_data_point min_prob="0.8" error="0.01"/></roc_error_data></peptideprophet_summary></analysis_summary>',
    '<msms_run_summary><search_summary>',
    '<aminoacid_modification aminoacid="C" massdiff="57.02146" mass="160.03065" variable="N"/>',
    '<aminoacid_modification aminoacid="M" massdiff="15.9949" mass="147.0354" variable="Y"/>',
    '<terminal_modification terminus="C" massdiff="-0.984016" mass="16.018724" variable="Y"/>',
    '<parameter name="decoy_prefix" value="rev_"/></search_summary>',
    '<search_hit hit_rank="1" peptide="TICSHVQNMIK" protein="P51410"><modification_info mod_cterm_mass="16.018724"><mod_aminoacid_mass position="3" mass="160.03065"/><mod_aminoacid_mass position="9" mass="147.0354"/></modification_info><analysis_result analysis="peptideprophet"><peptideprophet_result probability="0.9"/></analysis_result></search_hit>',
    '<search_hit hit_rank="1" peptide="TICSHVQNMIK" protein="rev_P51410"><modification_info><mod_aminoacid_mass position="3" mass="160.03065"/></modification_info><analysis_result analysis="peptideprophet"><peptideprophet_result probability="0.9"/></analysis_result></search_hit>',
    '</msms_run_summary></msms_pipeline_analysis>'), pep)
  mapped <- c("C:+57.0215" = "4", "M:+15.9949" = "35", "c:-0.9840" = "2")
  result <- import_lcms_proteomics_candidates(pep, modification_map = mapped)
  expect_equal(nrow(result), 1L)
  expect_equal(result$Protein, "P51410")
  expect_equal(result$Modification, "2;Carbamidomethyl;Oxidation")
  expect_error(import_lcms_proteomics_candidates(pep), "no UniMod identity")

  prot <- paste0(tempfile(), ".interact-prot.xml")
  writeLines(c(
    '<protein_summary xmlns="http://regis-web.systemsbiology.net/protXML">',
    '<protein_summary_data_filter min_probability="0.8" false_positive_error_rate="0.01"/>',
    '<protein_group group_number="1" probability="0.9"><protein protein_name="P51410">',
    '<peptide peptide_sequence="TICSHVQNMIK"><modification_info mod_cterm_mass="16.018724"><mod_aminoacid_mass position="3" mass="160.03065"/><mod_aminoacid_mass position="9" mass="147.0354"/></modification_info></peptide>',
    '</protein></protein_group></protein_summary>'), prot)
  grouped <- import_lcms_proteomics_candidates(prot, pepxml = pep,
    modification_map = mapped)
  expect_equal(grouped$Protein, "ProteinProphet_group_1")
  expect_equal(grouped$Protein_members, "P51410")
  expect_equal(grouped$Modification, "2;Carbamidomethyl;Oxidation")
})

test_that("paired real TPP XML samples use an explicit deamidation correction", {
  pep <- testthat::test_path("fixtures", "PXD001080_minimal.interact.pep.xml")
  prot <- testthat::test_path("fixtures", "PXD001080_minimal.interact-prot.xml")
  mapped <- c("C:+57.0215" = "4", "M:+15.9949" = "35",
    "Q:+0.9980" = "7", "N:+0.9980" = "7")
  expect_error(import_lcms_proteomics_candidates(pep, modification_map = mapped),
    "does not match source modification")
  expect_warning(peptide <- import_lcms_proteomics_candidates(pep,
    modification_map = mapped, modification_mass_tolerance = 0.02),
    "candidate formulas use the UniMod composition")
  expect_warning(grouped <- import_lcms_proteomics_candidates(prot, pepxml = pep,
    modification_map = mapped, modification_mass_tolerance = 0.02),
    "candidate formulas use the UniMod composition")
  expect_setequal(peptide$Peptide, c("VLGVYEILPK", "ATSTTQYPAFNAIVYNEKK"))
  expect_equal(peptide$Modification[peptide$Peptide == "ATSTTQYPAFNAIVYNEKK"],
    "2xDeamidated")
  expect_equal(peptide$neutral_mass[peptide$Peptide == "ATSTTQYPAFNAIVYNEKK"],
    2147.070235 - 2 * (0.998 - 0.984016), tolerance = 0.001)
  expect_equal(sort(peptide$neutral_mass), sort(grouped$neutral_mass))
})

test_that("FASTA preparation supplies the ionic formula used by MALDI scoring", {
  fasta <- tempfile(fileext = ".fasta")
  writeLines(c(">P1", "MPEPTIDEK"), fasta)
  diann <- data.frame(Stripped.Sequence = "PEPTIDE",
    Modified.Sequence = "PEPTIDE", Protein.Group = "P1")
  candidate <- import_lcms_proteomics_candidates(diann)
  prepared <- prepare_lcms_candidates_for_hitmap(candidate, fasta)
  expect_equal(prepared$neutral_formula, candidate$formula)
  expect_equal(prepared$mz,
    OrgMassSpecR::MonoisotopicMass(get_atoms(prepared$formula), charge = 0) - 0.0005485799,
    tolerance = 0.00001)
  expect_equal(prepared$start, 2L)
})
