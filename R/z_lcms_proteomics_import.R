#' Import LC-MS proteomics identifications as MALDI-IMS peptide candidates
#'
#' Converts peptide-level identification output from DIA-NN, FragPipe, MaxQuant,
#' Proteome Discoverer, or TPP into neutral-formula HiTMaP candidates. PTM sites
#' and LC-MS charge states are deliberately discarded because MALDI MS1 cannot
#' resolve them. Rows collapse by stripped sequence and PTM composition while
#' retaining all reported protein/group identifiers. PTM elemental compositions
#' and mass checks are resolved from HiTMaP's packaged UniMod data; this does
#' not run the legacy candidate-generation workflow.
#'
#' @param input A data frame, delimited export, TPP/FragPipe pepXML, or
#'   ProteinProphet protXML file.
#' @param format Producer: `"auto"`, `"diann"`, `"fragpipe"`, `"maxquant"`,
#'   `"proteome_discoverer"`, or `"tpp"`.
#' @param drop_decoys Exclude marked decoys/reverse hits and XML accessions
#'   matching the run's decoy prefix.
#' @param output_file Optional CSV path for the candidates.
#' @param fdr_threshold Maximum model-estimated FDR for TPP/FragPipe XML input
#'   (default 0.01). Set to `NULL` to skip XML confidence filtering. Tabular
#'   exports are assumed to have been filtered by their producer.
#' @param probability_cutoff Optional probability cutoff calibrated externally
#'   for the requested XML FDR, used when its error curve is absent.
#' @param modification_map Named character vector mapping source modification
#'   annotations to exact UniMod IDs (for example `c("M:+15.9949" = "35")`).
#'   Numeric shifts without an identity are rejected.
#' @param modification_mass_tolerance Maximum permitted difference in Da
#'   between an explicitly mapped source shift and the UniMod mass. Increase
#'   only for a known imprecise search setting; the output uses UniMod chemistry.
#' @param fragger_params Optional MSFragger `fragger.params` path. Required for
#'   FragPipe tables when the search used fixed modifications.
#' @param pepxml Optional accompanying pepXML path supplying search modification
#'   declarations for protXML.
#' @param mqpar MaxQuant `mqpar.xml` path containing fixed modifications.
#' @param decoy_prefix Optional protein accession prefix used for decoys.
#' @param drop_contaminants Exclude MaxQuant potential contaminants.
#' @return A data frame with a neutral `formula` and `neutral_mass`. `mz` is the
#'   calculated [M+H]+ matching coordinate; formula contains no adduct atoms.
#' @export
import_lcms_proteomics_candidates <- function(
    input, format = c("auto", "diann", "fragpipe", "maxquant",
                      "proteome_discoverer", "tpp"),
    drop_decoys = TRUE, output_file = NULL, fdr_threshold = 0.01,
    probability_cutoff = NULL, modification_map = NULL,
    modification_mass_tolerance = 0.001, fragger_params = NULL,
    pepxml = NULL, mqpar = NULL, decoy_prefix = NULL, drop_contaminants = TRUE) {
  format <- match.arg(format)
  if (!is.null(modification_map) && (!is.character(modification_map) ||
      is.null(names(modification_map)) || any(!nzchar(names(modification_map))) ||
      anyDuplicated(names(modification_map)) ||
      any(!grepl("^[0-9]+$", modification_map)))) {
    stop("`modification_map` must be a named character vector of UniMod IDs.", call. = FALSE)
  }
  if (!is.numeric(modification_mass_tolerance) || length(modification_mass_tolerance) != 1L ||
      !is.finite(modification_mass_tolerance) || modification_mass_tolerance < 0)
    stop("`modification_mass_tolerance` must be a non-negative finite number.", call. = FALSE)
  mapped_id <- function(key) if (key %in% names(modification_map)) modification_map[[key]] else NULL
  if (!is.null(fdr_threshold) && (!is.numeric(fdr_threshold) || length(fdr_threshold) != 1L ||
      is.na(fdr_threshold) || !is.finite(fdr_threshold) || fdr_threshold < 0 || fdr_threshold > 1)) {
    stop("`fdr_threshold` must be one number between 0 and 1, or NULL.", call. = FALSE)
  }
  if (!is.null(probability_cutoff) && (!is.numeric(probability_cutoff) ||
      length(probability_cutoff) != 1L || is.na(probability_cutoff) ||
      !is.finite(probability_cutoff) || probability_cutoff < 0 || probability_cutoff > 1)) {
    stop("`probability_cutoff` must be one number between 0 and 1, or NULL.", call. = FALSE)
  }
  if (is.null(fdr_threshold) && !is.null(probability_cutoff)) {
    stop("`probability_cutoff` requires a non-NULL `fdr_threshold`.", call. = FALSE)
  }
  first_column <- function(x, candidates, required = TRUE) {
    found <- candidates[candidates %in% names(x)]
    if (length(found)) return(found[[1]])
    if (required) stop("Missing required column. Expected one of: ",
                       paste(candidates, collapse = ", "), call. = FALSE)
    NULL
  }
  read_table <- function(path) {
    if (tolower(tools::file_ext(path)) == "parquet") {
      if (!requireNamespace("arrow", quietly = TRUE)) {
        stop("DIA-NN parquet input requires the optional arrow package. ",
             "Install arrow or export report.tsv.", call. = FALSE)
      }
      return(as.data.frame(arrow::read_parquet(path)))
    }
    if (tolower(tools::file_ext(path)) %in% c("tsv", "txt")) {
      utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    } else utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  unimod_environment <- new.env(parent = emptyenv())
  unimod_path <- system.file("data", "unimod.list.rda", package = "HiTMaP")
  if (nzchar(unimod_path)) withCallingHandlers(
    load(unimod_path, envir = unimod_environment),
    warning = function(w) {
      if (grepl("strings not representable in native encoding", conditionMessage(w), fixed = TRUE))
        invokeRestart("muffleWarning")
    })
  if (!exists("unimod.list", envir = unimod_environment, inherits = FALSE)) {
    stop("The packaged UniMod data are unavailable; reinstall HiTMaP.", call. = FALSE)
  }
  unimod_nodes <- unname(get("unimod.list", envir = unimod_environment)$modifications)
  unimod_attributes <- data.frame(do.call(rbind, lapply(unimod_nodes, function(node) {
    attributes <- node$.attrs
    attributes[match(names(unimod_nodes[[1]]$.attrs), names(attributes))]
  })), stringsAsFactors = FALSE)
  unimod_ids <- c(Carbamidomethyl = "4", Oxidation = "35", Acetyl = "1",
    Phospho = "21", Deamidated = "7", Methyl = "34", Dimethyl = "36",
    Trimethyl = "37", Formyl = "122", Biotin = "3", `Pyro-glu from Q` = "28",
    `Pyro-glu from E` = "27")
  unimod_for <- function(name) {
    id <- if (grepl("^[0-9]+$", name)) name else unimod_ids[[name]]
    entry <- unimod_attributes[unimod_attributes$record_id == id, , drop = FALSE]
    if (nrow(entry) != 1L || is.na(entry$composition) || !nzchar(entry$composition)) {
      stop("UniMod has no unambiguous elemental composition for ", name, ".", call. = FALSE)
    }
    entry
  }
  warned_mass_mappings <- new.env(parent = emptyenv())
  validate_mapped_mass <- function(id, value, key) {
    difference <- abs(as.numeric(unimod_for(id)$mono_mass) - value)
    if (difference > modification_mass_tolerance)
      stop("UniMod:", id, " does not match source modification ", key,
        " (difference ", signif(difference, 6), " Da).", call. = FALSE)
    if (difference > 0.001 && !exists(key, envir = warned_mass_mappings, inherits = FALSE)) {
      warning("Source modification ", key, " differs from UniMod:", id, " by ",
        signif(difference, 6), " Da; candidate formulas use the UniMod composition.", call. = FALSE)
      assign(key, TRUE, envir = warned_mass_mappings)
    }
  }
  tpp_table <- function(path) {
    if (!requireNamespace("XML", quietly = TRUE)) stop("TPP XML import requires XML.", call. = FALSE)
    document <- XML::xmlParse(path)
    root_node <- XML::xmlRoot(document)
    root <- XML::xmlName(root_node)
    if (root == "protein_summary" && is.null(pepxml))
      stop("protXML requires its accompanying `pepxml` to account for search modifications.", call. = FALSE)
    namespace_uri <- as.character(XML::xmlNamespace(root_node))
    namespaces <- if (length(namespace_uri) && nzchar(namespace_uri)) c(hitmap = namespace_uri) else character()
    tag <- function(name) if (length(namespaces)) paste0("hitmap:", name) else name
    children <- function(node, name) XML::getNodeSet(node, paste0("./", tag(name)), namespaces = namespaces)
    descendants <- function(node, name) XML::getNodeSet(node, paste0(".//", tag(name)), namespaces = namespaces)
    declaration_document <- if (!is.null(pepxml) && root == "protein_summary") XML::xmlParse(pepxml) else document
    declarations <- XML::getNodeSet(declaration_document,
      "//*[local-name()='aminoacid_modification' or local-name()='terminal_modification']")
    declaration_key <- function(node) {
      site <- XML::xmlGetAttr(node, "aminoacid", default = "")
      if (!nzchar(site)) site <- tolower(XML::xmlGetAttr(node, "terminus", default = ""))
      sprintf("%s:%+.4f", site, as.numeric(XML::xmlGetAttr(node, "massdiff")))
    }
    declared_id <- function(node) {
      key <- declaration_key(node)
      description <- XML::xmlGetAttr(node, "description", default = "")
      id <- if (nzchar(description) && grepl("unimod:[0-9]+", description, ignore.case = TRUE))
        sub(".*unimod:([0-9]+).*", "\\1", description, ignore.case = TRUE) else mapped_id(key)
      if (is.null(id) && nzchar(description)) id <- mapped_id(description)
      if (is.null(id)) stop("Search modification ", key,
        " has no UniMod identity. Supply `modification_map`.", call. = FALSE)
      validate_mapped_mass(id, as.numeric(XML::xmlGetAttr(node, "massdiff")), key)
      id
    }
    xml_decoy_prefix <- decoy_prefix
    if (is.null(xml_decoy_prefix)) {
      parameters <- XML::getNodeSet(declaration_document, "//*[local-name()='parameter' and @name='decoy_prefix']")
      if (length(parameters)) xml_decoy_prefix <- XML::xmlGetAttr(parameters[[1]], "value")
    }
    xml_probability_cutoff <- function(points, probability_field, error_field, label) {
      if (!is.null(probability_cutoff)) return(probability_cutoff)
      if (!length(points)) {
        stop(label, " has no usable error curve. Supply `probability_cutoff` from the run's ",
             100 * fdr_threshold, "% FDR report, or set `fdr_threshold = NULL` ",
             "to import without FDR filtering.", call. = FALSE)
      }
      probabilities <- suppressWarnings(as.numeric(vapply(points, XML::xmlGetAttr, "", probability_field)))
      errors <- suppressWarnings(as.numeric(vapply(points, XML::xmlGetAttr, "", error_field)))
      valid <- is.finite(probabilities) & probabilities >= 0 & probabilities <= 1 &
        is.finite(errors) & errors >= 0 & errors <= fdr_threshold
      if (!any(valid)) stop("No ", label, " cutoff achieves FDR <= ", fdr_threshold, ".", call. = FALSE)
      min(probabilities[valid])
    }
    modification_detail <- function(node, sequence) {
      info <- children(node, "modification_info")
      mods <- character()
      sites <- strsplit(sequence, "", fixed = TRUE)[[1]]
      observed <- if (length(info)) children(info[[1]], "mod_aminoacid_mass") else list()
      if (length(info) && !length(observed) &&
          grepl("[A-Z]\\[[0-9.]+\\]", XML::xmlGetAttr(info[[1]],
            "modified_peptide", default = "")))
        stop("XML modified_peptide has residue masses but no mod_aminoacid_mass records.", call. = FALSE)
      positions <- if (length(observed)) as.integer(vapply(observed, XML::xmlGetAttr, "", "position")) else integer()
      for (decl in declarations) {
        site <- XML::xmlGetAttr(decl, "aminoacid", default = "")
        term <- tolower(XML::xmlGetAttr(decl, "terminus", default = ""))
        expected <- as.numeric(XML::xmlGetAttr(decl, "mass"))
        fixed <- identical(XML::xmlGetAttr(decl, "variable", default = "Y"), "N")
        if (nzchar(site)) {
          matching <- which(sites == site)
          for (position in matching) {
            found <- which(positions == position & abs(as.numeric(vapply(observed,
              XML::xmlGetAttr, "", "mass")) - expected) < 0.02)
            if (fixed || length(found)) mods <- c(mods, paste0("UniMod:", declared_id(decl)))
          }
        } else if (term %in% c("n", "c")) {
          field <- if (term == "n") "mod_nterm_mass" else "mod_cterm_mass"
          value <- if (length(info)) as.numeric(XML::xmlGetAttr(info[[1]], field, default = NA_character_)) else NA_real_
          if (fixed || (is.finite(value) && abs(value - expected) < 0.02))
            mods <- c(mods, paste0("UniMod:", declared_id(decl)))
        }
      }
      if (length(observed)) {
        matched <- vapply(observed, function(mass_node) {
          position <- as.integer(XML::xmlGetAttr(mass_node, "position"))
          mass <- as.numeric(XML::xmlGetAttr(mass_node, "mass"))
          any(vapply(declarations, function(decl) {
            identical(XML::xmlGetAttr(decl, "aminoacid", default = ""), sites[[position]]) &&
              abs(as.numeric(XML::xmlGetAttr(decl, "mass")) - mass) < 0.02
          }, logical(1)))
        }, logical(1))
        if (any(!matched)) stop("pepXML modification_info has a mass absent from search_summary.", call. = FALSE)
      }
      if (length(info)) for (field in c("mod_nterm_mass", "mod_cterm_mass")) {
        value <- as.numeric(XML::xmlGetAttr(info[[1]], field, default = NA_character_))
        if (is.finite(value) && !any(vapply(declarations, function(decl)
          abs(as.numeric(XML::xmlGetAttr(decl, "mass")) - value) < 0.02, logical(1))))
          stop("pepXML terminal modification is absent from search_summary.", call. = FALSE)
      }
      paste(mods, collapse = ";")
    }
    if (root == "msms_pipeline_analysis") {
      hits <- descendants(root_node, "search_hit")
      if (!length(hits)) stop("No search_hit records found in pepXML.", call. = FALSE)
      ranks <- vapply(hits, XML::xmlGetAttr, "", "hit_rank")
      hits <- hits[!nzchar(ranks) | ranks == "1"]
      if (!is.null(fdr_threshold)) {
        has_interprophet <- any(vapply(hits, function(hit) {
          results <- children(hit, "analysis_result")
          any(vapply(results, XML::xmlGetAttr, "", "analysis") == "interprophet")
        }, logical(1)))
        analysis <- if (has_interprophet) "interprophet" else "peptideprophet"
        summaries <- descendants(root_node, "analysis_summary")
        summaries <- summaries[vapply(summaries, XML::xmlGetAttr, "", "analysis") == analysis]
        roc <- unlist(lapply(summaries, function(summary) {
          curves <- descendants(summary, "roc_error_data")
          curves <- curves[vapply(curves, XML::xmlGetAttr, "", "charge") == "all"]
          unlist(lapply(curves, children, "roc_data_point"), recursive = FALSE)
        }), recursive = FALSE)
        cutoff <- xml_probability_cutoff(roc, "min_prob", "error", paste0(analysis, " pepXML"))
        probabilities <- vapply(hits, function(hit) {
          results <- children(hit, "analysis_result")
          results <- results[vapply(results, XML::xmlGetAttr, "", "analysis") == analysis]
          if (!length(results)) return(NA_real_)
          scored <- children(results[[1]], paste0(analysis, "_result"))
          if (!length(scored)) return(NA_real_)
          suppressWarnings(as.numeric(XML::xmlGetAttr(scored[[1]], "probability", default = NA_character_)))
        }, numeric(1))
        if (anyNA(probabilities)) stop("pepXML search hits lack ", analysis,
          " probabilities needed for FDR filtering.", call. = FALSE)
        hits <- hits[probabilities >= cutoff]
      }
      if (!length(hits)) return(data.frame(Protein = character(), Peptide = character(),
        Modification_detail = character(), stringsAsFactors = FALSE))
      return(do.call(rbind, lapply(hits, function(hit) {
      alternatives <- children(hit, "alternative_protein")
      proteins <- c(XML::xmlGetAttr(hit, "protein"),
                    vapply(alternatives, XML::xmlGetAttr, "", "protein"))
      is_decoy <- !is.null(xml_decoy_prefix) && nzchar(xml_decoy_prefix) &&
        all(startsWith(proteins, xml_decoy_prefix))
      if (drop_decoys && !is.null(xml_decoy_prefix) && nzchar(xml_decoy_prefix))
        proteins <- proteins[!startsWith(proteins, xml_decoy_prefix)]
      data.frame(Protein = paste(unique(proteins[nzchar(proteins)]), collapse = ";"),
        Peptide = XML::xmlGetAttr(hit, "peptide"),
        Modification_detail = modification_detail(hit, XML::xmlGetAttr(hit, "peptide")),
        Decoy = if (is_decoy) "true" else "",
        stringsAsFactors = FALSE)
      })))
    }
    if (root != "protein_summary") {
      stop("TPP XML must have a pepXML msms_pipeline_analysis or ProteinProphet protein_summary root.", call. = FALSE)
    }
    groups <- descendants(root_node, "protein_group")
    if (!length(groups)) stop("No protein_group records found in protXML.", call. = FALSE)
    if (!is.null(fdr_threshold)) {
      points <- descendants(root_node, "protein_summary_data_filter")
      cutoff <- xml_probability_cutoff(points, "min_probability", "false_positive_error_rate",
        "ProteinProphet protXML")
      probabilities <- suppressWarnings(as.numeric(vapply(groups, XML::xmlGetAttr, "", "probability")))
      if (anyNA(probabilities)) stop("protXML protein groups lack probabilities needed for FDR filtering.", call. = FALSE)
      groups <- groups[probabilities >= cutoff]
    }
    if (!length(groups)) return(data.frame(Protein = character(), Protein_members = character(),
      Peptide = character(), Modification_detail = character(), stringsAsFactors = FALSE))
    rows <- list()
    for (group in groups) {
      proteins <- children(group, "protein")
      members <- unique(unlist(lapply(proteins, function(protein) {
        c(XML::xmlGetAttr(protein, "protein_name"),
          vapply(children(protein, "indistinguishable_protein"), XML::xmlGetAttr, "", "protein_name"))
      })))
      is_decoy <- !is.null(xml_decoy_prefix) && nzchar(xml_decoy_prefix) &&
        all(startsWith(members, xml_decoy_prefix))
      if (drop_decoys && !is.null(xml_decoy_prefix) && nzchar(xml_decoy_prefix))
        members <- members[!startsWith(members, xml_decoy_prefix)]
      group_id <- XML::xmlGetAttr(group, "group_number", default = "")
      for (protein in proteins) {
        peptides <- children(protein, "peptide")
        for (peptide in peptides) {
          variants <- c(list(peptide), children(peptide, "indistinguishable_peptide"))
          for (variant in variants) {
            sequence <- XML::xmlGetAttr(variant, "peptide_sequence")
            rows[[length(rows) + 1L]] <- data.frame(Protein = paste0("ProteinProphet_group_", group_id),
              Protein_members = paste(members[nzchar(members)], collapse = ";"),
              Peptide = sequence, Modification_detail = modification_detail(variant, sequence),
              Decoy = if (is_decoy) "true" else "",
              stringsAsFactors = FALSE)
          }
        }
      }
    }
    if (!length(rows)) stop("No peptide records found in protXML.", call. = FALSE)
    do.call(rbind, rows)
  }
  xml_source_format <- NULL
  if (is.character(input) && length(input) == 1L) {
    if (!file.exists(input)) stop("Input file does not exist: ", input, call. = FALSE)
    is_xml <- grepl("(?:[.-](?:pep|prot)\\.xml|\\.pepxml)$", input,
      ignore.case = TRUE, perl = TRUE)
    if (format == "tpp" || (is_xml && format %in% c("auto", "fragpipe"))) {
      xml_source_format <- if (format == "fragpipe") "fragpipe" else "tpp"
      table <- tpp_table(input); format <- "tpp"
    } else table <- read_table(input)
  } else if (is.data.frame(input)) table <- input
  else stop("`input` must be a data frame or a single file path.", call. = FALSE)

  if (format == "auto") {
    if (all(c("Stripped.Sequence", "Protein.Group") %in% names(table))) format <- "diann"
    else if (any(c("Annotated Sequence", "Annotated.Sequence") %in% names(table))) format <- "proteome_discoverer"
    else if (all(c("Sequence", "Modified sequence") %in% names(table))) format <- "maxquant"
    else if (any(c("Peptide Sequence", "Assigned Modifications") %in% names(table)) ||
             all(c("Peptide", "Protein") %in% names(table))) format <- "fragpipe"
    else stop("Could not detect the producer. Set `format` explicitly.", call. = FALSE)
  }
  make_standard <- function(x, format) {
    column <- function(options, required = TRUE) first_column(x, options, required)
    make <- function(protein, members, peptide, detail, decoy = NULL) data.frame(
      Protein = x[[protein]], Protein_members = if (is.null(members)) "" else x[[members]],
      Peptide = x[[peptide]], Modification_detail = if (is.null(detail)) "" else x[[detail]],
      Decoy = if (is.null(decoy)) "" else x[[decoy]], stringsAsFactors = FALSE)
    if (format == "diann") return(make(column(c("Protein.Group", "Protein.Ids")),
      column(c("Protein.Ids", "Protein.Group"), FALSE), column("Stripped.Sequence"),
      column("Modified.Sequence", FALSE), column("Decoy", FALSE)))
    if (format == "fragpipe") {
      modification <- column(c("Assigned Modifications", "Modified Sequence", "Modified Peptide"), FALSE)
      result <- make(column(c("Protein", "Protein ID", "Protein.Ids", "Mapped Proteins")),
        column(c("Mapped Proteins", "Protein", "Protein ID"), FALSE),
        column(c("Peptide Sequence", "Peptide", "Stripped.Sequence")), modification,
        column(c("Decoy", "Is Decoy"), FALSE))
      result$Protein_members <- gsub(",\\s*", ";", result$Protein_members)
      result$Modification_is_absolute <- identical(modification, "Modified Peptide")
      return(result)
    }
    if (format == "maxquant") {
      result <- make(column(c("Leading razor protein", "Leading proteins", "Proteins")),
        column(c("Proteins", "Leading proteins"), FALSE), column("Sequence"),
        column(c("Modifications", "Modified sequence"), FALSE), column("Reverse", FALSE))
      if (drop_contaminants && "Potential contaminant" %in% names(x))
        result <- result[!tolower(trimws(x[["Potential contaminant"]])) %in%
          c("+", "1", "true", "yes"), , drop = FALSE]
      return(result)
    }
    if (format == "proteome_discoverer") return(make(
      column(c("Master Protein Accessions", "Master.Protein.Accessions", "Protein Accessions", "Protein.Accessions")),
      column(c("Protein Accessions", "Protein.Accessions", "Master Protein Accessions"), FALSE),
      column(c("Annotated Sequence", "Annotated.Sequence", "Sequence")),
      column(c("Modifications", "Modification"), FALSE),
      column(c("Is Decoy", "Decoy", "Decoy PSM"), FALSE)))
    if (format == "tpp") return(x)
  }
  standard <- make_standard(table, format)
  for (field in c("Protein_members", "Decoy")) {
    if (!field %in% names(standard)) standard[[field]] <- ""
  }
  if (!"Modification_is_absolute" %in% names(standard)) standard$Modification_is_absolute <- format == "tpp"
  standard[] <- lapply(standard, as.character)
  if (drop_decoys) standard <- standard[!(tolower(trimws(standard$Decoy)) %in%
    c("1", "+", "true", "yes", "decoy")), , drop = FALSE]
  if (format == "fragpipe" && is.null(xml_source_format)) {
    if (is.null(fragger_params)) stop("FragPipe tables require `fragger_params` to account for fixed modifications.", call. = FALSE)
    settings <- readLines(fragger_params, warn = FALSE)
    settings <- trimws(sub("#.*$", "", settings))
    fixed <- settings[grepl("^add_(?:[A-Z]_|[NC]term_)", settings)]
    fixed <- strsplit(fixed, "=", fixed = TRUE)
    for (setting in fixed) {
      key <- trimws(setting[[1]])
      value <- as.numeric(trimws(setting[[2]]))
      if (!is.finite(value) || value == 0) next
      if (grepl("_protein$", key)) stop("Fixed protein-terminal FragPipe modification ",
        key, " requires protein-boundary mapping before import.", call. = FALSE)
      site <- if (grepl("^add_[A-Z]_", key)) substr(key, 5L, 5L) else
        if (grepl("Nterm", key)) "n" else "c"
      map_key <- sprintf("%s:%+.4f", site, value)
      id <- mapped_id(map_key)
      if (is.null(id)) stop("Fixed FragPipe modification ", map_key,
        " needs an exact UniMod ID in `modification_map`.", call. = FALSE)
      for (row in seq_len(nrow(standard))) {
        count <- if (site %in% c("n", "c")) 1L else
          lengths(regmatches(standard$Peptide[[row]], gregexpr(site, standard$Peptide[[row]], fixed = TRUE)))
        if (site %in% LETTERS && count) {
          annotated <- regmatches(standard$Modification_detail[[row]],
            gregexpr(paste0("[0-9]+", site, "\\([+-]?[0-9.]+\\)"),
              standard$Modification_detail[[row]], perl = TRUE))[[1]]
          if (length(annotated)) {
            shifts <- as.numeric(sub(".*\\(([+-]?[0-9.]+)\\)", "\\1", annotated))
            count <- count - sum(abs(shifts - value) < 0.001)
          }
        }
        if (count) standard$Modification_detail[[row]] <- paste(
          standard$Modification_detail[[row]], paste(rep(paste0("UniMod:", id), count), collapse = ";"), sep = ";")
      }
    }
  }
  if (format == "maxquant") {
    if (is.null(mqpar)) stop("MaxQuant evidence requires `mqpar` to account for fixed modifications.", call. = FALSE)
    if (!requireNamespace("XML", quietly = TRUE)) stop("MaxQuant mqpar import requires XML.", call. = FALSE)
    config <- XML::xmlParse(mqpar)
    fixed <- XML::getNodeSet(config, "//*[local-name()='fixedModifications']/*[local-name()='string']")
    for (node in fixed) {
      label <- XML::xmlValue(node)
      site <- sub(".*\\(([A-Z])\\)$", "\\1", label)
      if (!grepl("^[A-Z]$", site)) stop("Unsupported MaxQuant fixed modification site: ", label, call. = FALSE)
      name <- sub("\\s*\\([A-Z]\\)$", "", label)
      id <- if (name %in% names(unimod_ids)) unimod_ids[[name]] else mapped_id(label)
      if (is.null(id)) stop("MaxQuant fixed modification ", label,
        " needs a UniMod ID in `modification_map`.", call. = FALSE)
      unimod_for(id)
      for (row in seq_len(nrow(standard))) {
        count <- lengths(regmatches(standard$Peptide[[row]], gregexpr(site, standard$Peptide[[row]], fixed = TRUE)))
        if (count) standard$Modification_detail[[row]] <- paste(
          standard$Modification_detail[[row]], paste(rep(paste0("UniMod:", id), count), collapse = ";"), sep = ";")
      }
    }
  }
  empty_candidates <- function() data.frame(
    Protein = character(), Protein_members = character(), Peptide = character(),
    Modification = character(), formula = character(), neutral_mass = numeric(),
    pepmz = numeric(), mz = numeric(), adduct = character(), charge = integer(),
    isdecoy = integer(), moleculeNames = character(), source_rows = integer(),
    source_format = character(), stringsAsFactors = FALSE)

  # UniMod is the formula authority for imported PTMs.  Keep source-format
  # aliases below, but do not duplicate UniMod elemental compositions here.
  unimod_atoms <- function(composition) {
    compact <- gsub("[[:space:]]+", "", composition)
    if (grepl("(?:^|\\))\\d+[A-Z]", compact, perl = TRUE))
      stop("Isotope-labelled UniMod composition needs isotope-aware formula scoring: ",
        composition, call. = FALSE)
    terms <- regmatches(compact, gregexpr("[A-Z][a-z]?(?:\\(-?[0-9]+\\))?", compact, perl = TRUE))[[1]]
    if (!length(terms) || !identical(paste0(terms, collapse = ""), compact)) {
      stop("Unsupported UniMod elemental composition: ", composition, call. = FALSE)
    }
    elements <- sub("([A-Z][a-z]?).*", "\\1", terms)
    counts <- vapply(terms, function(term) {
      if (!grepl("\\(", term)) return(1)
      as.numeric(sub(".*\\((-?[0-9]+)\\)", "\\1", term))
    }, numeric(1))
    stats::setNames(as.list(counts), elements)
  }

  modification_counts <- function(detail, peptide, absolute_mass = FALSE) {
    detail <- ifelse(is.na(detail), "", detail)
    if (!nzchar(detail) || tolower(detail) %in% c("unmodified", "none", "nan") ||
        identical(toupper(detail), peptide)) return(integer())
    aliases <- c(Carbamidomethyl = "carbamidomethyl", Oxidation = "oxidation",
      Acetyl = "(?<![a-z])acetyl", Phospho = "phospho",
      Deamidated = "deamidat", Methyl = "(?<![A-Za-z])methyl",
      Dimethyl = "dimethyl", Trimethyl = "trimethyl", Formyl = "formyl",
      Biotin = "biotin(?=\\s*\\[)")
    reported_unimod <- regmatches(detail, gregexpr("unimod:[0-9]+", detail, ignore.case = TRUE))[[1]]
    reported_ids <- sub(".*:", "", reported_unimod)
    for (id in reported_ids) unimod_for(id)
    answer <- integer()
    for (name in names(aliases)) {
      found <- gregexpr(aliases[[name]], detail, perl = TRUE, ignore.case = TRUE)[[1]]
      if (length(found) && found[[1]] != -1L) {
        prefixes <- substr(detail, pmax(1L, found - 6L), found - 1L)
        counts <- suppressWarnings(as.integer(sub(".*?([0-9]+)(?:\\s*x)?\\s*$", "\\1", prefixes,
          perl = TRUE)))
        counts[is.na(counts)] <- 1L
        answer[[name]] <- sum(counts)
      }
    }
    if (grepl("pyro[- ]?glu|gln.?pyro|glu.?pyro|unimod:28|unimod:27", detail, ignore.case = TRUE)) {
      if (substr(peptide, 1L, 1L) == "Q") answer[["Pyro-glu from Q"]] <- 1L
      else if (substr(peptide, 1L, 1L) == "E") answer[["Pyro-glu from E"]] <- 1L
      else stop("Pyro-glu is reported for a peptide not starting with Q/E: ", peptide, call. = FALSE)
    }
    for (id in reported_ids) {
      name <- names(unimod_ids)[match(id, unimod_ids)]
      if (is.na(name)) name <- id
      answer[[name]] <- if (name %in% names(answer)) answer[[name]] + 1L else 1L
    }
    numeric_positions <- gregexpr("(?<=[\\[\\(])[+-]?[0-9]+(?:\\.[0-9]+)?(?=[\\]\\)])",
      detail, perl = TRUE)[[1]]
    values <- if (numeric_positions[[1]] == -1L) character() else
      regmatches(detail, list(numeric_positions))[[1]]
    if (absolute_mass) {
      residue_mass <- c(A=71.037114,R=156.101111,N=114.042927,D=115.026943,C=103.009185,E=129.042593,Q=128.058578,G=57.021464,H=137.058912,I=113.084064,L=113.084064,K=128.094963,M=131.040485,F=147.068414,P=97.052764,S=87.032028,T=101.047679,W=186.079313,Y=163.063329,V=99.068414)
      tokens <- regmatches(detail, gregexpr("[ncA-Z]\\[[0-9.]+\\]", detail, perl = TRUE))[[1]]
      values <- vapply(tokens, function(token) {
        mass <- as.numeric(gsub(".*\\[|\\]", "", token))
        residue <- substr(token, 1, 1)
        if (residue == "n") mass - 1.007825 else if (residue == "c") mass - 17.002740 else mass - residue_mass[[residue]]
      }, numeric(1))
    } else values <- as.numeric(values)
    if (!length(values)) {
      if (length(answer)) return(answer[order(names(answer))])
      stop("Cannot calculate a formula for modification annotation: ", detail, call. = FALSE)
    }
    for (i in seq_along(values)) {
      value <- values[[i]]
      site <- if (absolute_mass) substr(tokens[[i]], 1L, 1L) else {
        preceding <- substr(detail, 1L, numeric_positions[[i]] - 2L)
        if (grepl("(?:N-term|Nterm)$", preceding, ignore.case = TRUE)) "n" else
          if (grepl("(?:C-term|Cterm)$", preceding, ignore.case = TRUE)) "c" else
            if (grepl("[A-Z]$", preceding)) substring(preceding, nchar(preceding)) else ""
      }
      key <- sprintf("%s:%+.4f", site, value)
      id <- mapped_id(key)
      if (is.null(id)) id <- mapped_id(sprintf("%+.4f", value))
      if (is.null(id)) stop("Modification mass shift ", signif(value, 7),
        " has no UniMod identity. Supply `modification_map` for ", key, ".", call. = FALSE)
      name <- names(unimod_ids)[match(id, unimod_ids)]
      if (is.na(name)) name <- id
      validate_mapped_mass(id, value, key)
      answer[[name]] <- if (name %in% names(answer)) answer[[name]] + 1L else 1L
    }
    answer[order(names(answer))]
  }
  formula_for <- function(peptide, modifications) {
    atoms <- ConvertPeptide(peptide, Element_tbl = BuildElement(uniqueAA = unique(strsplit(peptide, "")[[1]])))
    for (name in names(modifications)) atoms <- merge_atoms(atoms,
      unimod_atoms(unimod_for(name)$composition), mode="add",
      multiplier=c(1, modifications[[name]]), check_merge=TRUE)
    get_formula_from_atom(atoms[unlist(atoms) != 0])
  }
  if (!nrow(standard)) {
    candidates <- empty_candidates()
    if (!is.null(output_file)) utils::write.csv(candidates, output_file, row.names = FALSE)
    return(candidates)
  }
  if (any(is.na(standard$Peptide) | !nzchar(trimws(standard$Peptide)))) {
    stop("The input contains row(s) with a missing peptide sequence.", call. = FALSE)
  }
  missing_protein <- (is.na(standard$Protein) | !nzchar(trimws(standard$Protein))) &
    (is.na(standard$Protein_members) | !nzchar(trimws(standard$Protein_members)))
  if (any(missing_protein)) {
    stop("The input contains row(s) without a protein or protein-group identifier.", call. = FALSE)
  }
  standard$Peptide <- sub("^\\[?[A-Za-z]\\]?\\.([A-Za-z]+)\\.\\[?[A-Za-z]\\]?$", "\\1", standard$Peptide)
  standard$Peptide <- toupper(gsub("[^A-Za-z]", "", standard$Peptide))
  invalid <- !grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", standard$Peptide)
  if (any(invalid)) stop("Only canonical amino-acid sequences can be converted: ",
    paste(unique(standard$Peptide[invalid]), collapse=", "), call.=FALSE)
  standard$mods <- Map(modification_counts, standard$Modification_detail, standard$Peptide,
                       absolute_mass = standard$Modification_is_absolute)
  standard$Modification <- vapply(standard$mods, function(x) if (!length(x)) "" else
    paste(ifelse(x == 1L, names(x), paste0(x, "x", names(x))), collapse=";"), character(1))
  standard$formula <- mapply(formula_for, standard$Peptide, standard$mods, USE.NAMES=FALSE)
  groups <- split(seq_len(nrow(standard)), paste(standard$Peptide, standard$Modification, sep="\r"))
  candidates <- do.call(rbind, lapply(groups, function(index) {
    row <- standard[index[[1]], , drop=FALSE]
    values <- function(x) paste(unique(x[!is.na(x) & nzchar(x)]), collapse=";")
    members <- standard$Protein_members[index]
    if (!all(grepl("^ProteinProphet_group_", standard$Protein[index]))) {
      members <- c(standard$Protein[index], members)
    }
    neutral_mass <- OrgMassSpecR::MonoisotopicMass(get_atoms(row$formula), charge=0)
    data.frame(Protein=values(standard$Protein[index]), Protein_members=values(members),
      Peptide=row$Peptide, Modification=row$Modification, formula=row$formula, neutral_mass=neutral_mass,
      pepmz=neutral_mass, mz=neutral_mass + 1.007276466621, adduct="M+H", charge=1L, isdecoy=0L,
      moleculeNames=paste(row$Peptide,row$Modification,sep=";"), source_rows=length(index),
      source_format=if (is.null(xml_source_format)) format else xml_source_format,
      stringsAsFactors=FALSE)
  }))
  rownames(candidates) <- NULL
  if (!is.null(output_file)) utils::write.csv(candidates, output_file, row.names=FALSE)
  candidates
}

#' Expand imported candidates into FASTA-mapped HiTMaP scoring rows
#'
#' The LC-MS import table deliberately stores one neutral formula per peptide/PTM
#' identity. HiTMaP protein scoring instead requires one row for every peptide to
#' protein mapping, including protein-relative peptide positions. This function
#' performs that expansion against the exact FASTA used for the LC-MS search.
#' It also retains `neutral_formula` and sets `formula` to the [M+H]+ atom
#' composition expected by legacy PMF isotope scoring.
#'
#' @param candidates Output from [import_lcms_proteomics_candidates()].
#' @param database FASTA used for the LC-MS identification search.
#' @return Candidate rows suitable for `protein_scoring()`. The global
#'   `Index_of_protein_sequence` is also set for the legacy scoring function.
#' @export
prepare_lcms_candidates_for_hitmap <- function(candidates, database) {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("FASTA mapping requires Biostrings.", call. = FALSE)
  }
  required <- c("Protein", "Protein_members", "Peptide", "Modification",
                "formula", "mz", "isdecoy")
  missing <- setdiff(required, names(candidates))
  if (length(missing)) stop("Candidate table is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!is.character(database) || length(database) != 1L || !file.exists(database)) {
    stop("`database` must be an existing FASTA file.", call. = FALSE)
  }

  accession <- function(value) {
    value <- trimws(value)
    token <- strsplit(value, "\\s+")[[1]][1]
    parts <- strsplit(token, "\\|", fixed = FALSE)[[1]]
    if (length(parts) >= 2L && parts[[1]] %in% c("sp", "tr", "ref", "gb", "dbj", "emb")) parts[[2]] else token
  }
  index <- Biostrings::fasta.index(database, seqtype = "AA")
  index$recno <- seq_len(nrow(index))
  index$desc <- as.character(index$desc)
  fasta_tokens <- sub("\\s.*$", "", index$desc)
  sequence_accessions <- ifelse(
    grepl("^(sp|tr|ref|gb|dbj|emb)\\|", fasta_tokens),
    sub("^[^|]*\\|([^|]+)\\|.*$", "\\1", fasta_tokens), fasta_tokens
  )
  assign("Index_of_protein_sequence", index, envir = .GlobalEnv)

  mapped <- lapply(seq_len(nrow(candidates)), function(candidate_id) {
    candidate <- candidates[candidate_id, , drop = FALSE]
    members <- unique(unlist(strsplit(paste(candidate$Protein, candidate$Protein_members, sep = ";"), ";", fixed = TRUE)))
    members <- unique(vapply(members[nzchar(members)], accession, character(1)))
    sequence_ids <- which(sequence_accessions %in% members)
    if (!length(sequence_ids)) return(NULL)
    sequences <- Biostrings::readAAStringSet(index[sequence_ids, , drop = FALSE])
    do.call(rbind, lapply(seq_along(sequence_ids), function(position) {
      sequence_id <- sequence_ids[[position]]
      starts <- gregexpr(candidate$Peptide, as.character(sequences[[position]]), fixed = TRUE)[[1]]
      if (identical(starts, -1L)) return(NULL)
      result <- candidate[rep(1L, length(starts)), , drop = FALSE]
      result$Candidate_ID <- candidate_id
      result$Protein_source <- candidate$Protein
      result$Protein <- index$recno[[sequence_id]]
      result$start <- as.integer(starts)
      result$end <- as.integer(starts) + nchar(candidate$Peptide) - 1L
      result$pro_end <- nchar(as.character(sequences[[position]]))
      result$desc <- index$desc[[sequence_id]]
      result
    }))
  })
  mapped <- Filter(Negate(is.null), mapped)
  if (!length(mapped)) {
    stop("No imported peptide/protein identifiers could be mapped to the supplied FASTA.", call. = FALSE)
  }
  result <- do.call(rbind, mapped)
  rownames(result) <- NULL
  result$neutral_formula <- result$formula
  result$formula <- vapply(result$neutral_formula, function(formula)
    get_formula_from_atom(merge_atoms(get_atoms(formula), list(H = 1))), character(1))
  mapped_ids <- unique(result$Candidate_ID)
  attr(result, "unmapped_candidates") <- candidates[setdiff(seq_len(nrow(candidates)), mapped_ids), , drop = FALSE]
  result
}
