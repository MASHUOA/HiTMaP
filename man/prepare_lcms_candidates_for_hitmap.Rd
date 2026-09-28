% Generated manually from R/z_lcms_proteomics_import.R
\name{prepare_lcms_candidates_for_hitmap}
\alias{prepare_lcms_candidates_for_hitmap}
\title{Expand imported candidates into FASTA-mapped HiTMaP scoring rows}
\usage{
prepare_lcms_candidates_for_hitmap(candidates, database)
}
\arguments{
\item{candidates}{Output from \code{import_lcms_proteomics_candidates}.}
\item{database}{The FASTA used for the LC-MS identification search.}
}
\value{Rows mapped to FASTA proteins with numeric \code{Protein}, \code{start},
\code{end}, and \code{pro_end} columns. Unmapped candidates are available in
the \code{unmapped_candidates} attribute. \code{neutral_formula} preserves the
imported neutral composition; \code{formula} contains the [M+H]+ atoms used by
legacy PMF isotope scoring.}
\description{Expands collapsed LC-MS candidates to the peptide-to-protein rows
required by legacy HiTMaP protein scoring.}
