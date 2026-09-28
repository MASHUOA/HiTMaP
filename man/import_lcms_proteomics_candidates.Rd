% Generated manually from R/z_lcms_proteomics_import.R
\name{import_lcms_proteomics_candidates}
\alias{import_lcms_proteomics_candidates}
\title{Import LC-MS proteomics identifications as MALDI-IMS peptide candidates}
\usage{
import_lcms_proteomics_candidates(
  input,
  format = c("auto", "diann", "fragpipe", "maxquant", "proteome_discoverer", "tpp"),
  drop_decoys = TRUE,
  output_file = NULL,
  fdr_threshold = 0.01,
  probability_cutoff = NULL,
  modification_map = NULL,
  modification_mass_tolerance = 0.001,
  fragger_params = NULL,
  pepxml = NULL,
  mqpar = NULL,
  decoy_prefix = NULL,
  drop_contaminants = TRUE
)
}
\arguments{
\item{input}{A data frame, delimited export, DIA-NN parquet report, TPP or
FragPipe pepXML (for example \code{interact.pep.xml}), or ProteinProphet protXML
(for example \code{interact.prot.xml}) file.}
\item{format}{The input producer, or \code{"auto"} to detect it from columns.}
\item{drop_decoys}{Exclude marked decoys and XML accessions matching the run's decoy prefix.}
\item{output_file}{Optional CSV destination for candidates.}
\item{fdr_threshold}{Maximum model-estimated FDR for TPP or FragPipe XML input.
Defaults to 0.01; set to \code{NULL} to skip XML confidence filtering.
Tabular exports are assumed to have been filtered by their producer.}
\item{probability_cutoff}{Optional probability threshold calibrated externally
for the requested XML FDR. Use when the XML error curve is absent.}
\item{modification_map}{Named character vector mapping source annotations to exact
UniMod IDs, such as \code{c("M:+15.9949" = "35")}. Unidentified mass shifts are rejected.}
\item{modification_mass_tolerance}{Maximum difference in Da between an explicitly
mapped source shift and UniMod mass. Increase only for an imprecise legacy
search setting; the candidate formula and mass use UniMod chemistry.}
\item{fragger_params}{MSFragger \code{fragger.params} path, required for FragPipe tables
to account for fixed modifications.}
\item{pepxml}{Accompanying pepXML path, required for protXML search modifications.}
\item{mqpar}{MaxQuant \code{mqpar.xml} path, required for evidence fixed modifications.}
\item{decoy_prefix}{Optional decoy protein accession prefix for XML input.}
\item{drop_contaminants}{Exclude MaxQuant rows marked \code{Potential contaminant}.}
}
\value{A collapsed peptide/PTM candidate data frame. \code{formula} and
\code{neutral_mass} are neutral; \code{mz} is calculated as [M+H]+ for
MALDI-MS1 matching.}
\description{Imports conventional global LC-MS proteomics identification
outputs and calculates neutral peptide formulas. For TPP, pepXML preserves
peptide-level evidence and protXML provides ProteinProphet protein-group
membership. PTM localisations and LC-MS charge states are intentionally not
retained. The base peptide formula is calculated directly and PTM formula
composition and mass validation use HiTMaP's packaged UniMod data; this does
not invoke legacy candidate generation. For XML, the FDR filter uses the
embedded PeptideProphet/iProphet error curve for pepXML or ProteinProphet
error table for protXML. These represent PSM and protein-group confidence,
respectively; neither guarantees 1\% FDR for the collapsed peptide/PTM list.
XML without a usable curve requires an externally calibrated
\code{probability_cutoff} or \code{fdr_threshold = NULL}. Numeric mass shifts
require explicit UniMod identities. Search modification declarations are
checked against the mapped UniMod masses. Isotope-labelled compositions are
rejected until the formula and isotope-scoring path supports labeled atoms.}
\examples{
\dontrun{
candidates <- import_lcms_proteomics_candidates("report.tsv", format = "diann")
}
}
