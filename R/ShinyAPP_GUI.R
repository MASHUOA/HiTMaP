#' HiTMaP shiny-based graphical user interface (GUI)
#' 
#'
#' @docType function
#' @name HiTMaP_GUI
#' @rdname HiTMaP_GUI
#'
#' @examples
#' ## Start the GUI
#' \dontrun{
#' HiTMaP_GUI()
#' }
#'
#'

#'
HiTMaP_GUI <- function(wd="~/", port = 3838) {
  
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("The HiTMaP GUI requires the 'shiny' package. Install it before starting the GUI.", call. = FALSE)
  }
  appDir <- system.file("HiTMaP_GUI", package = "HiTMaP")
  if (appDir == "") {
    stop("Could not find GUI directory. Try re-install `HiTMaP`.",
         call. = FALSE
    )
  }
  if (wd=="") wd=appDir
  
  WorkingDir_global<<-wd
  message(WorkingDir_global)
  shiny::runApp(appDir, port = port)
}
