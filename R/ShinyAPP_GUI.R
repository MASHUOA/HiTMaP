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
  
  if(!(require(shiny))) install.packages("shiny")
  appDir <- system.file("HiTMaP_GUI", package = "HiTMaP")
  if (appDir == "") {
    stop("Could not find GUI directory. Try re-install `HiTMaP`.",
         call. = FALSE
    )
  }
  if (wd=="") wd=appDir
  
  WorkingDir_global<<-wd
  message(WorkingDir_global)
  runApp(appDir,port = port)
}