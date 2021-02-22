#' HiTMaP shiny-based graphical user interface (GUI)
#' 
#'
#' @docType function
#' @name shinyWYSIWYG-GUI
#' @rdname shinyWYSIWYG-GUI
#'
#' @examples
#' ## Start the GUI
#' \dontrun{
#' shinyWYSIWYG()
#' }
#'
#' @import shiny
#'
#' @export HiTMaP_GUI
#'
HiTMaP_GUI <- function() {
  appDir <- system.file("shiny", "HiTMaP_GUI", package = "HiTMaP")
  if (appDir == "") {
    stop("Could not find GUI directory. Try re-install `HiTMaP`.",
         call. = FALSE
    )
  }
  
  runApp(appDir,port = 3838, launch.browser = T )
}