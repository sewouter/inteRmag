#' @title Convert an .AZ file to a pmob
#'
#' @description Convert an .AZ file (from the Centro Interuniversitario di
#' Magnetismo Naturale "Roberto Lanza" - Alpine Laboratory of Paleomagnetism;
#' CIMAN-ALP format) into a pmob (\strong{P}aleo\strong{M}agnetic
#' \strong{OB}ject)
#'
#' @param file the directory of the file to convert
#'
#' @examples
#' file <- system.file("tests",
#'                     "CIMAN-ALP_example.AZ",
#'                     package = "inteRmag")
#'
#' pmob_from_AZ(file)
#'
#' @export
#' @importFrom utils read.table

pmob_from_AZ <- function(file){
  tab <- read.table(file, sep = "")
  colnames(tab) <- c("id", "az", "dip", "bedaz", "beddip")

  res <- as.pmob(sampid = tab$id,
                 sampaz = tab$az,
                 sampdip = tab$dip,
                 samprot = 0,
                 bedaz = tab$bedaz,
                 beddip = tab$beddip)
  return(res)
}
