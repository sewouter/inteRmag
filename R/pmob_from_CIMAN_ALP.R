#' @title Convert a CIMAN-ALP .txt data file to a pmob
#'
#' @description Convert a CIMAN-ALP .txt data file (from the Centro
#' Interuniversitario di Magnetismo Naturale "Roberto Lanza" - Alpine
#' Laboratory of Paleomagnetism; CIMAN-ALP format) into a pmob
#' (\strong{P}aleo\strong{M}agnetic \strong{OB}ject)
#'
#' @param file the directory of the file to convert
#' @param treat the treatment used for samples: either thermal ("TH") or
#' alternating field demagnetization ("AF")
#' @param units the units of the treatment: can be "Celsius" (default for
#' thermal demagnetization), "Kelvin",  "Millitesla" (default for alternating
#' field demagnetization) or "Tesla".
#'
#' @examples
#' file <- system.file("tests",
#'                     "CIMAN-ALP_example.txt",
#'                     package = "inteRmag")
#'
#' pmob_from_CIMAN_ALP(file, treat = "TH", units = "Celsius")
#'
#' @export
#' @importFrom stringr str_extract
#' @importFrom utils read.table

pmob_from_CIMAN_ALP <- function(file, treat = c("TH", "AF"), units = NULL)
{
  if(!is.null(units)){
    if(units != "Celsius" & units != "Kelvin" &
       units != "Tesla" & units != "Millitesla"){
      stop("The 'units' parameter should be NULL, 'Celsius', ",
           "'Kelvin', 'Tesla' or 'Millitesla'")
    }
  }

  treat <- treat[1]

  if(treat != "TH" & treat != "AF") {
    stop("The 'treat' parameter should be 'TH' or 'AF'")
  }

  if(treat == "TH" & is.null(units)) units <- "Celsius"

  if(treat == "AF" & is.null(units)) units <- "Millitesla"

  tab <- read.table(file, skip = 6)
  colnames(tab) <- c("id", "treat", "ic", "cd", "j", "cdecl", "cincl",
                     "gdecl", "gincl", "bdecl", "bincl", "susc", "v/m")

  add.info <- readLines(file, 2)[2]

  raw.lat  <- str_extract(add.info, "^LAT:[ 0-9\\.]+")
  raw.long <- str_extract(add.info, "LON:[ 0-9\\.]+")

  raw.rec <- str_extract(add.info, "REC:[ 0-9\\.]+")
  raw.sam <- str_extract(add.info, "SAM:[ 0-9\\.]+")

  if(length(raw.lat) != 1 | is.na(raw.lat))  {
    stop("Unreckognized latitude data formatting")
  }
  if(length(raw.long) != 1  | is.na(raw.long)) {
    stop("Unreckognized longitude data formatting")
  }

  if(length(raw.lat) != 1 | is.na(raw.lat))  {
    stop("Unreckognized rec data formatting")
  }
  if(length(raw.long) != 1  | is.na(raw.long)) {
    stop("Unreckognized sam data formatting")
  }

  lat  <- as.numeric(substr(raw.lat,  5, nchar(raw.lat)))
  long <- as.numeric(substr(raw.long, 6, nchar(raw.lat)))

  rec <- as.numeric(substr(raw.rec, 5, nchar(raw.rec)))
  sam <- as.numeric(substr(raw.sam, 5, nchar(raw.sam)))

  int.Am2 <- tab$j * 1e-3

  cart.int <- transphere(dec = tab$cdecl, inc = tab$cincl, int = int.Am2)

  pmob1 <- as.pmob(sampid = tab$id,
                   xint = cart.int$x,
                   yint = cart.int$y,
                   zint = cart.int$z,
                   lat = lat,
                   long = long,
                   susc = tab$susc,
                   v.m = tab$`v/m`,
                   rec = rec,
                   sam = sam)

  if(treat == "TH"){

    if(units == "Celsius"){

      treattemp <- tab$treat + 273.15

    } else if(units == "Kelvin"){

      treattemp <- tab$treat

    }

    pmob2 <- pmob.add(pmob = pmob1, treattemp = treattemp)

  } else if(treat == "AF"){

    if(units == "Tesla"){

      treataf <- tab$treat

    } else if(units == "Millitesla"){

      treataf <- tab$treat/1000

    }

    pmob2 <- pmob.add(pmob = pmob1,
                      treatafx = treataf,
                      treatafy = treataf,
                      treatafz = treataf)

  }

  return(pmob2)

}
