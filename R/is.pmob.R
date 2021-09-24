#' @title Check if object is a pmob
#'
#' @description Check if object is a pmob (\strong{P}aleo\strong{M}agnetic
#' \strong{OB}ject).
#'
#' @param pmob Object to check.
#'
#' @examples
#' pmob <- as.pmob(sampid = "Bosso1", exc = 1, measmin = NA)
#'
#' is.pmob(pmob)
#'
#' @export
#' @import lubridate


is.pmob <- function(pmob)
{

  # exit.fun <- function(x) if(!x) print(x)
  # exit.fun <- function(x) if(!x) return(x)

  if(!inherits(pmob, "data.frame")) {
    warning("To be a pmob, the object should be of class 'data.frame'")
    return(FALSE)
  }

  cn <- colnames(pmob)

  if(!any(cn %in% "initialorder")){
    warning("To be a pmob, the object should have a column 'initial order'")
    return(FALSE)
  }

  if(!any(cn %in% "version")){
    warning("To be a pmob, the object should have a column 'version'")
    return(FALSE)
  }

  # Elements class check ----

  cn.char <- c("sampid", "specid", "slotid","measid", "siteid",
               "measdevice", "magsusdevice",
               "blankmethod", "holdermethod",
               "pcacomponent", "pcacomponentsingle",
               "circlecomponent","circlecomponentsingle",
               "comments", "originalfile",
               "version")

  cn.int <- c("initialorder",
              "measmin", "meashour",
              "measday", "measmonth", "measyear",
              "samphour", "sampday",
              "sampmonth", "sampyear",
              "samptimezonemin", "samptimezonehour",
              "originalfileline")

  cn.num <- c("depth",
              "meassec",
              "xint", "yint", "zint",
              "xinterror", "yinterror", "zinterror",
              "xvol", "yvol", "zvol",
              "xmass", "ymass", "zmass",
              "totmagsus", "totmagsuserror",
              "volmagsus","massmagsus",
              "vol", "mass", "dens",
              "area",
              "xefflength", "yefflength", "zefflength",
              "xeffvol", "yeffvol", "zeffvol",
              "xeffmass", "yeffmass", "zeffmass",
              "magsuseffvol", "magsuseffmass",
              "xintconv", "yintconv", "zintconv",
              "sampaz", "sampdip", "samprot",
              "coraz", "cordip", "corrot",
              "bedaz", "beddip",
              "foldtrend", "foldplunge",
              "magaz",
              "solaraz", "solarangle",
              "long", "lat", "height",
              "sampmin",
              "magvar",
              "treatafx", "treatafy", "treatafz",
              "treattemp",
              "treatirmx", "treatirmy", "treatirmz",
              "treatarmafx", "treatarmafy", "treatarmafz",
              "treatarmbiasx", "treatarmbiasy", "treatarmbiasz")

  cn.log <- c("blankcor", "isblank",
              "holdercor", "isholder",
              "exactvol", "exactmass", "exactdens",
              "defaultarea",
              "discrete",
              "usemagaz",
              "magazvarcor", "bedazvarcor", "foldtrendcarcor",
              "pcaanchor")

  p.char <- pmob[,cn %in% cn.char, drop = F]
  p.int  <- pmob[,cn %in% cn.int,  drop = F]
  p.num  <- pmob[,cn %in% cn.num,  drop = F]
  p.log  <- pmob[,cn %in% cn.log,  drop = F]

  i.char <- sapply(as.list(p.char), inherits, "character")
  i.int  <- sapply(as.list(p.int),  inherits, "integer")
  i.num  <- sapply(as.list(p.num),  inherits, "numeric")
  i.log  <- sapply(as.list(p.log),  inherits, "logical")

  class.type <- function(l, type){
    if(any(!l)) {
      prob <- names(l)[!l]
      warning("The following pmob column(s) should be of class '", type,
              "':\n - ", paste(prob, collapse = "\n - "), call. = F)
      return(FALSE)
    }

    return(TRUE)

  }

  if(!class.type(i.char, "character")) return(F)
  if(!class.type(i.int,  "integer"))   return(F)
  if(!class.type(i.num,  "numeric"))   return(F)
  if(!class.type(i.log,  "logical"))   return(F)

  if(any(duplicated(pmob$initialorder)) & !any(is.na(pmob$initialorder))){
    warning("There should be no duplicated or NA value in the",
            " pmob column 'initialorder'")
    return(FALSE)
  }

  # Lower and upper limits check ----

  between.check <- function(arg, header,
                            xmin = NULL, min.in = TRUE,
                            xmax = NULL, max.in = TRUE)
  {
    if(!is.null(arg)){
      if(!is.null(xmin)) {

        if(isTRUE(any(arg < xmin))) {
          warning("Values in the pmob '", header,
                  "' column should be higher than ", xmin,
                  call. = FALSE)
          return(FALSE)
        }

        if(!min.in & isTRUE(any(xmin == arg))){
          warning("Values in the pmob '", header,
                  "' column should not be equal to ", xmin,
                  call. = FALSE)
          return(FALSE)
        }

      }

      if(!is.null(xmax)) {

        if(isTRUE(any(arg > xmax))) {
          warning("Values in the pmob '", header,
                  "' column should be lower than ", xmax,
                  call. = FALSE)
          return(FALSE)
        }

        if(!max.in & isTRUE(any(xmax == arg))){
          warning("Values in the pmob '", header,
                  "' column should not be equal to ", xmax,
                  call. = FALSE)
          return(FALSE)
        }
      }
    }

    return(TRUE)

  }

  if(!between.check(pmob$initialorder, "initialorder",
                    xmin = 1L)) return(F)

  if(!between.check(pmob$meassec, "meassec",
                    xmin = 0,
                    xmax = 60, max.in = FALSE)) return(F)

  if(!between.check(pmob$measmin, "measmin",
                    xmin = 0,
                    xmax = 60, max.in = FALSE)) return(F)

  if(!between.check(pmob$meashour, "meashour",
                    xmin = 0L,
                    xmax = 24L, max.in = FALSE)) return(F)

  if(!between.check(pmob$measday, "measday",
                    xmin = 0L, xmax = 31L)) return(F)

  if(!between.check(pmob$measmonth, "measmonth",
                    xmin = 0L, xmax = 12L)) return(F)

  if(!between.check(pmob$vol,  "vol",xmin = 0L)) return(F)
  if(!between.check(pmob$mass, "mass",xmin = 0L)) return(F)
  if(!between.check(pmob$dens, "dens",xmin = 0L)) return(F)

  if(!between.check(pmob$area, "area",xmin = 0L)) return(F)

  if(!between.check(pmob$xefflength , "xefflength ",xmin = 0L)) return(F)
  if(!between.check(pmob$yefflength , "yefflength ",xmin = 0L)) return(F)
  if(!between.check(pmob$zefflength , "zefflength ",xmin = 0L)) return(F)

  if(!between.check(pmob$xeffvol, "xeffvol",xmin = 0L)) return(F)
  if(!between.check(pmob$yeffvol, "yeffvol",xmin = 0L)) return(F)
  if(!between.check(pmob$zeffvol, "zeffvol",xmin = 0L)) return(F)

  if(!between.check(pmob$xeffmass, "xeffmass",xmin = 0L)) return(F)
  if(!between.check(pmob$yeffmass, "yeffmass",xmin = 0L)) return(F)
  if(!between.check(pmob$zeffmass, "zeffmass",xmin = 0L)) return(F)

  if(!between.check(pmob$magsuseffvol,  "magsuseffvol", xmin = 0L)) return(F)
  if(!between.check(pmob$magsuseffmass, "magsuseffmass",xmin = 0L)) return(F)

  if(!between.check(pmob$sampleaz, "sampaz",
                    xmin = 0,
                    xmax = 360, max.in = F)) return(F)

  if(!between.check(pmob$sampledip, "sampdip",
                    xmin = -90,
                    xmax = 90)) return(F)

  if(!between.check(pmob$samplerot, "samprot",
                    xmin = 0,
                    xmax = 360, max.in = F)) return(F)

  if(!between.check(pmob$correctionaz, "coraz",
                    xmin = 0,
                    xmax = 360, max.in = F)) return(F)

  if(!between.check(pmob$correctiondip, "cordip",
                    xmin = -90,
                    xmax = 90)) return(F)

  if(!between.check(pmob$correctionrot, "corrot",
                    xmin = 0,
                    xmax = 360, max.in = F)) return(F)

  if(!between.check(pmob$bedaz, "bedaz",
                    xmin = 0,
                    xmax = 360, max.in = F)) return(F)

  if(!between.check(pmob$beddip, "beddip",
                    xmin = 0,
                    xmax = 180)) return(F)

  if(!between.check(pmob$foldtrend, "foldtrend",
                    xmin = 0,
                    xmax = 360, max.in = F)) return(F)

  if(!between.check(pmob$foldplunge, "foldplunge",
                    xmin = 0,
                    xmax = 180)) return(F)

  if(!between.check(pmob$magaz, "magaz",
                    xmin = 0,
                    xmax = 360, max.in = F)) return(F)

  if(!between.check(pmob$solaz, "solaraz",
                    xmin = 0,
                    xmax = 360, max.in = F)) return(F)

  if(!between.check(pmob$solangle, "solarangle",
                    xmin = 0,
                    xmax = 360, max.in = F)) return(F)

  if(!between.check(pmob$long, "long",
                    xmin = -180, min.in = F,
                    xmax = 180)) return(F)

  if(!between.check(pmob$lat, "lat",
                    xmin = -90,
                    xmax = 90)) return(F)

  if(!between.check(pmob$samplingmin, "sampmin",
                    xmin = 0,
                    xmax = 60, max.in = FALSE)) return(F)

  if(!between.check(pmob$samplinghour, "samphour",
                    xmin = 0L,
                    xmax = 24L, max.in = FALSE)) return(F)

  if(!between.check(pmob$samplingday, "sampday",
                    xmin = 0L, xmax = 31L)) return(F)

  if(!between.check(pmob$samplingmonth, "sampmonth",
                    xmin = 0L, xmax = 12L)) return(F)

  if(!between.check(pmob$samplingtimezonemin,
                    "samptimezonemin",
                    xmin = 0L,
                    xmax = 60L, max.in = FALSE)) return(F)

  if(!between.check(pmob$magvar, "magvar",
                    xmin = 0,
                    xmax = 360, max.in = F)) return(F)

  if(!between.check(pmob$treattemp, "treattemp", xmin = 0L)) return(F)

  if(!between.check(pmob$originalfileline,
                    "originalfileline", xmin = 1L)) return(F)

  # Validity of measurement date ----

  if(!is.null(pmob$measday) &
     !is.null(pmob$measonth) &
     !is.null(pmob$measyear)){

    meas.rawdate <- paste(pmob$measday,
                          pmob$measmonth,
                          pmob$measyear, sep = "-")

    meas.redate <- as.Date(meas.rawdate,
                           format("%d-%m-%Y"))

    meas.na.date <- is.na(meas.redate)

    if(any(meas.na.date &
           !is.na(pmob$measday) &
           ! is.na(pmob$measmonth) &
           !is.na(pmob$measyear))){

      meas.invalid.pos <- which(meas.na.date &
                                  !is.na(pmob$measday) &
                                  ! is.na(pmob$measmonth) &
                                  !is.na(pmob$measyear))

      warning("The following measurement dates (day-month-year) ",
              "are not valid:\n",
              paste(unique(meas.rawdate[meas.invalid.pos]), collapse = "\n"),
              call. = F)
      return(FALSE)
    }
  }

  # Validity of sampling date ----

  if(!is.null(pmob$sampday) &
     !is.null(pmob$sampmonth) &
     !is.null(pmob$sampyear)){

    meas.rawdate <- paste(pmob$sampday,
                          pmob$sampmonth,
                          pmob$sampyear, sep = "-")

    meas.redate <- as.Date(meas.rawdate,
                           format("%d-%m-%Y"))

    meas.na.date <- is.na(meas.redate)

    if(any(meas.na.date &
           !is.na(pmob$sampday) &
           ! is.na(pmob$sampmonth) &
           !is.na(pmob$sampyear))){

      meas.invalid.pos <- which(meas.na.date &
                                  !is.na(pmob$sampday) &
                                  ! is.na(pmob$sampmonth) &
                                  !is.na(pmob$sampyear))

      warning("The following sampling dates (day-month-year) are not valid:\n",
              paste(unique(meas.rawdate[meas.invalid.pos]), collapse = "\n"),
              call. = F)
      return(FALSE)
    }
  }

  return(TRUE)

}


