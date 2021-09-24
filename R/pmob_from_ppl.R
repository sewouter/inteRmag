#' @title Convert an .ppl file to a pmob
#'
#' @description Convert .ppl file (from the PuffinPlot format) into a pmob
#' (\strong{P}aleo\strong{M}agnetic \strong{OB}ject)
#'
#' @param file the directory of the file to convert
#' @param depth.factor the factor by which the depth values in the .ppl file
#' need to be multiplied to convert to meters (e.g. if the units are
#' centimeters, the depth.factor parameter needs to set at 100). By default the
#' function assumes that the units are in meters (depth.factor = 1).
#' @param time.format how the time in the TIMESTAMP column of the .ppl file is
#' formatted. For more details, check the \code{strptime} function in R base.
#' @param magsus.type whether the magnetic susceptibility units are normalised
#' by volume (magsus.type = "vol"), by mass (magsus.type = "mass"), or are not
#' normalised (i.e. total magnetic susceptibility; magsus.type = "tot"). By
#' default the units are considered to be normalised by volume.
#' @param magsus.factor the factor by which the magnetic susceptibility values
#' in the .ppl file need to be multiplied to convert to the proper magnetic
#' susceptibility units. If the magnetic susceptibility is normalised by volume,
#' the values need to be converted into dimensionless SI units, if normalised
#' by mass, into cubic meters by kilograms, and if not normalised (total
#' magnetic susceptibility), into cubic meters. By default magsus.factor is set
#' to 1.
#' @param irm.axis the axis along which the Isothermal Remanent Magnetisation
#' (IRM) was applied (can be "x", "y", "z", "-x", "-y" or "-z")
#' @param arm.bias.axis the axis along which the bias field of the Anhysteretic
#' Remanent Magnetisation (ARM) was applied (can be "x", "y", "z", "-x", "-y" or
#' "-z")
#' @param defaults a list of the default values of the columns in the .ppl file
#' (they will be translated to NA values; or in the case of the volume and area,
#' they will be signaled as being default/not exact values).
#'
#' @examples
#' file <- system.file("tests",
#'                     "PuffinPlot_example.ppl",
#'                     package = "inteRmag")
#'
#' pmob_from_ppl(file)
#'
#' @export
#' @importFrom StratigrapheR merge_list

pmob_from_ppl <- function(file,
                          depth.factor = 1,
                          time.format = "%m/%d/%y %H%M",
                          magsus.type = c("vol", "mass", "tot"),
                          magsus.factor = 1,
                          irm.axis = c("x", "y", "z", "-x", "-y", "-z"),
                          arm.bias.axis = c("x", "y", "z", "-x", "-y", "-z"),
                          defaults = list(DISCRETE_ID = "UNSET",
                                          DEPTH = "null",
                                          RUN_NUMBER = "-1",
                                          TIMESTAMP = "UNSET",
                                          SLOT_NUMBER = "-1",
                                          MEAS_TYPE = "UNSET",
                                          X_MOMENT = NULL,
                                          Y_MOMENT = NULL,
                                          Z_MOMENT = NULL,
                                          MAG_SUS = "NaN",
                                          VOLUME = "10.8",
                                          AREA = "4.0",
                                          SAMPLE_AZ = "NaN",
                                          SAMPLE_DIP = "NaN",
                                          FORM_AZ = "NaN",
                                          FORM_DIP = "NaN",
                                          MAG_DEV = "0.0",
                                          TREATMENT = "UNKNOWN",
                                          AF_X = "NaN",
                                          AF_Y = "NaN",
                                          AF_Z = "NaN",
                                          TEMPERATURE = "NaN",
                                          IRM_FIELD = "NaN",
                                          ARM_FIELD = "NaN",
                                          ARM_AXIS = "UNKNOWN",
                                          PP_SELECTED = NULL,
                                          PP_ANCHOR_PCA = NULL,
                                          PP_HIDDEN = NULL,
                                          PP_ONCIRCLE = NULL,
                                          PP_INPCA = NULL))
{

  ppt <- read.ppl(file)

  # Default values ----

  defaults <- merge_list(defaults,
                         list(DISCRETE_ID = "UNSET",
                              DEPTH = "null",
                              RUN_NUMBER = "-1",
                              TIMESTAMP = "UNSET",
                              SLOT_NUMBER = "-1",
                              MEAS_TYPE = "UNSET",
                              X_MOMENT = NULL,
                              Y_MOMENT = NULL,
                              Z_MOMENT = NULL,
                              MAG_SUS = "NaN",
                              VOLUME = "10.8",
                              AREA = "4.0",
                              SAMPLE_AZ = "NaN",
                              SAMPLE_DIP = "NaN",
                              FORM_AZ = "NaN",
                              FORM_DIP = "NaN",
                              MAG_DEV = "0.0",
                              TREATMENT = "UNKNOWN",
                              AF_X = "NaN",
                              AF_Y = "NaN",
                              AF_Z = "NaN",
                              TEMPERATURE = "NaN",
                              IRM_FIELD = "NaN",
                              ARM_FIELD = "NaN",
                              ARM_AXIS = "UNKNOWN",
                              PP_SELECTED = NULL,
                              PP_ANCHOR_PCA = NULL,
                              PP_HIDDEN = NULL,
                              PP_ONCIRCLE = NULL,
                              PP_INPCA = NULL))

  defaults <- defaults[!sapply(defaults, is.null)]

  mtch <- match(names(defaults), colnames(ppt))

  y <- as.list(ppt)

  for(i in seq(length(defaults))){

    y[[mtch[i]]][y[[mtch[i]]] == defaults[i]] <- NA

  }

  ppt <- as.data.frame(y)

  # Time of measurement ----

  ntime <- strptime(ppt$TIMESTAMP, format = time.format)

  meassec   <- as.integer(format(ntime, "%S"))
  measmin   <- as.integer(format(ntime, "%M"))
  meashour  <- as.integer(format(ntime, "%H"))
  measday   <- as.integer(format(ntime, "%d"))
  measmonth <- as.integer(format(ntime, "%m"))
  measyear  <- as.integer(format(ntime, "%Y"))

  # Discrete or continuous measurement ----

  discrete <- rep(NA, nrow(ppt))

  discrete[ppt$MEAS_TYPE == "DISCRETE"]   <- T
  discrete[ppt$MEAS_TYPE == "CONTINUOUS"] <- F

  # Volume exactitude ----

  vol <- as.numeric(ppt$VOLUME)

  exactvol <- rep(T, nrow(ppt))

  if(any(is.na(ppt$VOLUME))){

    defvol <- as.numeric(defaults$VOLUME)

    if(is.na(defvol)) {
      stop("Default volume cannot be converted into numerical value")
    }

    exactvol[is.na(ppt$VOLUME)] <- F

    vol[is.na(ppt$VOLUME)] <- defvol

  }

  vol <- vol * 10^-6

  # Area exactitude ----

  area <- as.numeric(ppt$AREA)

  defaultarea <- rep(F, nrow(ppt))

  if(any(is.na(ppt$AREA))){

    defarea <- as.numeric(defaults$AREA)

    if(is.na(defarea)) {
      stop("Default area cannot be converted into numerical value")
    }

    defaultarea[is.na(ppt$AREA)] <- T

    area[is.na(ppt$AREA)] <- defarea

  }

  area <- area * 10^-4

  # Magnetic Deviation ----

  magvar <- as.numeric(ppt$MAG_DEV)

  if(any(is.na(magvar))){

    defmagvar <- as.numeric(defaults$MAG_DEV)

    if(is.na(defmagvar)) {
      stop("Default magnetic deviation/variation cannot be converted",
           " into numerical value")
    }

    magvar[is.na(ppt$MAG_DEV)] <- defmagvar

  }

  sampaz   <- rep(NA, nrow(ppt))

  magaz       <- rep(NA, nrow(ppt))
  usemagaz    <- rep(NA, nrow(ppt))
  magazvarcor <- rep(NA, nrow(ppt))
  bedazvarcor <- rep(NA, nrow(ppt))

  sampaz[is.na(ppt$MAG_DEV)] <- ppt$SAMPLE_AZ[is.na(ppt$MAG_DEV)]

  magaz[!is.na(ppt$MAG_DEV)]       <- ppt$SAMPLE_AZ[!is.na(ppt$MAG_DEV)]
  usemagaz[!is.na(ppt$MAG_DEV)]    <- T
  magazvarcor[!is.na(ppt$MAG_DEV)] <- T
  bedazvarcor[!is.na(ppt$MAG_DEV)] <- T

  # Treatment ----

  # Thermal

  treattemp <- as.numeric(ppt$TEMPERATURE) + 273.15

  # Conflicts

  mult.treat <- as.integer(!is.na(ppt$AF_X) | !is.na(ppt$AF_Y) | !is.na(ppt$AF_Z)) +
    as.integer(!is.na(ppt$TEMPERATURE)) +
    as.integer(!is.na(ppt$IRM_FIELD))

  if(any(mult.treat > 1)){

    warning("Multiple treatment methodologies were applied to unique samples")

  }

  af.pos1 <- ppt$TREATMENT == "DEGAUSS_XYZ" &
    (is.na(ppt$AZ_X) | is.na(ppt$AZ_Y)  | is.na(ppt$AZ_Z))

  if(any(af.pos1)){

    warning("Some AF information is lacking to have a XYZ Degaussing procedure",
            " as expected form the type of treatment documented")

  }

  af.pos2 <- ppt$TREATMENT == "DEGAUSS_Z" & is.na(ppt$AZ_Z)

  if(any(af.pos2)){

    warning("Some AF information is lacking to have a Z Degaussing procedure",
            " as expected form the type of treatment documented")

  }

  th.pos <- ppt$TREATMENT == "THERMAL" & is.na(ppt$TEMPERATURE)

  if(any(th.pos)){

    warning("Some temperature information is lacking")

  }

  irm.pos <- ppt$TREATMENT == "IRM" & is.na(ppt$IRM_FIELD)

  if(any(irm.pos)){

    warning("Some IRM information is lacking")

  }

  arm.pos <- ppt$TREATMENT == "ARM" &
    (is.na(ppt$ARM_FIELD) | is.na(ppt$ARM_AXIS) |
       is.na(ppt$AZ_X) | is.na(ppt$AZ_Y)  |  is.na(ppt$AZ_Z))

  if(any(arm.pos)){

    warning("Some ARM information is lacking")

  }

  # IRM ----

  if(!is.null(ppt$IRM_FIELD)){

    if(irm.axis[1] == "x") {

      treatirmx <- as.numeric(ppt$IRM_FIELD)
      treatirmy <- NULL
      treatirmz <- NULL

    } else if(irm.axis[1] == "y"){

      treatirmx <- NULL
      treatirmy <- as.numeric(ppt$IRM_FIELD)
      treatirmz <- NULL

    } else if(irm.axis[1] == "z"){

      treatirmx <- NULL
      treatirmy <- NULL
      treatirmz <- as.numeric(ppt$IRM_FIELD)

    } else if(irm.axis[1] == "-x"){

      treatirmx <- -as.numeric(ppt$IRM_FIELD)
      treatirmy <- NULL
      treatirmz <- NULL

    } else if(irm.axis[1] == "-y"){

      treatirmx <- NULL
      treatirmy <- -as.numeric(ppt$IRM_FIELD)
      treatirmz <- NULL

    } else if(irm.axis[1] == "-z"){

      treatirmx <- NULL
      treatirmy <- NULL
      treatirmz <- -as.numeric(ppt$IRM_FIELD)

    } else {

      stop("irm.axis should be 'x', 'y', 'z', '-x', '-y', or '-z'")

    }

  } else {

    treatirmx <- NULL
    treatirmy <- NULL
    treatirmz <- NULL

  }


  # AF or ARM ----

  if(!is.null(ppt$ARM_FIELD)){

    treatarmafx <- as.numeric(ppt$AF_X)
    treatarmafy <- as.numeric(ppt$AF_Y)
    treatarmafz <- as.numeric(ppt$AF_Z)

    if(arm.bias.axis[1] == "x") {

      treatarmbiasx <- ppt$ARM_FIELD
      treatarmbiasy <- NULL
      treatarmbiasz <- NULL

    } else if(arm.bias.axis[1] == "y"){

      treatarmbiasx <- NULL
      treatarmbiasy <- ppt$ARM_FIELD
      treatarmbiasz <- NULL

    } else if(arm.bias.axis[1] == "z"){

      treatarmbiasx <- NULL
      treatarmbiasy <- NULL
      treatarmbiasz <- ppt$ARM_FIELD

    } else if(arm.bias.axis[1] == "-x"){

      treatarmbiasx <- -ppt$ARM_FIELD
      treatarmbiasy <- NULL
      treatarmbiasz <- NULL

    } else if(arm.bias.axis[1] == "-y"){

      treatarmbiasx <- NULL
      treatarmbiasy <- -ppt$ARM_FIELD
      treatarmbiasz <- NULL

    } else if(arm.bias.axis[1] == "-z"){

      treatarmbiasx <- NULL
      treatarmbiasy <- NULL
      treatarmbiasz <- -ppt$ARM_FIELD

    } else {

      stop("arm.bias.axis should be 'x', 'y', 'z', '-x', '-y', or '-z'")

    }

    treatafx <- NULL
    treatafy <- NULL
    treatafz <- NULL

  } else {

    treatarmafx <- NULL
    treatarmafy <- NULL
    treatarmafz <- NULL

    treatarmbiasx <- NULL
    treatarmbiasy <- NULL
    treatarmbiasz <- NULL

    treatafx <- as.numeric(ppt$AF_X)
    treatafy <- as.numeric(ppt$AF_Y)
    treatafz <- as.numeric(ppt$AF_Z)

  }

  # Magnetic susceptibility ----

  magsus.unit <- magsus.factor * as.numeric(ppt$MAG_SUS)

  if(magsus.type[1] == "vol"){

    volmagsus  <- magsus.unit
    massmagsus <- NULL
    totmagsus  <- NULL

  } else if(magsus.type[1] == "mass"){

    volmagsus  <- NULL
    massmagsus <- magsus.unit
    totmagsus  <- NULL


  } else if(magsus.type[1] == "tot"){

    volmagsus  <- NULL
    massmagsus <- NULL
    totmagsus  <- magsus.unit

  } else {

    stop("magsus.type should be 'vol', 'mass', or 'tot'")

  }


  # PCA and Great Circles ----

  pcacomponent                           <- rep("P0", nrow(ppt))
  pcacomponent[as.logical(ppt$PP_INPCA)] <- "P1"
  pcacomponentsingle                     <-rep("P1", nrow(ppt))

  by_sample_pca <- split(pcacomponent, ppt$DISCRETE_ID)

  pca_out <- lapply(by_sample_pca, function(x) all(x == "PO"))
  pca_na  <- unsplit(pca_out, ppt$DISCRETE_ID)

  pcacomponent[pca_na]       <- NA
  pcacomponentsingle[pca_na] <- NA

  circlecomponent                              <- rep("C0", nrow(ppt))
  circlecomponent[as.logical(ppt$PP_ONCIRCLE)] <- "C1"
  circlecomponentsingle                        <- rep("C1", nrow(ppt))

  by_sample_circle <- split(circlecomponent, ppt$DISCRETE_ID)

  circle_out <- lapply(by_sample_circle, function(x) all(x == "C0"))
  circle_na  <- unsplit(circle_out, ppt$DISCRETE_ID)

  circlecomponent[circle_na]       <- NA
  circlecomponentsingle[circle_na] <- NA


  # Divide columns in the ones convertibles into explicit pmob ones, and others

  treated <- c("DISCRETE_ID",
               "DEPTH",
               "RUN_NUMBER",
               "TIMESTAMP",
               "SLOT_NUMBER",
               "MEAS_TYPE",
               "X_MOMENT",
               "Y_MOMENT",
               "Z_MOMENT",
               "VOLUME",
               "AREA",
               "SAMPLE_AZ",
               "SAMPLE_DIP",
               "FORM_AZ",
               "FORM_DIP",
               "MAG_DEV",
               "AF_X",
               "AF_Y",
               "AF_Z",
               "TEMPERATURE",
               "PP_ANCHOR_PCA",
               "PP_ONCIRCLE",
               "PP_INPCA")

  in.treated <- match(colnames(ppt), treated)

  others <- as.list(ppt[,which(is.na(in.treated))])

  # Generate pmob object ----

  out <- as.pmob(
    others,
    sampid = ppt$DISCRETE_ID,
    slotid = ppt$SLOT_NUMBER,
    measid = ppt$RUN_NUMBER,
    depth = as.numeric(ppt$DEPTH) * depth.factor,
    meassec = meassec,
    measmin =measmin,
    meashour = meashour,
    measday = measday,
    measmonth = measmonth,
    measyear = measyear,
    discrete = discrete,
    xvol = as.numeric(ppt$X_MOMENT),
    yvol = as.numeric(ppt$Y_MOMENT),
    zvol = as.numeric(ppt$Z_MOMENT),
    vol = vol,
    exactvol = exactvol,
    area = area,
    defaultarea = defaultarea,
    sampaz = as.numeric(sampaz),
    sampdip = as.numeric(ppt$SAMPLE_DIP),
    samprot = 0,
    bedaz = as.numeric(ppt$FORM_AZ),
    beddip = as.numeric(ppt$FORM_DIP),
    magaz = as.numeric(magaz),
    usemagaz = usemagaz,
    magvar = magvar,
    magazvarcor = magazvarcor,
    bedazvarcor = bedazvarcor,
    treatafx = treatafx,
    treatafy = treatafx,
    treatafz = treatafx,
    treattemp = treattemp,
    treatarmafx	= treatarmafx,
    treatarmafy = treatarmafy,
    treatarmafz = treatarmafz,
    treatarmbiasx = treatarmbiasx,
    treatarmbiasy = treatarmbiasy,
    treatarmbiasz = treatarmbiasz,
    pcaanchor = as.logical(ppt$PP_ANCHOR_PCA),
    pcacomponent = pcacomponent,
    pcacomponentsingle = pcacomponentsingle,
    circlecomponent = circlecomponent,
    circlecomponentsingle = circlecomponentsingle,
    totmagsus = totmagsus,
    volmagsus = volmagsus,
    massmagsus = massmagsus
  )


  rem <- apply(out, 2, function(x) all(is.na(x)))

  pm1 <- out[,!rem]

  return(pm1)

}






