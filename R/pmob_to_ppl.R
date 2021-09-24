#' @title Convert a pmob to an .ppl file
#'
#' @description Convert a pmob (\strong{P}aleo\strong{M}agnetic \strong{OB}ject)
#' into a.ppl file (from the PuffinPlot format).
#'
#' @param pmob the data in \strong{P}aleo\strong{M}agnetic \strong{OB}ject form
#' to convert into a .ppl file
#' @param file the name of the .ppl file
#' @param time.format how the time in the TIMESTAMP column of the .ppl file
#' needs to be formatted. For more details, check the \code{strptime} function
#' in R base.
#' @param magsus.type whether the magnetic susceptibility units need to be the
#' ones normalised by volume (magsus.type = "vol"), by mass (magsus.type =
#' "mass"), or are not normalised (i.e. total magnetic susceptibility;
#' magsus.type = "tot"). By default the units are considered to be normalised
#' by volume.
#' @param other.columns if TRUE, all columns that have headers that are not
#' regular pmob or .ppl ones will be added into the .ppl file as flags (if
#' the values are "false" or "true" or as notes if the values are characters).
#' The parameter can also be a vector of all the column names to use as flags or
#' notes.
#'
#' @examples
#' file <- system.file("tests",
#'                     "PuffinPlot_example.ppl",
#'                     package = "inteRmag")
#'
#' pmob_Bosso <- pmob_from_ppl(file)
#'
#'\dontrun{
#' pmob_to_ppl(pmob_Bosso, "Bosso_Example")}
#'
#' @export

pmob_to_ppl <- function(pmob, file,
                        time.format = "%m/%d/%y %H%M",
                        magsus.type = c("none", "vol", "mass", "tot"),
                        other.columns = TRUE)
{

  if(!is.pmob(pmob)) stop("Incorrect pmob object")

  if(!is.null(pmob$isblank))  pmob <- pmob[!pmob$isblank, , drop = F]
  if(!is.null(pmob$isholder)) pmob <- pmob[!pmob$isholder, , drop = F]

  if(is.null(pmob$xvol)){
    stop("pmob$xvol is missing")
  } else if (any(is.na(pmob$xvol))){
    stop("missing values in pmob$xvol")
  }

  if(is.null(pmob$yvol)){
    stop("pmob$yvol is missing")
  } else if (any(is.na(pmob$yvol))){
    stop("missing values in pmob$yvol")
  }

  if(is.null(pmob$zvol)){
    stop("pmob$zvol is missing")
  } else if (any(is.na(pmob$zvol))){
    stop("missing values in pmob$zvol")
  }

  defaults <- list(DISCRETE_ID = "UNSET",
                   DEPTH = "null",
                   RUN_NUMBER = "-1",
                   TIMESTAMP = "UNSET",
                   SLOT_NUMBER = "-1",
                   MEAS_TYPE = "UNSET",
                   X_MOMENT = "0",
                   Y_MOMENT = "0",
                   Z_MOMENT = "0",
                   MAG_SUS = "NaN",
                   VOLUME = "10.8",
                   AREA = "4.0",
                   SAMPLE_AZ = "NaN",
                   SAMPLE_DIP = "NaN",
                   FORM_AZ = "NaN",
                   FORM_DIP = "NaN",
                   MAG_DEV = "0",
                   TREATMENT = "UNKNOWN",
                   AF_X = "NaN",
                   AF_Y = "NaN",
                   AF_Z = "NaN",
                   TEMPERATURE = "NaN",
                   IRM_FIELD = "NaN",
                   ARM_FIELD = "NaN",
                   ARM_AXIS = "UNKNOWN",
                   PP_SELECTED = "false",
                   PP_ANCHOR_PCA = "false",
                   PP_HIDDEN = "false",
                   PP_ONCIRCLE = "false",
                   PP_INPCA = "false")

  populated_column <- function(x, rep = nrow(pmob))
  {

    if(!is.null(x)){

      out <- rep(T, length(x))

      out[is.na(x)] <- F

      out[x == 0] <- F

    } else {

      out <- rep(F, rep)

    }

    # print(out)
    return(out)

  }

  pmob.cols <- c(names(formals(as.pmob))[-1], "version")

  match.pmob <- match(names(pmob), pmob.cols)

  regs <- pmob[,!is.na(match.pmob),drop = F]

  others <- pmob[,is.na(match.pmob), drop = F]

  # Convert regs ----

  output <- data.frame(matrix(ncol = 0, nrow = nrow(pmob)))

  output$DISCRETE_ID <- pmob$sampid

  output$DEPTH <- sprintf("%.3f", pmob$depth)

  if(!is.null(pmob$measid)) {

    if(!inherits(pmob$pmob$measid, "integer") |
       !inherits(pmob$pmob$measid, "numeric")){

      stop("To convert measid into RUN_NUMBER, it should be convertible ",
           "into an integer. Otherwise, suppress measid.")

    }

    output$RUN_NUMBER <- as.integer(pmob$measid + 0.1)

  }

  to.translate <- strptime(paste0(pmob$measyear, "-",
                                  pmob$measmonth, "-",
                                  pmob$measday, " ",
                                  pmob$meashour, ":",
                                  pmob$measmin, ":",
                                  "00"),
                           "%Y-%m-%d %H:%M:%S")

  output$TIMESTAMP <- format(to.translate, time.format)

  if(!is.null(pmob$slotid)) {

    if(!inherits(pmob$pmob$slotid, "integer") |
       !inherits(pmob$pmob$slotid, "numeric")){

      stop("To convert slotid into SLOT_NUMBER, it should be convertible ",
           "into an integer. Otherwise, suppress slotid.")

    }

    output$SLOT_NUMBER <- as.integer(pmob$slotsid + 0.1)

  }

  if(!is.null(pmob$discrete)){

    output$MEAS_TYPE <- rep("CONTINUOUS", nrow(pmob))

    output$MEAS_TYPE[pmob$discrete] <- "DISCRETE"

  } else {

    output$MEAS_TYPE <- rep("DISCRETE", nrow(pmob))

  }

  # ----

  x.moment <- sprintf("%.15E", pmob$xvol)

  x.rep.plus <- grepl("E+0", x.moment)
  x.rep.minus <- grepl("E-0", x.moment)

  x.moment[x.rep.plus] <- sub("E+0", "E+", x.moment[x.rep.plus])
  x.moment[x.rep.minus] <- sub("E-0", "E-", x.moment[x.rep.minus])

  output$X_MOMENT <- x.moment

  # ----

  y.moment <- sprintf("%.15E", pmob$yvol)

  y.rep.plus <- grepl("E+0", y.moment)
  y.rep.minus <- grepl("E-0", y.moment)

  y.moment[y.rep.plus] <- sub("E+0", "E+", y.moment[y.rep.plus])
  y.moment[y.rep.minus] <- sub("E-0", "E-", y.moment[y.rep.minus])

  output$Y_MOMENT <- y.moment

  # ----

  z.moment <- sprintf("%.15E", pmob$zvol)

  z.rep.plus <- grepl("E+0", z.moment)
  z.rep.minus <- grepl("E-0", z.moment)

  z.moment[z.rep.plus] <- sub("E+0", "E+", z.moment[z.rep.plus])
  z.moment[z.rep.minus] <- sub("E-0", "E-", z.moment[z.rep.minus])

  output$Z_MOMENT <- z.moment

  # ----

  if(magsus.type[1] == "vol") {

    if(is.null(pmob$volmagsus)) {
      stop("The volume normalised magnetic susceptibility is absent in the pmob")
    }

    ms <- pmob$volmagsus

  } else if(magsus.type[1] == "mass"){

    if(is.null(pmob$massmagsus)) {
      stop("The mass normalised magnetic susceptibility is absent in the pmob")
    }

    ms <- pmob$massmagsus

  } else if(magsus.type[1] == "tot"){

    if(is.null(pmob$totmagsus)) {
      stop("The total magnetic susceptibility is absent in the pmob")
    }

    ms <- pmob$totmagsus

  } else {

    ms <- pmob$MAG_SUS

  }

  if(!is.null(ms)){

    magsus <- sprintf("%.15E", ms)

    ms.rep.plus <- grepl("E+0", magsus)
    ms.rep.minus <- grepl("E-0", magsus)

    magsus[ms.rep.plus] <- sub("E+0", "E+", magsus[ms.rep.plus])
    magsus[ms.rep.minus] <- sub("E-0", "E-", magsus[ms.rep.minus])

    output$MAG_SUS <- magsus

  }

  if(!is.null(pmob$vol))  output$VOLUME <- sprintf("%.3f", pmob$vol * 10^6)
  if(!is.null(pmob$area)) output$AREA   <- sprintf("%.3f", pmob$area * 10^4)

  if(!is.null(pmob$samprot)){
    if(any(pmob$samprot != 0)) stop("samprot should be 0")
  }

  if(!is.null(pmob$sampaz))  output$SAMPLE_AZ  <- sprintf("%.0f", pmob$sampaz)
  if(!is.null(pmob$sampdip)) output$SAMPLE_DIP <- sprintf("%.0f", pmob$sampdip)

  if(!is.null(pmob$coraz))  if(any(pmob$coraz != 0)) stop("coraz should be 0")
  if(!is.null(pmob$cordip)) if(any(pmob$cordip != 90)) stop("cordip should be 90")
  if(!is.null(pmob$corrot)) if(any(pmob$corrot != 0)) stop("corrot should be 0")

  if(!is.null(pmob$bedaz)) output$FORM_AZ  <- sprintf("%.0f", pmob$bedaz)
  if(!is.null(pmob$bedaz)) output$FORM_DIP <- sprintf("%.0f", pmob$beddip)

  if(any(pmob$usemagaz & pmob$magazvarcor)){
    stop("The magnetic azimuth needs to be corrected by the magnetic variation")
  }

  if(any(pmob$bedazvarcor)){
    stop("The bedding azimuth needs to be corrected by the magnetic variation")
  }

  output$MAG_DEV <- defaults$MAG_DEV

  if(!is.null(pmob$TREATMENT)) {

    ppl.treats <- c("NONE", "DEGAUSS_XYZ", "DEGAUSS_Z",
                    "ARM", "IRM",  "THERMAL", "UNKNOWN", NA)

    test.treat <- pmob$TREATMENT %in% ppl.treats

    if(any(!test.treat)) {
      stop("If provided, the TREATMENT parameter should be made of the",
           " following entries: ", paste(ppl.treats, collapse = ", "))
    }

    output$TREATMENT <- pmob$TREATMENT

  } else {

    t.temp <- !is.null(pmob$treattemp)
    t.irm  <- !is.null(pmob$treatirmx)
    t.arm  <- !is.null(pmob$treatarmafx)

    t.afxyz <- !is.null(pmob$treatafx) &
      !is.null(pmob$treatafy)&
      !is.null(pmob$treatafz)

    t.afz <- !is.null(pmob$treatafz) &
      is.null(pmob$treatafx) &
      is.null(pmob$treatafy)

    if(t.afxyz){

      tt.afxyz <- T
      tt.afz   <- F

    } else if(t.afz){

      tt.afxyz <- F
      tt.afz   <- T

    } else {

      tt.afxyz <- F
      tt.afz   <- F

    }

    sum.treat <- as.integer(t.temp) +
      as.integer(tt.afxyz) +
      as.integer(tt.afz) +
      as.integer(t.irm) +
      as.integer(t.arm)

    if(sum.treat == 1){

      if(t.temp) {

        output$TREATMENT <- "THERMAL"

      } else if(tt.afxyz){

        output$TREATMENT <- "DEGAUSS_XYZ"

      } else if(tt.afz){

        output$TREATMENT <- "DEGAUSS_Z"

      } else if(t.irm){

        output$TREATMENT <- "IRM"

      } else if(t.arm){

        output$TREATMENT <- "ARM"

      }

    } else {

      output$TREATMENT <- "UNKNOWN"

    }

  }

  if(!is.null(pmob$treattemp)) {
    output$TEMPERATURE <- sprintf("%.0f", pmob$treattemp - 273.15)
  }

  if((!is.null(pmob$treatafx) |
      !is.null(pmob$treatafy) |
      !is.null(pmob$treatafz)) &
     (!is.null(pmob$treatarmafx) |
      !is.null(pmob$treatarmafy) |
      !is.null(pmob$treatarmafz) |
      !is.null(pmob$treatarmbiasx) |
      !is.null(pmob$treatarmbiasy) |
      !is.null(pmob$treatarmbiasz)
     )){

    warning("The .ppl format combines AF and ARM treatment into indentical",
            " columns; they cannot be accomodated together (and quite frankly,",
            " they should not be, they are 2 completely different procedures)")

  }

  if(!is.null(pmob$treatarmafx) |
     !is.null(pmob$treatarmafy) |
     !is.null(pmob$treatarmafz) |
     !is.null(pmob$treatarmbiasx) |
     !is.null(pmob$treatarmbiasy) |
     !is.null(pmob$treatarmbiasz)){

    if(!is.null(pmob$treatarmafx)) output$AF_X <- sprintf("%.0f",pmob$treatarmafx)
    if(!is.null(pmob$treatarmafy)) output$AF_Y <- sprintf("%.0f",pmob$treatarmafy)
    if(!is.null(pmob$treatarmafz)) output$AF_Z <- sprintf("%.0f",pmob$treatarmafz)

    nz <- any(populated_column(pmob$treatarmbiasx) |
                populated_column(pmob$treatarmbiasy))

    if(nz){
      stop(".ppl files only support ARM bias field applied along the z axis")
    }


    if(!is.null(pmob$pmob$treatarmbiasz)) {
      output$ARM_FIELD <- sprintf("%.0f", pmob$treatarmbiasz)
    }

    output$ARM_AXIS <- "AXIAL"


  } else if(!is.null(pmob$treatafx) |
            !is.null(pmob$treatafy) |
            !is.null(pmob$treatafz)){

    if(!is.null(pmob$treatafx)) output$AF_X <- sprintf("%.0f", pmob$treatafx)
    if(!is.null(pmob$treatafy)) output$AF_Y <- sprintf("%.0f", pmob$treatafy)
    if(!is.null(pmob$treatafz)) output$AF_Z <- sprintf("%.0f", pmob$treatafz)

  }

  output$PP_SELECTED <- pmob$PP_SELECTED
  output$PP_HIDDEN   <- pmob$PP_HIDDEN

  output$PP_ANCHOR_PCA <- "false"
  output$PP_ONCIRCLE   <- "false"
  output$PP_INPCA      <- "false"

  if(!is.null(pmob$pcaanchor) ){
    output$PP_ANCHOR_PCA[pmob$pcaanchor] <- "true"
  }

  if(!is.null(pmob$pcacomponent) & !is.null(pmob$pcacomponentsingle)){

    output$PP_ONCIRCLE[pmob$pcacomponent == pmob$pcacomponentsingle] <- "true"

  }

  if(!is.null(pmob$circlecomponent) & !is.null(pmob$circlecomponentsingle)){

    circle.pos <- pmob$circlecomponent == pmob$circlecomponentsingle

    output$PP_ONCIRCLE[circle.pos] <- "true"

  }

  output$CREATION_DATE      <- pmob$CREATION_DATE

  if(is.null(pmob$MODIFICATION_DATE)){

    current.time <- strftime(as.POSIXlt(Sys.time()), "%Y-%m-%dT%H:%M:%S%z")

    output$MODIFICATION_DATE <- sub('(+[0-9]{2})([0-9]{2}$)','\\1:\\2',
                                    current.time , fixed=FALSE)

  } else {

    output$MODIFICATION_DATE <- pmob$MODIFICATION_DATE

  }

  output$ORIGINAL_FILE_TYPE <- pmob$ORIGINAL_FILE_TYPE

  if(is.null(pmob$ORIGINAL_CREATOR_PROGRAM)){
    output$ORIGINAL_CREATOR_PROGRAM <- "inteRmag 0.1.0 pmob format 0.0.0.9020"
  } else {
    output$ORIGINAL_CREATOR_PROGRAM <- pmob$ORIGINAL_CREATOR_PROGRAM
  }

  if(is.null(pmob$SAVED_BY_PROGRAM)){
    output$SAVED_BY_PROGRAM <- "inteRmag 0.1.0 pmob format 0.0.0.9020"
  } else {
    output$SAVED_BY_PROGRAM <- pmob$SAVED_BY_PROGRAM
  }

  output$HEIGHT <- pmob$HEIGHT

  output$latitude  <- pmob$lat
  output$longitude <- pmob$long

  output$SITE        <- pmob$SITE
  output$declination <- pmob$declination
  output$inclination <- pmob$inclination

  ppl.cols <- c("DISCRETE_ID",
                "DEPTH",
                "RUN_NUMBER",
                "TIMESTAMP",
                "SLOT_NUMBER",
                "MEAS_TYPE",
                "X_MOMENT",
                "Y_MOMENT",
                "Z_MOMENT",
                "MAG_SUS",
                "VOLUME",
                "AREA",
                "SAMPLE_AZ",
                "SAMPLE_DIP",
                "FORM_AZ",
                "FORM_DIP",
                "MAG_DEV",
                "TREATMENT",
                "AF_X",
                "AF_Y",
                "AF_Z",
                "TEMPERATURE",
                "IRM_FIELD",
                "ARM_FIELD",
                "ARM_AXIS",
                "PP_SELECTED",
                "PP_ANCHOR_PCA",
                "PP_HIDDEN",
                "PP_ONCIRCLE",
                "PP_INPCA",
                "MEASUREMENT_TYPE",
                "CREATION_DATE",
                "MODIFICATION_DATE",
                "ORIGINAL_FILE_TYPE",
                "ORIGINAL_CREATOR_PROGRAM",
                "SAVED_BY_PROGRAM",
                "HEIGHT",
                "latitude",
                "longitude","SITE",
                "declination",
                "inclination")

  if(isTRUE(other.columns)){

    additional <- others[,is.na(match(colnames(others), ppl.cols)), drop = F]

    output <- cbind(output, additional)

  } else if(is.character(other.columns)){

    if(any(other.columns %in% c(pmob.cols, ppl.cols))){

      stop("The columns selected by others.columns should not be regular .ppl ",
           "columns or regular pmob columns, e.g., \nnames(formals(as.pmob))\n")

    }

    additional <- pmob[,match(other.columns, colnames(pmob)), drop = F]

    output <- cbind(output, additional)

  }

  # Default values ----

  redefault.pos <- match(colnames(output), names(defaults))

  redefault <- defaults[redefault.pos, drop = F]

  names(redefault) <- colnames(output)

  as.is <- sapply(redefault, is.null)

  redefault[as.is]  <- as.list(output)[as.is]
  redefault[!as.is] <- lapply(redefault[!as.is],
                              function(x) rep(x, nrow(output)))

  default.vector <- as.vector(as.matrix(data.frame(redefault)))
  output.vector  <- as.vector(as.matrix(output))

  output.vector[is.na(output.vector)] <- default.vector[is.na(output.vector)]

  output2 <- data.frame(matrix(output.vector, ncol = ncol(output)))

  colnames(output2) <- colnames(output)

  # Write the .ppl file ----

  write.ppl(output2, file)

}


