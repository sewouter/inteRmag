#' @title Check and normalize measurements values by blank and holder
#'  measurements
#'
#' @description Check and normalize sample (and holder) measurements values by
#' blank and holder measurements. This function comes with features allowing to
#' detect abnormalities such as excessive time between sample measurements and
#' blank/holder measurements, excessive blank holder/values, etc. Problematic
#' data lines in the pmob can be outputted if needed.
#'
#' @param pmob a pmob object, for which the measurements hae not been corrected
#' (for blank and/or holder)
#' @param id which column to use as a measurement id (by default, "sampid")
#' @param correct.blank whether to correct for blank measurement (TRUE by
#' default)
#' @param correct.holder whether to correct for the holder (TRUE by default);
#' the measurements should already be corrected for blanks.
#' @param compress.blank whether to remove all the blank measurements from the
#' pmob after correction. This is not advised, to avoid information loss (FALSE
#' by default).
#' @param compress.holder whether to remove all the holder measurements from the
#' pmob after correction. This is not advised, to avoid information loss (FALSE
#' by default).
#' @param method.blank The methodology used for computing the blank value at
#' the time of the sample measurement for correction. Supported values are
#' "approx" for a linear approximation between the previous and the next blank
#' values, using the time of the measurement, "mean" stands for a mean of the
#' previous and the next blank values, "previous" takes the previous blank
#' value, and "next", the next one.
#' @param method.holder The methodology used for computing the holder value at
#' the time of the sample measurement for correction. Supported values are
#' "approx" for a linear approximation between the previous and the next holder
#' values, using the time of the measurement, "mean" stands for a mean of the
#' previous and the next holder values, "previous" takes the previous holder
#' value, and "next", the next one.
#' @param plot whether to plot the behaviors of the blank measurements (x, y
#' and z) and of the holder measurements (x, y and z).
#' @param time.format the format of the plot's time axis; the default,
#' "day 1" means month-day (the month in abbreviated textual form), "day2"
#' year-month-day (the month in abbreviated textual form), and "hour", "min",
#' "sec" the respective time units starting at zero with the first measurement.
#' @param center.blank Whether to center the xint, yint and zint of the blank
#' measurement plots on zero.
#' @param center.holder Whether to center the xint, yint and zint of the holder
#' measurement plots on zero.
#' @param contact.blank Whether to check if the blank measurements are
#' directly contiguous with the sample or holder measurements they
#' correct.
#' @param contact.holder Whether to check if the holder measurements are
#' directly contiguous with the sample measurements they correct.
#' @param id.contact.blank Whether to check that the blank measurements values
#' used to correct the sample measurements are directly contiguous with the
#' holder/sample measurements \bold{of identical id} they correct.
#' @param id.contact.holder Whether to check that the holder measurements values
#' used to correct the sample measurements are directly contiguous with the
#' sample measurements \bold{of identical id} they correct.
#' @param tlim.blank The maximum time interval allowed for measurements implied
#' in the same sample blank correction
#' @param tlim.holder The maximum time interval allowed for measurements implied
#' in the same sample holder correction
#' @param blank.x.lim The range of x blank values allowed. This is found in the
#' plots as red lines. If one value is provided, lim = c(-lim, lim), otherwise,
#' lim = range(lim).
#' @param blank.y.lim The range of y blank values allowed. This is found in the
#' plots as red lines. If one value is provided, lim = c(-lim, lim), otherwise,
#' lim = range(lim).
#' @param blank.z.lim The range of z blank values allowed. This is found in the
#' plots as red lines. If one value is provided, lim = c(-lim, lim), otherwise,
#' lim = range(lim).
#' @param holder.x.lim The range of x holder values allowed. This is found in
#' the plots as red lines. If one value is provided, lim = c(-lim, lim),
#' otherwise, lim = range(lim).
#' @param holder.y.lim The range of y holder values allowed. This is found in
#' the plots as red lines. If one value is provided, lim = c(-lim, lim),
#' otherwise, lim = range(lim).
#' @param holder.z.lim The range of z holder values allowed. This is found in
#' the plots as red lines. If one value is provided, lim = c(-lim, lim),
#' otherwise, lim = range(lim).
#' @param output.suspect.lines Whether to have the suspect lines as an output,
#' rather than the resulting pmob.
#'
#' @seealso The \code{\link{background.correction}} function
#'
#' @examples
#' frontone14.pmob
#'
#' pmob.correct.measure(frontone14.pmob, plot = F)
#'
#' @importFrom dplyr lead lag
#' @export

pmob.correct.measure <- function(pmob, id = "sampid",
                                 correct.blank  = TRUE,
                                 correct.holder  = TRUE,
                                 compress.blank = FALSE,
                                 compress.holder = FALSE,
                                 method.blank  = c("mean", "approx",
                                                   "previous", "next"),
                                 method.holder = c("previous", "next",
                                                   "mean", "approx"),
                                 plot = TRUE,
                                 time.format = c("day1", "day2",
                                                 "hour", "min", "sec"),
                                 center.blank = FALSE,
                                 center.holder = FALSE,
                                 contact.blank = FALSE,
                                 contact.holder = FALSE,
                                 id.contact.blank = FALSE,
                                 id.contact.holder = FALSE,
                                 tlim.blank = NULL,
                                 tlim.holder = NULL,
                                 blank.x.lim = NULL,
                                 blank.y.lim = NULL,
                                 blank.z.lim = NULL,
                                 holder.x.lim = NULL,
                                 holder.y.lim = NULL,
                                 holder.z.lim = NULL,
                                 output.suspect.lines = FALSE)
{

  if(!is.pmob(pmob)) stop("Incorrect pmob object")

  # Provide id ----

  if(length(id) == 0) stop("The 'id' parameter should not be of length 0.")

  id.pos <- match(id, colnames(pmob))

  if(any(is.na(id.pos))) stop("Unidentified 'id' column in the pmob object.")

  id.code <- apply(pmob[,id.pos,drop = F], 1, paste, collapse = "/")

  # Check blank and holder parameters ----

  if(!isTRUE(correct.blank) & !isTRUE(correct.holder)){

    stop("If correct.blank & correct.holder are both FALSE,",
         " this function is useless.")

  }

  if(isTRUE(correct.blank)){

    if(is.null(pmob$isblank)) stop("Missing pmob$isblank column")
    if(is.null(pmob$blankcor)) stop("Missing pmob$blankcor column")

    if(!is.null(pmob$blankmethod)) {
      warning("pmob$blankmethod will be overwritten by the function")
    }

  }

  if(isTRUE(correct.holder)){

    if(is.null(pmob$isholder)) stop("Missing pmob$isholder column")
    if(is.null(pmob$holdercor)) stop("Missing pmob$holdercor column")

    if(!is.null(pmob$holdermethod)) {
      warning("pmob$holdermethod will be overwritten by the function")
    }

  }

  # Check for blankcor and holdercor,
  # and the relation to correct.blank and correct.holder ----

  if(isTRUE(correct.blank)){

    if(any(pmob$blankcor[!pmob$isblank])){

      stop('Some measurements were already corrected for blank background')

    }

  }

  if(isTRUE(correct.holder)){

    if(any(pmob$holdercor[!pmob$isholder & !pmob$isblank])){

      stop('Some measurements were already corrected for the holder')

    }

  }


  if(!isTRUE(correct.blank) & isTRUE(correct.holder)){

    if(any(pmob$blankcor[!pmob$isblank])) {

      stop("If no blank correction is required, all the non-blank ",
           "measurements should be blank corrected")

    }

  }

  if(any(pmob$holdercor)){

    stop("No measurement should already be corrected for the holder")

  }

  # Check for complete date information,
  # and define index in function of that ----

  null.time <- is.null(pmob$measyear) | is.null(pmob$measmonth) |
    is.null(pmob$measday) | is.null(pmob$meashour) |
    is.null(pmob$measmin) | is.null(pmob$meassec)

  if(!null.time){

    na.time <- is.na(pmob$measyear) | is.na(pmob$measmonth) |
      is.na(pmob$measday) | is.na(pmob$meashour) |
      is.na(pmob$measmin) | is.na(pmob$meassec)

  }

  if(null.time) {

    if(method.blank[1] == "approx" | method.holder[1] == "approx"){

      stop("Missing measurement time information in the pmob",
           " (meassec, measmin, meashour, measday, measmonth or measyear)")

    } else {

      time <- seq(nrow(pmob))

      index <- T

    }

  } else {

    if(any(na.time)){

      if(method.blank[1] == "approx" | method.holder[1] == "approx"){

        na.time.pos <- which(na.time)

        stop("NA values in the measurement time information of the pmob",
             " (meassec, measmin, meashour, measday, measmonth or measyear), ",
             "at lines ",
             paste(na.time.pos, collapse = ", "))

      } else {

        time <- seq(nrow(pmob))

        index <- T

      }

    } else {

      format.time <- paste0(pmob$measyear, "-",
                            pmob$measmonth, "-",
                            pmob$measday, " ",
                            pmob$meashour, ":",
                            pmob$measmin, ":",
                            round(pmob$meassec,0))

      date.time <- as.POSIXlt(format.time, format = c("%Y-%m-%d %H:%M:%S"))

      cont.time <- as.numeric(date.time) + pmob$meassec - trunc(pmob$meassec)

      time <- cont.time - cont.time[1]


      if(time.format[1] == "day1" | time.format[1] == "day2"){

        zero.time <- paste0(pmob$measyear[1], "-",
                            pmob$measmonth[1], "-",
                            pmob$measday[1], " ",
                            0, ":", 0, ":", 0)

        date.zero.time <- as.POSIXlt(zero.time, format = c("%Y-%m-%d %H:%M:%S"))

        num.zero.time <- as.numeric(date.zero.time)

        axis.sec <- seq(num.zero.time, max(cont.time) + 60*60*24, 60*60*24)

        if(time.format[1] == "day1") {
          dateFormat <- "%b %d"
        } else if(time.format[1] == "day2") {
          dateFormat <- "%Y %b %d"
        }

        osys <- Sys.getlocale(category = "LC_TIME")

        NOTUSED <- Sys.setlocale("LC_TIME", "C")

        axis.date <- format(as.Date(as.POSIXlt(axis.sec,
                                               origin = "1970-01-01")),
                            dateFormat)

        NOTUSED <- Sys.setlocale("LC_TIME", osys)

      }

      if(is.unsorted(cont.time)){
        stop("Measurements should be sorted in chronological ",
             "order of measurements")
      }

      if(time.format[1] == "hour") {
        time.factor <- 60*60
      } else if(time.format[1] == "min"){
        time.factor <- 60
      } else if(time.format[1] == "sec" |
                time.format[1] == "day1" |
                time.format[1] == "day2" ){
        time.factor <- 1
      }

      index <- F

    }

  }

  # Plot blank behaviour ----

  if(isTRUE(correct.blank)){

    if(length(blank.x.lim) == 1) {
      blank.x.lim <- c(-blank.x.lim, blank.x.lim)
    } else if(!is.null(blank.x.lim)) {
      blank.x.lim <- range(blank.x.lim)
    }

    if(length(blank.y.lim) == 1) {
      blank.y.lim <- c(-blank.y.lim, blank.y.lim)
    } else if(!is.null(blank.y.lim)) {
      blank.y.lim <- range(blank.y.lim)
    }

    if(length(blank.z.lim) == 1) {
      blank.z.lim <- c(-blank.z.lim, blank.z.lim)
    } else if(!is.null(blank.z.lim)) {
      blank.z.lim <- range(blank.z.lim)
    }

    if(isTRUE(plot)){

      opar <- par("mfrow")

      par(mfrow = c(3,1))

      xylim <- range(c(pmob$xint[pmob$isblank], blank.x.lim))
      yylim <- range(c(pmob$yint[pmob$isblank], blank.y.lim))
      zylim <- range(c(pmob$zint[pmob$isblank], blank.z.lim))

      if(isTRUE(center.blank)){

        xylim <- c(-max(abs(xylim)), max(abs(xylim)))
        yylim <- c(-max(abs(yylim)), max(abs(yylim)))
        zylim <- c(-max(abs(zylim)), max(abs(zylim)))

      }

      if(index){

        plot(time[pmob$isblank], pmob$xint[pmob$isblank],
             ylim = xylim, pch = 19, type = "l",
             ylab = "Blank X intensity (A per squared m)",
             xlab = "Index")
        abline(h = blank.x.lim, col = "red")
        if(isTRUE(center.blank)) abline(h = 0)

        plot(time[pmob$isblank], pmob$yint[pmob$isblank],
             ylim = yylim, pch = 19, type = "l",
             ylab = "Blank Y intensity (A per squared m)",
             xlab = "Index")
        abline(h = blank.y.lim, col = "red")
        if(isTRUE(center.blank)) abline(h = 0)

        plot(time[pmob$isblank], pmob$zint[pmob$isblank],
             ylim = zylim, pch = 19, type = "l",
             ylab = "Blank Z intensity (A per squared m)",
             xlab = "Index")
        abline(h = blank.z.lim, col = "red")
        if(isTRUE(center.blank)) abline(h = 0)


      } else {

        if(time.format[1] == "sec" |
           time.format[1] == "min" |
           time.format[1] == "hour"){

          plot(time[pmob$isblank]/time.factor, pmob$xint[pmob$isblank],
               ylim = xylim, pch = 19, type = "l",
               ylab = "Blank X intensity (A per squared m)",
               xlab = paste0("Time (", time.format[1], ")"))
          abline(h = blank.x.lim, col = "red")
          if(isTRUE(center.blank)) abline(h = 0)

          plot(time[pmob$isblank]/time.factor, pmob$yint[pmob$isblank],
               ylim = yylim, pch = 19, type = "l",
               ylab = "Blank Y intensity (A per squared m)",
               xlab = paste0("Time (", time.format[1], ")"))
          abline(h = blank.y.lim, col = "red")
          if(isTRUE(center.blank))abline(h = 0)

          plot(time[pmob$isblank]/time.factor, pmob$zint[pmob$isblank],
               ylim = zylim, pch = 19, type = "l",
               ylab = "Blank Z intensity (A per squared m)",
               xlab = paste0("Time (", time.format[1], ")"))
          abline(h = blank.z.lim, col = "red")
          if(isTRUE(center.blank)) abline(h = 0)

        } else if (time.format[1] == "day1" | time.format[1] == "day2" ){

          plot(time[pmob$isblank]/time.factor, pmob$xint[pmob$isblank],
               ylim = xylim, pch = 19, type = "l",
               ylab = "Blank X intensity (A per squared m)",
               xlab = "Date (start of the day)", axes = F)

          axis(2)
          axis(1, at = axis.sec - cont.time[1], labels = axis.date)
          box()

          abline(h = blank.x.lim, col = "red")
          if(isTRUE(center.blank)) abline(h = 0)

          plot(time[pmob$isblank]/time.factor, pmob$yint[pmob$isblank],
               ylim = yylim, pch = 19, type = "l",
               ylab = "Blank Y intensity (A per squared m)",
               xlab = "Date (start of the day)", axes = F)

          axis(2)
          axis(1, at = axis.sec - cont.time[1], labels = axis.date)
          box()

          abline(h = blank.y.lim, col = "red")
          if(isTRUE(center.blank))  abline(h = 0)

          plot(time[pmob$isblank]/time.factor, pmob$zint[pmob$isblank],
               ylim = zylim, pch = 19, type = "l",
               ylab = "Blank Z intensity (A per squared m)",
               xlab = "Date (start of the day)", axes = F)

          axis(2)
          axis(1, at = axis.sec - cont.time[1], labels = axis.date)
          box()

          abline(h = blank.z.lim, col = "red")
          if(isTRUE(center.blank)) abline(h = 0)

        }

      }


      par(mfrow = opar)

    }

  }

  # blank background correction ----

  if(isTRUE(correct.blank)){

    xblank <- background.correction(x = pmob$xint,
                                    bg = pmob$isblank,
                                    t = time,
                                    id = id.code,
                                    method = method.blank[1],
                                    contact = contact.blank,
                                    id.contact = contact.blank,
                                    warn = F,
                                    tlim = tlim.blank,
                                    bglim = blank.x.lim)

    yblank <- background.correction(x = pmob$yint,
                                    bg = pmob$isblank,
                                    t = time,
                                    id = id.code,
                                    method = method.blank[1],
                                    contact = F,
                                    id.contact = F,
                                    warn = F,
                                    tlim = NULL,
                                    bglim = blank.y.lim)

    zblank <- background.correction(x = pmob$zint,
                                    bg = pmob$isblank,
                                    t = time,
                                    id = id.code,
                                    method = method.blank[1],
                                    contact = F,
                                    id.contact = F,
                                    warn = F,
                                    tlim = NULL,
                                    bglim = blank.z.lim)

    nline <- xblank$line

    pmob2 <- pmob[nline, ]
    time2 <- time[nline]

    pmob2$xint <- xblank$corrected.x
    pmob2$yint <- yblank$corrected.x
    pmob2$zint <- zblank$corrected.x

    pmob2$blankcor <- T

  } else {

    nline <- seq(nrow(pmob))[!pmob$isblank,]

    pmob2 <- pmob[!pmob$isblank,]

  }

  # Plot holder ----

  if(isTRUE(correct.holder)){

    if(length(holder.x.lim) == 1) {
      holder.x.lim <- c(-holder.x.lim, holder.x.lim)
    } else if(!is.null(holder.x.lim)) {
      holder.x.lim <- range(holder.x.lim)
    }

    if(length(holder.y.lim) == 1) {
      holder.y.lim <- c(-holder.y.lim, holder.y.lim)
    } else if(!is.null(holder.y.lim)) {
      holder.y.lim <- range(holder.y.lim)
    }

    if(length(holder.z.lim) == 1) {
      holder.z.lim <- c(-holder.z.lim, holder.z.lim)
    } else if(!is.null(holder.z.lim)) {
      holder.z.lim <- range(holder.z.lim)
    }

    if(isTRUE(plot)){

      opar <- par("mfrow")

      par(mfrow = c(3,1))

      xylim2 <- range(c(pmob2$xint[pmob2$isholder], holder.x.lim))
      yylim2 <- range(c(pmob2$yint[pmob2$isholder], holder.y.lim))
      zylim2 <- range(c(pmob2$zint[pmob2$isholder], holder.z.lim))

      if(isTRUE(center.holder)){

        xylim2 <- c(-max(abs(xylim2)), max(abs(xylim2)))
        yylim2 <- c(-max(abs(yylim2)), max(abs(yylim2)))
        zylim2 <- c(-max(abs(zylim2)), max(abs(zylim2)))

      }

      if(index){

        plot(time2[pmob2$isholder], pmob2$xint[pmob2$isholder],
             ylim = xylim2, pch = 19, type = "l",
             ylab = "Holder X intensity (A per squared m)",
             xlab = "Index")
        abline(h = holder.x.lim, col = "red")
        if(isTRUE(center.holder)) abline(h = 0)

        plot(time2[pmob2$isholder], pmob2$yint[pmob2$isholder],
             ylim = yylim2, pch = 19, type = "l",
             ylab = "Holder Y intensity (A per squared m)",
             xlab = "Index")
        abline(h = holder.y.lim, col = "red")
        if(isTRUE(center.holder)) abline(h = 0)

        plot(time2[pmob2$isholder], pmob2$zint[pmob2$isholder],
             ylim = zylim2, pch = 19, type = "l",
             ylab = "Holder Z intensity (A per squared m)",
             xlab = "Index")
        abline(h = holder.z.lim, col = "red")
        if(isTRUE(center.holder)) abline(h = 0)

      } else {

        if(time.format[1] == "sec" |
           time.format[1] == "min" |
           time.format[1] == "hour"){

          plot(time2[pmob2$isholder]/time.factor, pmob2$xint[pmob2$isholder],
               ylim = xylim2, pch = 19, type = "l",
               ylab = "Holder X intensity (A per squared m)",
               xlab = paste0("Time (", time.format[1], ")"))
          abline(h = holder.x.lim, col = "red")
          if(isTRUE(center.holder)) abline(h = 0)

          plot(time2[pmob2$isholder]/time.factor, pmob2$yint[pmob2$isholder],
               ylim = yylim2, pch = 19, type = "l",
               ylab = "Holder Y intensity (A per squared m)",
               xlab = paste0("Time (", time.format[1], ")"))
          abline(h = holder.y.lim, col = "red")
          if(isTRUE(center.holder)) abline(h = 0)

          plot(time2[pmob2$isholder]/time.factor, pmob2$zint[pmob2$isholder],
               ylim = zylim2, pch = 19, type = "l",
               ylab = "Holder Z intensity (A per squared m)",
               xlab = paste0("Time (", time.format[1], ")"))
          abline(h = holder.z.lim, col = "red")
          if(isTRUE(center.holder)) abline(h = 0)

        } else if (time.format[1] == "day1" | time.format[1] == "day2" ){

          plot(time2[pmob2$isholder]/time.factor, pmob2$xint[pmob2$isholder],
               ylim = xylim2, pch = 19, type = "l",
               ylab = "Holder X intensity (A per squared m)",
               xlab = "Date (start of the day)", axes = F)

          axis(2)
          axis(1, at = axis.sec - cont.time[1], labels = axis.date)
          box()

          abline(h = holder.x.lim, col = "red")
          if(isTRUE(center.holder)) abline(h = 0)

          plot(time2[pmob2$isholder]/time.factor, pmob2$yint[pmob2$isholder],
               ylim = yylim2, pch = 19, type = "l",
               ylab = "Holder Y intensity (A per squared m)",
               xlab = "Date (start of the day)", axes = F)

          axis(2)
          axis(1, at = axis.sec - cont.time[1], labels = axis.date)
          box()

          abline(h = holder.y.lim, col = "red")
          if(isTRUE(center.holder)) abline(h = 0)

          plot(time2[pmob2$isholder]/time.factor, pmob2$zint[pmob2$isholder],
               ylim = zylim2, pch = 19, type = "l",
               ylab = "Holder Z intensity (A per squared m)",
               xlab = "Date (start of the day)", axes = F)

          axis(2)
          axis(1, at = axis.sec - cont.time[1], labels = axis.date)
          box()

          abline(h = holder.z.lim, col = "red")
          if(isTRUE(center.holder)) abline(h = 0)

        }

      }

      par(mfrow = opar)

    }

  }

  out.pmob <- pmob

  if(isTRUE(correct.holder)){

    xholder <- background.correction(x = pmob2$xint,
                                     bg = pmob2$isholder,
                                     t = time2,
                                     id = id.code[nline],
                                     method = method.holder[1],
                                     contact = contact.holder,
                                     id.contact = id.contact.holder,
                                     warn = F,
                                     tlim = tlim.holder,
                                     bglim = holder.x.lim)

    yholder <- background.correction(x = pmob2$yint,
                                     bg = pmob2$isholder,
                                     t = time2,
                                     id = id.code[nline],
                                     method = method.holder[1],
                                     contact = F,
                                     id.contact = F,
                                     warn = F,
                                     tlim = NULL,
                                     bglim = holder.y.lim)

    zholder <- background.correction(x = pmob2$zint,
                                     bg = pmob2$isholder,
                                     t = time2,
                                     id = id.code[nline],
                                     method = method.holder[1],
                                     contact = F,
                                     id.contact = F,
                                     warn = F,
                                     tlim = NULL,
                                     bglim = holder.z.lim)

    meas.lines   <- xblank$line[xholder$line]
    holder.lines <- xblank$line[-xholder$line]

    out.pmob$xint[meas.lines] <- xholder$corrected.x
    out.pmob$yint[meas.lines] <- yholder$corrected.x
    out.pmob$zint[meas.lines] <- zholder$corrected.x

    out.pmob$holdercor[meas.lines] <- T

    if(isTRUE(correct.blank)){

      out.pmob$xint[holder.lines] <- pmob2$xint[pmob2$isholder]
      out.pmob$yint[holder.lines] <- pmob2$yint[pmob2$isholder]
      out.pmob$zint[holder.lines] <- pmob2$zint[pmob2$isholder]

      out.pmob$blankcor[meas.lines]   <- T
      out.pmob$blankcor[holder.lines] <- T

      out.pmob <- pmob.add(pmob = out.pmob,
                           blankmethod = method.blank[1],
                           holdermethod = method.holder[1])

    } else {

      out.pmob <- pmob.add(pmob = out.pmob,
                           holdermethod = method.holder[1])

    }



  } else {

    out.pmob$xint[!pmob$isblank] <- pmob2$xint
    out.pmob$yint[!pmob$isblank] <- pmob2$yint
    out.pmob$zint[!pmob$isblank] <- pmob2$zint

    out.pmob$blankcor[!pmob$isblank] <- pmob2$blankcor

    out.pmob <- pmob.add(pmob = out.pmob,
                         blankmethod = method.blank[1])

  }

  if(isTRUE(compress.blank) & isTRUE(compress.holder)){

    out.pmob <- out.pmob[!out.pmob$isblank & !out.pmob$isholder,]

  } else if(isTRUE(compress.blank) & !isTRUE(compress.holder)){

    out.pmob <- out.pmob[!out.pmob$isblank,]

  } else if(isTRUE(compress.holder) & !isTRUE(compress.blank)) {

    stop("If 'compress.holder' is TRUE, 'compress.blank' should also be TRUE.")

  }

  if(isTRUE(correct.blank) & isTRUE(correct.holder)){

    suspect.lines <- list(blank.incorrect.contact =
                            xblank$suspect.lines$incorrect.contact,
                          blank.incorrect.id.contact =
                            xblank$suspect.lines$incorrect.id.contact,
                          blank.excessive.time =
                            xblank$suspect.lines$excessive.time,
                          blank.excessive.x.background =
                            xblank$suspect.lines$excessive.background,
                          blank.excessive.y.background =
                            yblank$suspect.lines$excessive.background,
                          blank.excessive.z.background =
                            zblank$suspect.lines$excessive.background,
                          holder.incorrect.contact =
                            nline[xholder$suspect.lines$incorrect.contact],
                          holder.incorrect.id.contact =
                            nline[xholder$suspect.lines$incorrect.id.contact],
                          holder.excessive.time =
                            nline[xholder$suspect.lines$excessive.time],
                          holder.excessive.x.background =
                            nline[xholder$suspect.lines$excessive.background],
                          holder.excessive.y.background =
                            nline[yholder$suspect.lines$excessive.background],
                          holder.excessive.z.background =
                            nline[zholder$suspect.lines$excessive.background])

  } else if(isTRUE(correct.blank) & !isTRUE(correct.holder)){

    suspect.lines <- list(blank.incorrect.contact = NULL,
                          blank.incorrect.id.contact = NULL,
                          blank.excessive.time = NULL,
                          blank.excessive.x.background = NULL,
                          blank.excessive.y.background = NULL,
                          blank.excessive.z.background = NULL,
                          holder.incorrect.contact =
                            nline[xholder$suspect.lines$incorrect.contact],
                          holder.incorrect.id.contact =
                            nline[xholder$suspect.lines$incorrect.id.contact],
                          holder.excessive.time =
                            nline[xholder$suspect.lines$excessive.time],
                          holder.excessive.x.background =
                            nline[xholder$suspect.lines$excessive.background],
                          holder.excessive.y.background =
                            nline[yholder$suspect.lines$excessive.background],
                          holder.excessive.z.background =
                            nline[zholder$suspect.lines$excessive.background])

  } else if(!isTRUE(correct.blank) & isTRUE(correct.holder)){

    suspect.lines <- list(blank.incorrect.contact =
                            xblank$suspect.lines$incorrect.contact,
                          blank.incorrect.id.contact =
                            xblank$suspect.lines$incorrect.id.contact,
                          blank.excessive.time =
                            xblank$suspect.lines$excessive.time,
                          blank.excessive.x.background =
                            xblank$suspect.lines$excessive.background,
                          blank.excessive.y.background =
                            yblank$suspect.lines$excessive.background,
                          blank.excessive.z.background =
                            zblank$suspect.lines$excessive.background,
                          holder.incorrect.contact = NULL,
                          holder.incorrect.id.contact = NULL,
                          holder.excessive.time = NULL,
                          holder.excessive.x.background = NULL,
                          holder.excessive.y.background = NULL,
                          holder.excessive.z.background = NULL)

  }

  is.suspect <- !sapply(suspect.lines, function(x) length(x) == 0)

  if(any(is.suspect)){

    warning("Problematic data lines detected")
    warning("")

    w1 <- names(suspect.lines)
    w2 <- unlist(lapply(suspect.lines, paste, collapse = ", "), use.names = F)

    w3 <- paste(w1, w2, sep = " at lines: \n")

    warning(paste(w3, collapse = "\n"))

    warning("")

    warning("Set the output.suspect.lines parameter to TRUE for more info.")

  }

  if(isTRUE(output.suspect.lines)){

    return(suspect.lines)

  } else {

    return(out.pmob)

  }

}





