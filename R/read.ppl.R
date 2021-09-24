#' @title Read .ppl files
#'
#' @description Convert a .ppl file (from the PuffinPlot format) into a table
#'
#' @param file the directory of the file to convert
#'
#' @references Lurcock, P. C., and G. S., Wilson. 2012. "PuffinPlot: A versatile,
#' user-friendly program for paleomagnetic analysis". Geochemistry, Geophysics,
#' Geosystems, 13, Q06Z45, https://doi.org/10.1029/2012GC004098.
#'
#' @examples
#' file <- system.file("tests",
#'                     "PuffinPlot_example.ppl",
#'                     package = "inteRmag")
#'
#' read.ppl(file)
#'
#' @export
#' @importFrom readr read_delim
#' @importFrom dplyr full_join
#' @importFrom stringr str_count

read.ppl <- function(file)
{
  x <- readLines(file)

  lf <- length(x)

  if(!grepl("PuffinPlot file\\. Version 3", x[1])){
    stop("The first line in the file should be:\n",
         "\tPuffinPlot file. Version 3")
  }

  sep <- which(nchar(x) == 0)

  if(length(sep) != 1) stop(paste("There should be one empty line",
                                  " to delimit measurements and site",
                                  sep = ""))

  meas <- x[2:(sep - 1)]

  ncol1 <- str_count(meas[1], "\t") + 1

  m <- as.data.frame(read_delim(I(c(meas, "")), delim = "\t",
                                col_types = strrep("c", ncol1)))

  m <- m[-nrow(m),]

  # ----

  smt <- which(grepl("^SUITE\tMEASUREMENT_TYPE\t", x))

  if(length(smt) != 0){
    mt <- sub("^SUITE\tMEASUREMENT_TYPE\t","" ,x[smt])
    if(!(mt == "DISCRETE" | mt == "continuous")) {
      stop("Invalid Suite Measurement Type: should be DISCRETE or CONTINUOUS")
    }
  } else {
    mt <- NULL
  }

  scf <- which(grepl("^SUITE\tCUSTOM_FLAG_NAMES", x))
  scn <- which(grepl("^SUITE\tCUSTOM_NOTE_NAMES", x))

  if(length(scf) != 0){
    cf <- strsplit(sub("^SUITE\tCUSTOM_FLAG_NAMES\t","" ,x[scf]), split = "\t")
    cf <- unlist(cf)
  } else {
    cf <- NULL
  }

  if(length(scn) != 0){
    cn <- strsplit(sub("^SUITE\tCUSTOM_NOTE_NAMES\t","" ,x[scn]), split = "\t")
    cn <- unlist(cn)
  } else {
    cn <- NULL
  }

  scd <- which(grepl("^SUITE\tCREATION_DATE", x))
  smd <- which(grepl("^SUITE\tMODIFICATION_DATE", x))

  if(length(scn) != 0){
    cd <- sub("^SUITE\tCREATION_DATE\t","" ,x[scd])
  } else {
    cd <- NULL
  }

  if(length(smd) != 0){
    md <- sub("^SUITE\tMODIFICATION_DATE\t","" ,x[smd])
  } else {
    md <- NULL
  }

  sof <- which(grepl("^SUITE\tORIGINAL_FILE_TYPE", x))
  soc <- which(grepl("^SUITE\tORIGINAL_CREATOR_PROGRAM", x))
  ssp <- which(grepl("^SUITE\tSAVED_BY_PROGRAM", x))

  if(length(sof) != 0){

    acc <- c("TWOGEE", "ZPLOT", "PUFFINPLOT_OLD", "PUFFINPLOT_NEW", "CALTECH",
             "IAPD", "UCDAVIS", "DIRECTIONS", "CUSTOM_TABULAR", "PMD_ENKIN",
             "JR6" ,"UNKNOWN")

    of <- sub("^SUITE\tORIGINAL_FILE_TYPE\t","" ,x[sof])

    if(!(of %in% acc)){
      stop("The Original File Type should be one of the following:\n",
           paste(acc, collapse = "\n"))
    }

  } else {
    of <- NULL
  }

  if(length(soc) != 0){
    oc <- sub("^SUITE\tORIGINAL_CREATOR_PROGRAM\t","" ,x[soc])
  } else {
    oc <- NULL
  }

  if(length(ssp) != 0){
    sp <- sub("^SUITE\tSAVED_BY_PROGRAM\t","" ,x[ssp])
  } else {
    sp <- NULL
  }


  # ----

  cus.flag <- which(grepl("^SAMPLE\t.+\tCUSTOM_FLAGS", x))
  cus.note <- which(grepl("^SAMPLE\t.+\tCUSTOM_NOTES", x))

  if(length(cus.flag) != 0){

    ncol2 <- str_count(x[cus.flag][1], "\t") + 1

    flag <- as.data.frame(read_delim(I(c(x[cus.flag], "")), delim = "\t",
                                     col_names = F,
                                     col_types = strrep("c", ncol2)))[,-c(1,3)]

    # flag <- flag[-nrow(flag),]

    if(ncol(flag) - 1L != length(cf)){
      stop("There should as many Custom Flag Names defined in Suite",
           " than found in Sample")
    }

    names(flag) <- c("DISCRETE_ID", cf)

  } else {

    flag <- NULL

  }

  if(length(cus.note) != 0){

    ncol3 <- str_count(x[cus.note][1], "\t") + 1

    note <- as.data.frame(read_delim(I(c(x[cus.note], "")), delim = "\t",
                                     col_names = F,
                                     col_types = strrep("c", ncol3)))[,-c(1,3)]

    # note <- note[-nrow(note),]

    if(ncol(note) - 1L != length(cn)){
      stop("There should as many Custom Note Names defined in Suite",
           " than found in Sample")
    }

    names(note) <- c("DISCRETE_ID", cn)

  } else {

    note <- NULL

  }

  # ----

  site.id <- which(grepl("^SAMPLE\t.+\tSITE", x))

  if(length(site.id) != 0){

    ncol4 <- str_count(x[site.id][1], "\t") + 1

    id.site <- as.data.frame(read_delim(I(c(x[site.id], "")), delim = "\t",
                                        col_names = F,
                                        col_types = strrep("c",
                                                           ncol4)))[,-c(1,3)]

    # id.site <- id.site[-nrow(id.site),]

    names(id.site) <- c("DISCRETE_ID", "SITE")

  } else {

    id.site <- NULL

  }

  imp.dir <- which(grepl("^SAMPLE\t.+\tIMPORTED_DIRECTION", x))

  if(length(imp.dir) != 0){

    ncol5 <- str_count(x[imp.dir][1], "\t") + 1

    dir.imp <- as.data.frame(read_delim(I(c(x[imp.dir], "")), delim = "\t",
                                        col_names = F,
                                        col_types = strrep("c",
                                                           ncol5)))[,-c(1,3)]

    # dir.imp <- dir.imp[-nrow(dir.imp),]

    names(dir.imp) <- c("DISCRETE_ID", "declination", "inclination")

  } else {

    dir.imp <- NULL

  }

  site.h <- which(grepl("^SITE\t.+\tHEIGHT", x))

  if(length(site.h) != 0){

    ncol6 <- str_count(x[site.h][1], "\t") + 1

    h.site <- as.data.frame(read_delim(I(c(x[site.h], "")), delim = "\t",
                                       col_names = F,
                                       col_types = strrep("c",
                                                          ncol6)))[,-c(1,3)]

    # h.site <- h.site[-nrow(h.site),]

    names(h.site) <- c("SITE", "height")

  } else {

    h.site <- NULL

  }

  site.l <- which(grepl("^SITE\t.+\tLOCATION", x))

  if(length(site.l) != 0){

    ncol7 <- str_count(x[site.l][1], "\t") + 1

    l.site <- as.data.frame(read_delim(I(c(x[site.l], "")), delim = "\t",
                                       col_names = F,
                                       col_types = strrep("c",
                                                          ncol7)))[,-c(1,3)]

    # l.site <- l.site[-nrow(l.site),]

    names(l.site) <- c("SITE", "latitude", "longitude")

  } else {

    l.site <- NULL

  }

  # ----

  all <- c(seq_len(sep),
           smt, scf, scn, scd, smd, sof, soc, ssp,
           cus.flag, cus.note,
           site.id, imp.dir, site.h, site.l)

  check <- which(!(seq_len(length(x)) %in% all))

  if(length(check) != 0){
    warning("The following lines in the .ppl files are not recognised:\n",
            paste(check, collapse = " "))
  }

  accu <- m

  if(!is.null(note)) accu <- full_join(accu, note, by = "DISCRETE_ID")
  if(!is.null(flag)) accu <- full_join(accu, flag, by = "DISCRETE_ID")

  if(!is.null(id.site)) accu <- full_join(accu, id.site, by = "DISCRETE_ID")
  if(!is.null(dir.imp)) accu <- full_join(accu, dir.imp, by = "DISCRETE_ID")

  if(!is.null(h.site)) accu <- full_join(accu, h.site, by = "SITE")
  if(!is.null(l.site)) accu <- full_join(accu, l.site, by = "SITE")

  sing <- list(MEASUREMENT_TYPE = mt,
               CREATION_DATE = cd, MODIFICATION_DATE = md,
               ORIGINAL_FILE_TYPE = of, ORIGINAL_CREATOR_PROGRAM = oc,
               SAVED_BY_PROGRAM = sp)

  sl <- nrow(accu)

  sli <- lapply(sing, rep, each = sl)

  info <- sapply(sing, is.null)

  more <- data.frame(matrix(unlist(sli), nrow=sl),stringsAsFactors=FALSE)

  names(more) <- names(sing)[!info]

  accu <- cbind(accu, more)

  return(accu)

}

