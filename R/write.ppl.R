#' @title Write .ppl files
#'
#' @description Convert a table into a .ppl file (from the PuffinPlot format)
#'
#' @param x the table to convert.
#' @param file the directory of the file to convert
#'
#' @references Lurcock, P. C., and G. S., Wilson. 2012. "PuffinPlot: A versatile,
#' user-friendly program for paleomagnetic analysis". Geochemistry, Geophysics,
#' Geosystems, 13, Q06Z45, https://doi.org/10.1029/2012GC004098.
#'
#' @examples
#' \dontrun{
#' write.ppl(data.frame(), "empty_file")}
#'
#' @export
#' @importFrom StratigrapheR seq_mult

write.ppl <- function(x, file = "")
{

  filename <- paste0(file, ".ppl")

  # Expected headers ----

  titles <- c("DISCRETE_ID",
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
              "PP_INPCA")

  suite.titles <- c("MEASUREMENT_TYPE",
                    "CREATION_DATE",
                    "MODIFICATION_DATE",
                    "ORIGINAL_FILE_TYPE",
                    "ORIGINAL_CREATOR_PROGRAM",
                    "SAVED_BY_PROGRAM")

  site.titles <- c("HEIGHT",
                   "latitude",
                   "longitude")

  sample.titles <- c("SITE",
                     "declination",
                     "inclination")

  # ----

  dupl <- duplicated(x$DISCRETE_ID)

  # Prepare for default: lack of measurement type

  if(nrow(x) == 0) x <- data.frame(MEASUREMENT_TYPE = "DISCRETE")

  if(is.null(x$MEASUREMENT_TYPE)) x$MEASUREMENT_TYPE <- "DISCRETE"

  if(is.null(x$MEAS_TYPE)) x$MEAS_TYPE <- "DISCRETE"

  # Part 1 : first line----

  part1 <- "PuffinPlot file. Version 3"

  # Part 2 : main data ----

  p2.pos <- match(colnames(x), titles)

  repos.p2 <- p2.pos[!is.na(p2.pos)]

  part2 <- x[,!is.na(p2.pos), drop = F]

  part2 <- part2[,order(repos.p2), drop = F]

  part2.header <- paste(colnames(part2), collapse = "\t")

  part2.text <- apply(part2, 1, function(x) paste(x, collapse = "\t"))


  # Part 4 : Sites --------------------------------------------------------------------------------

  if(any(colnames(x) == "SITE")){

    p4.pos <- match(colnames(x), site.titles)

    repos.p4 <- p4.pos[!is.na(p4.pos)]

    part4 <- x[,!is.na(p4.pos), drop = F]

    part4 <- part4[,order(repos.p4), drop = F]

    part4.collapsed <- apply(part4, 1, function(x) paste(x, collapse = " "))

    to.test.part4 <- split(part4.collapsed, x$SITE)

    test.part4 <- sapply(to.test.part4, function(x) length(unique(x)))

    if(any(test.part4 != 1)){

      prob.names <- paste(names(test.part4)[which(test.part4 != 1)],
                          collapse = ", ")

      stop("Height, latitude or longitude show variability ",
           "inside sites ", prob.names)

    }

    dupl.sites <- duplicated(x$SITE)

    p4.comp <- part4[!dupl.sites, ,drop = F]

    sites.comp <- x$SITE[!dupl.sites]

    if(!is.null(x$HEIGHT)){

      df2.height <- data.frame(site = "SITE", si = sites.comp,
                               height = "HEIGHT", he = p4.comp$HEIGHT)

      df2.height <- df2.height[!is.na(df2.height$si),]

      coll.height <- apply(df2.height, 1, function(x) paste(x, collapse = "\t"))

    } else {

      coll.height <- NULL

    }

    if(!is.null(x$latitude) & !is.null(x$longitude)){

      df2.loc <- data.frame(site = "SITE", si = sites.comp,
                            location = "LOCATION",
                            lat = p4.comp$latitude,
                            long = p4.comp$longitude)


      df2.loc <- df2.loc[!is.na(df2.loc$si),]


      coll.loc <- apply(df2.loc, 1, function(x) paste(x, collapse = "\t"))

    } else {

      coll.loc <- NULL

    }

    if(!is.null(coll.height) & !is.null(coll.loc)){

      part4.text <- c(coll.height, coll.loc)[seq_mult(length(coll.height)*2,
                                                      2, inv = T)]

    } else {

      part4.text <- c(coll.height, coll.loc)

    }

  } else {

    part4.text <- NULL

  }

  # Part 5: Flags, Notes, Sites and Imported Directions ----

  # Sample to site connection ----

  if(!is.null(x$SITE)){

    site.df <- data.frame(sample = x$DISCRETE_ID, site = x$SITE)

    site.test <- split(site.df$site, site.df$sample)

    test.site <- sapply(site.test, function(x) length(unique(x)))

    if(any(test.site != 1)){

      prob.names.site <- paste(names(test.site)[which(test.site != 1)],
                               collapse = ", ")

      stop("Different sites provided for the samples ", prob.names.sites)

    }

    site.comp <- site.df[!dupl,]

    site.df2 <- data.frame(samp = "SAMPLE", site.comp$sample,
                           si = "SITE", site = site.comp$site)

    partconn <- apply(site.df2, 1, function(x) paste(x, collapse = "\t"))

    na.conn <- is.na(site.df2$site)

  } else {

    partconn <- NULL

    na.conn <- NULL
    # na.site <- NULL

  }

  # declination and inclination ----

  if(any(colnames(x) == "declination") &
     any(colnames(x) == "inclination")){

    import.df <- data.frame(declination = x$declination,
                            inclination = x$inclination)

    import.collapsed <- apply(import.df, 1, function(x) paste(x, collapse = " "))

    to.test.import <- split(import.collapsed, x$DISCRETE_ID)

    test.import <- sapply(to.test.import, function(x) length(unique(x)))

    if(any(test.import != 1)){

      prob.names <- paste(names(test.import)[which(test.import != 1)],
                          collapse = ", ")

      stop("Imported direction shows variability ",
           "inside samples ", prob.names)

    }

    df2.import <- cbind(data.frame(samp = "SAMPLE",
                                   sample = x$DISCRETE_ID[!dupl],
                                   import = "IMPORTED_DIRECTION"),
                        import.df[!dupl,])

    import.text <- apply(df2.import, 1, function(x) paste(x, collapse = "\t"))

    na.import <- is.na(df2.import$declination) & is.na(df2.import$inclination)


  } else {

    import.text <- NULL

    na.import <- NULL

  }

  # Flag and note ----


  known <- c(titles, suite.titles, site.titles, sample.titles)

  others <- x[,is.na(match(colnames(x), known)), drop = F]

  if(ncol(others) == 0){

    flags.text <- NULL
    notes.text <- NULL

    st2 <- NULL
    st3 <- NULL

    na.flag <- NULL
    na.note <- NULL


  } else {

    collapsed <- apply(others, 1, function(x) paste(x, collapse = " "))

    spl.collapsed <- split(collapsed, x$DISCRETE_ID)

    test.others <- sapply(spl.collapsed, function(x) length(unique(x)))

    if(any(test.others != 1)){

      prob.names <- paste(names(test.others)[which(test.others != 1)],
                          collapse = ", ")

      stop("Additional information to be added shows variability ",
           "inside samples ", prob.names)

    }

    comp.others <- others[!dupl,,drop = F]

    test.others.2 <- apply(comp.others, 2,
                           function(x) all(x == "true" | x == "false"))

    flags.pos <- test.others.2 == T

    df.flags <- comp.others[,flags.pos, drop = F]
    df.notes <- comp.others[,!flags.pos, drop = F]

    # Empty values for notes ----

    inter.note <- as.vector(as.matrix(df.notes))

    inter.note[is.na(inter.note)] <- ""

    df.notes2 <- as.data.frame(matrix(inter.note, ncol = ncol(df.notes)))

    colnames(df.notes2) <- colnames(df.notes)

    # ----

    if(ncol(df.notes2) != 0){

      df2.notes <- cbind(data.frame(SAMPLE = "SAMPLE",
                                    samp = unique(x$DISCRETE_ID),
                                    FLAGS = "CUSTOM_NOTES"), df.notes2)

      notes.text <- apply(df2.notes, 1, function(x) paste(x, collapse = "\t"))

      s3 <- colnames(df.notes2)

      st3 <- paste(c("SUITE", "CUSTOM_NOTE_NAMES", s3), collapse = "\t")

      na.note <- rep(F, length(notes.text))

    } else {

      st3        <- NULL
      notes.text <- NULL
      na.note    <- NULL

    }

    if(ncol(df.flags) != 0){

      df2.flags <- cbind(data.frame(SAMPLE = "SAMPLE",
                                    samp = unique(x$DISCRETE_ID),
                                    FLAGS = "CUSTOM_FLAGS"), df.flags)

      flags.text <- apply(df2.flags, 1, function(x) paste(x, collapse = "\t"))

      s2 <- colnames(df.flags)

      st2 <- paste(c("SUITE", "CUSTOM_FLAG_NAMES", s2), collapse = "\t")

      na.flag <- rep(F, length(flags.text))

    } else {

      st2        <- NULL
      flags.text <- NULL
      na.flag    <- NULL

    }

  }

  # Remix ----

  part5.dis <- c(flags.text, notes.text, partconn, import.text)

  if(!is.null(part5.dis)){
    mult <- as.integer(!is.null(flags.text)) +
      as.integer(!is.null(notes.text)) +
      as.integer(!is.null(partconn)) +
      as.integer(!is.null(import.text))

    recombine <- seq_mult(part5.dis, mult, inv = T)

    na.all <- c(na.flag, na.note, na.conn, na.import)[recombine]

    part5 <- part5.dis[recombine]

    part5 <- part5[!na.all]

  } else {

    part5 <- NULL

  }



  # Part 6 : suite ----

  suite <- x[,!is.na(match(colnames(x), suite.titles)), drop = F]

  test.suite <- apply(suite, 2, function(x) length(unique(x)))

  if(any(test.suite != 1)){

    stop("Suite information should be identical for all samples")

  }


  comp.suite <- suite[1,,drop =F]

  s1 <- comp.suite[names(comp.suite) == "MEASUREMENT_TYPE"]

  st1 <- paste(c("SUITE", "MEASUREMENT_TYPE", s1), collapse = "\t")

  # ----

  s4.pos <- match(names(comp.suite), suite.titles[-1])

  s4.keep <- !is.na(s4.pos)

  s4 <- comp.suite[s4.keep]

  if(ncol(s4) != 0){

    s4 <- s4[order(s4.pos[s4.keep])]

    t4 <- cbind(data.frame(suite = "SUITE", names = colnames(s4)), t(s4))

    st4 <- apply(t4, 1, function(x) paste(x, collapse = "\t"))

  } else {

    st4 <- NULL

  }

  # write the file -----

  fileConn <- file(filename)

  writeLines(c(part1,
               part2.header,
               part2.text,
               "",
               part5,
               part4.text,
               st1,
               st2,
               st3,
               st4),
             fileConn)

  close(fileConn)

}

