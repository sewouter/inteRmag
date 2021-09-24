#' @title Merge two pmobs together
#'
#' @description Merge two pmobs (\strong{P}aleo\strong{M}agnetic
#' \strong{OB}jects) together based on the columns taken account as standard
#' (see parameters in \code{\link{as.pmob}}) and found in both pmobs. This
#' function acts like a left join: the first (left) pmob will be completed by
#' the second (right) pmob.
#'
#' @param pmob1 pmob object to complement
#' @param pmob2 pmob object complementing pmob1
#' @param byinitialorder whether to use the \code{initialorder} parameter to
#' merge the two pmobs. Is set as FALSE as default, as initialorder is
#' automatically generated for all pmob objects, and should be considered
#' carefully
#'
#' @examples
#' pmob1 <- as.pmob(sampid = c("Bosso1","Bosso1", "Bosso2", "Bosso3"),
#'                  specid = rep("b", 4),
#'                  slotid = c(1,1,2,3),
#'                  treattemp = c(20,150, 20, 20))
#'
#' pmob2 <- as.pmob(sampid = c("Bosso1", "Bosso2", "Bosso3", "Bosso4"),
#'                  specid = rep("b", 4),
#'                  slotid = c(1,1,3,4),
#'                  sampip = c(30, 50, 45, 20))
#'
#' pmob.merge(pmob1 = pmob1, pmob2 = pmob2)
#'
#' @export
#' @importFrom dplyr left_join

pmob.merge <- function(pmob1, pmob2, byinitialorder = FALSE){

  if(!is.pmob(pmob1)) stop("The 'pmob1' parameter should be a pmob")
  if(!is.pmob(pmob2)) stop("The 'pmob2' parameter should be a pmob")

  nr1 <- nrow(pmob1)
  nr2 <- nrow(pmob2)

  if(nr1 < nr2) {
    stop("'pmob1' should have more (or the same",
         " amount of) elements than 'pmob2'")
  }

  pmob1 <- pmob.col.sort(pmob1)
  pmob2 <- pmob.col.sort(pmob2)

  mergeid <- intersect(names(pmob1), names(pmob2))

  target <- c(names(formals(as.pmob))[-1], "version")

  if(any(!(mergeid %in% target))){

    remove <- mergeid[!(mergeid %in% target)]

    warning("The following columns, found in both pmobs,",
            " will not be used for merging:\n - ",
            paste(remove, collapse = "\n - "),
            "\n Use the merge() function if you want to use them for merging")

    mergeid <- mergeid[-which(mergeid %in% remove)]

  }

  if(isTRUE(byinitialorder)){
    mergeid <- mergeid[-which(mergeid %in% c("version"))]
  } else {
    mergeid <- mergeid[-which(mergeid %in% c("initialorder", "version"))]
  }

  pmob1.mi <- pmob1[,names(pmob1) %in% mergeid, drop = F]
  pmob2.mi <- pmob2[,names(pmob2) %in% mergeid, drop = F]

  pmob1.un <- do.call(paste, c(pmob1.mi, sep="-"))
  pmob2.un <- do.call(paste, c(pmob2.mi, sep="-"))

  if(length(unique(pmob2.un)) != nr2){

    dupl <- which(duplicated(pmob2.un) | duplicated(pmob2.un, fromLast = T))

    stop("The following lines in pmob2 are duplicated: \n    ",
         paste(dupl, collapse = ", "),
         "\ni.e. all of these lines have at least one duplicate having the ",
         "same values \nin the columns used for merging the two pmobs, ",
         "which are the following:\n - ",
         paste(mergeid, collapse = "\n - "))

  }

  not.in <- !(pmob1.un %in% pmob2.un)

  if(any(not.in)){

    w.not.in <- which(not.in)

    warning("The following line(s) in pmob1 is (are) not found in pmob2: \n    ",
            paste(w.not.in, collapse = ", "),
            "\ni.e. all of these lines have no equivalent in pmob2 in the",
            " columns\nused for merging the two pmobs, which are the ",
            "following:\n - ",
            paste(mergeid, collapse = "\n - "),
            "\nThis will generate NA (Not Available) values.")

  }

  if(isTRUE(byinitialorder)){
    add.pmob2 <- pmob2[,-which(names(pmob2) %in% c("version")),
                       drop = F]
  } else {
    add.pmob2 <- pmob2[,-which(names(pmob2) %in% c("version", "initialorder")),
                       drop = F]
  }

  res <- left_join(pmob1, add.pmob2, by = mergeid,
                   suffixes = c(".pmob1", ".pmob2"))

  res <- pmob.col.sort(res)

  return(res)

}


