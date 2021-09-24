#' @title Attribute certain pmob values as NA
#'
#' @description Attribute certain values in a pmob
#' (\strong{P}aleo\strong{M}agnetic \strong{OB}ject) as NA (Not Available), and
#' allow the columns entirely made of NA values to be removed
#'
#' @param pmob pmob object
#' @param l list of parameters and the values that are defined as NA
#' @param remove whether to remove the columns made only of NA values (TRUE), or
#' to keep them (FALSE; is the default)
#'
#' @examples
#' pmob1 <- as.pmob(sampid = c("holder", "Bosso1", "Bosso2", "Bosso3"),
#'                  specid = c("holder", "b", "c", "NA"),
#'                  slotid = c("holder", "1", "2", "3"),
#'                  treattemp = c(20, 20, 20, 20),
#'                  long = rep(0, 4), lat = rep(0,4))
#'
#' pmob1
#'
#' pmob.na(pmob1,
#'         l = list(specid = c("holder", "NA"),
#'                  slotid = "holder",
#'                  long = 0, lat = 0),
#'         remove = TRUE)
#'
#' @export
#' @importFrom StratigrapheR merge_list

pmob.na <- function(pmob, l = list(), remove = FALSE)
{

  if(!is.pmob(pmob)) stop("The 'pmob' parameter should be a pmob")

  if(!inherits(l, "list")) stop("The 'l' parameter should be a list")

  check.l <- l[names(l) %in% names(pmob)]

  check.p <- as.list(pmob[,match(names(check.l), names(pmob)), drop = F])

  na.loc <- mapply(FUN = function(x,y) x %in% y,
                   x = check.p,
                   y = check.l,
                   SIMPLIFY = F)

  out <- as.data.frame(check.p, stringsAsFactors = F)

  pos.na <- which(as.data.frame(na.loc) == TRUE,
                  arr.ind=TRUE)

  out[pos.na] <- NA

  output <- data.frame(merge_list(as.list(out), as.list(pmob)),
                       stringsAsFactors = F)

  output <- pmob.col.sort(output)

  if(isTRUE(remove)){

    all.na <- apply(output, 2, function(x) all(is.na(x)))

    output <- output[, !all.na, drop = FALSE]

  }

  return(output)

}





