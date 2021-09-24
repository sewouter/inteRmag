#' @title Order the columns of a pmob
#'
#' @description Order the columns of a pmob (\strong{P}aleo\strong{M}agnetic
#' \strong{OB}ject).
#'
#' @param pmob Object for which the columns need to be ordered (same order as
#' the arguments in \code{\link{as.pmob}}.
#'
#' @examples
#' pmob <- as.pmob(sampid = "Bosso1", exc = 1, measmin = NA)
#' pmob <- pmob[,c(1,5,3,2,4)]
#'
#' pmob
#'
#' pmob <- pmob.col.sort(pmob)
#'
#' pmob
#'
#' @export

pmob.col.sort <- function(pmob){

  if(!is.pmob(pmob)) stop("The 'pmob' parameter should be a pmob")

  target <- c(names(formals(as.pmob))[-1], "version")

  pos <- match(names(pmob), target)

  res <- pmob[,order(pos), drop = FALSE]

  return(res)
}
