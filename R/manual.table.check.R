#' @title Check a table for rare occurrences of values in its columns.
#'
#' @description If you expect a few mistakes in a table, this function will find
#' the amount of unique elements in each column, and if the amount is beyond the
#' 'n_elements' threshold, it will show the frequency of occurrence of the
#' elements. This helps identify rare values, which could be mistakes. Then, if
#' the frequencies are below a certain threshold (the 'frequency' parameter),
#' it will track down and show the individual suspect (potentially problematic)
#' lines. A bar plot is also integrated to better visualise how values (in
#' alphabetic order) are distributed.
#'
#' @param x the data frame (made entirely of character class columns) to check.
#' @param n_elements the threshold of the amount of unique elements below which
#' further investigation is required to check for rare occurrences. This can
#' also be a ratio of the total number of rows in the table (by default 0.01,
#' i.e. 1 percent).
#' @param frequency the frequency of occurrence of unique elements below which
#' individual lines in the table are to be shown. This can also be a ratio of
#' the total number of rows in the table (by default 0.01, i.e. 1 percent).
#' @param plot whether to show a barplot of each occurrence of a value in each
#' column. This plot helps to visualise anomalies.
#' @param indicator_columns a character vector of the headers of the columns
#' to be shown in the suspect lines. The default (TRUE) shows all the columns.
#'
#' @examples
#' set.seed(42)
#'
#' # Evenly distributed values ----
#'
#' binary <- sample(0:1, 1000, replace = T)
#'
#' letters <- sample(c("a", "b", "c"), 1000, replace = T)
#'
#' random  <- round(runif(1000), 2)
#'
#' # Anomalies ----
#'
#' binary[90:100] <- 2
#'
#' letters[6]      <- "GG"
#' letters[60:400] <- "EZ"
#'
#' random[400:450] <- 0.5
#'
#' # ----
#'
#' x <- data.frame(binary = binary, letters = letters, random = random)
#'
#' opar <- par("mfrow")
#'
#' par(mfrow = c(1,3))
#'
#' manua.table.check(x) # The plot shows anomalies that are more
#' # complex than detected only by the function
#'
#' par(mfrow = opar)
#'
#' manual.table.check(x, frequency = 100, plot = F)
#'
#' @export

manual.table.check <- function(x,
                               n_elements = 0.01, frequency = 0.01,
                               plot = TRUE, indicator_columns = TRUE)
{

  if(!inherits(x, "data.frame")) stop("'x' should be of class 'data.frame'")

  test.character <- apply(x, 2, function(x) inherits(x, "character"))

  if(any(!test.character)) {

    stop("All columns should be of class 'character'.")

  }

  nrows <- nrow(x)

  if(n_elements < 1) n_elements <- as.integer(nrows * n_elements)

  if(frequency < 1) frequency <- as.integer(nrows * frequency)

  unique_length <- apply(x, 2, function(x) length(unique(x)))

  out <- list(n_rows = nrows,
              n_unique_elements = unique_length,
              show_frequencies_if_n_unique_elements_less_than = n_elements)

  small_x <- x[, unique_length <= n_elements, drop = F]

  occurrences <- apply(small_x, 2, table)

  l_occurrences <- lapply(occurrences,
                          function(x) as.data.frame(x, stringsAsFactors = F))

  l_occurrences <- lapply(l_occurrences, function(x) {
    colnames(x) <- c("element", "frequency")
    return(x)
  })

  out$frequencies <- l_occurrences

  warning.list <- lapply(l_occurrences,
                         function(x) x[x$frequency <= frequency,])

  warning.list        <- warning.list[lapply(warning.list, nrow) != 0]

  accu.line <- c()

  accu.value <- c()

  accu.column <- c()

  for(i in seq(length(warning.list))){

    i.col <- names(warning.list)[i]

    all.values <- x[,match(i.col, colnames(x))]

    i.line <- which(all.values %in% warning.list[[i]]$element)

    accu.line <- c(accu.line, i.line)

    accu.column <- c(accu.column, rep(i.col, length(i.line)))

    accu.value <- c(accu.value, all.values[i.line])

  }

  suspect.df <- data.frame(column = accu.column,
                           value = accu.value,
                           line = accu.line,
                           fill = "")

  colnames(suspect.df) <- c("Problem Column", "Problem Value",
                            "line", "........")


  if(isTRUE(indicator_columns)){

    add <- x[suspect.df$line, ]

  } else {

    indicator_pos <- match(indicator_columns, colnames(x))

    if(any(is.na(indicator_pos))) {
      stop("indicator_columns not matched in data frame column header")
    }

    add <- x[suspect.df$line, indicator_pos]

  }

  suspect.df <- suspect.df[order(suspect.df$line),]

  suspect.df2 <- cbind(suspect.df, add)

  rownames(suspect.df2) <- NULL

  out$show_suspect_lines_if_frequency_less_than = frequency

  out$suspect.lines <- suspect.df2

  # Plot more difficult occurrences ----

  if(isTRUE(plot)){

    oparmar <- par("mar")

    on.exit(par(mar = oparmar))

    complete.table <- apply(x, 2, table)

    par(mar = c(10,4,4,2) + 0.1)

    for(i in  seq_len(ncol(x))){

      barplot(complete.table[[i]], las = 2, main = names(complete.table)[i])

    }

    par(mar = oparmar)

  }

  # print(out)

  return(out)

}


