#' @title Correct measurements values by background values taken before and/or
#' after
#'
#' @description Correct measurements values by background values taken before
#' and/or after, assuming a linear relationship with time. This function can warn
#' when background values exceed given values, or if the time elapsed between
#' measurement of the sample and of the background exceeds a given threshold.
#'
#' @param x The measurement and background values, in chronological order of
#' their measurement.
#' @param bg Whether the value is a background measurement (T), or a sample
#' measurement.
#' @param t A numerical value of time.
#' @param id Sample id.
#' @param method The method for computing the background value at the time of the
#' sample measurement; "approx" stands for al inear approximation between the
#' previous and the next background values, using the time of the measurement;
#' "mean" stands for a mean of the previous and the next background values;
#' "previous" takes the previous background value; "next", the next one.
#' @param compress Whether to remove the background measurements from the output.
#' @param contact Whether to check if the background measurements are directly
#' contiguous with the sample measurements they correct.
#' @param id.contact Whether to check that the background measruements used to
#' correct the sample measurements are directly contiguous with the sample
#' measurements of identical id they correct.
#' @param tlim The maximum time interval allowed for measurements implied in the
#' same sample correction.
#' @param bglim The range of background values allowed. If one value is provided,
#' bglim = c(-bglim, bglim), otherwise, bglim = range(bglim).
#' @param warn Whether to provide a warning when background values exceed given
#' values, or if the time elapsed between measurement of the sample and of the
#' background exceeds a given threshold (otherwise this information is only
#' provided in the output list).
#'
#' @return a list of the initial line of the measurement ($line), the sample id
#' ($id), the sample measurement time ($t), the original sample measurement value
#' ($original.x), the corrected sample measurement value ($corrected.x), the
#' computed background value used for the correction ($background.x), the time of
#' measurement of the previous background value ($previous.t), the
#' background value preceding the sample measurement ($previous.x), the time of
#' measurement of the next background value ($next.t), and the
#' background value following the sample measurement ($next.x). Additionally, a
#' list of the problematic lines of data is provided ($suspect.lines),
#' indicating the lines of incorrect contact ($suspect.lines$incorrect.contact),
#' of incorrect id contact ($suspect.lines$incorrect.id.contact), of excessive
#' background values ($suspect.lines$excessive.background) and of excessive time
#' (suspect.lines$excessive.time).
#'
#' @examples
#' x  <- c(2, 10, 3, 4, 25, 3)
#' bg <- c(T, F, T, T, F, T)
#' t  <- c(1:3, 20:22)
#' id <- c("S1", "S1", "S1", "S2", "S2", "S2")
#'
#' output <- background.correction(x = x, bg = bg, t = t, id = id)
#'
#' plot(t, x, pch = 19, ylim = c(0, 30))
#'
#' lines(t[bg], x[bg])
#'
#' points(output$t, output$background.x, col = "red", pch = 19)
#'
#' @importFrom dplyr lead lag
#' @export

background.correction <- function(x, bg, t = NULL, id = NULL,
                                  method = c("approx", "mean",
                                             "previous", "next"),
                                  compress = T,
                                  contact = F, id.contact = F,
                                  tlim = NULL, bglim = NULL,
                                  warn = T)
{

  if(!is.null(t)){

    if(is.unsorted(t)) stop("The values in 't' should be sorted.")

  }

  # Check parameters ----

  if(!(method[1] %in% c("approx", "mean", "previous", "next"))){

    stop("The first element of the 'method' parameter should be either",
         " 'approx', 'mean', 'previous' or 'next'.")

  }

  if(isTRUE(id.contact) & !isTRUE(contact)){

    stop("If the 'id.contact' parameter is TRUE, the 'contact' parameter",
         " should be TRUE too.")

  }

  if(!inherits(bg, "logical")){
    stop("The 'bg' parameter should be of class 'logical'.")
  }

  if(!is.null(t)){

    if(is.unsorted(t)){

      stop("The values in the 't' parameter should be sorted.")

    }

  }

  if(method[1] == "approx" & is.null(t)) {
    stop("If 'method' is 'approx', the 't' parameter should not be NULL")
  }

  if(is.null(id) & isTRUE(id.contact)){

    stop("If the 'id.contact' parameter is TRUE, 'id' should not be NULL")

  }

  # Error list ----

  suspect.lines <- list(incorrect.contact = NULL,
                 incorrect.id.contact = NULL,
                 excessive.background = NULL,
                 excessive.time = NULL)

  # Tests of contact ----

  if(isTRUE(contact)){

    if(!is.null(id)){

      test <- data.frame(id = id,            bg = bg,
                         lag.id = lag(id),   lag.bg = lag(bg),
                         lead.id = lead(id), lead.bg = lead(bg))

    } else {

      test <- data.frame(bg = bg,
                         lag.bg = lag(bg),
                         lead.bg = lead(bg))

    }

    if(method[1] == "approx" | method[1] == "mean"){

      test$testcontact <- test$lag.bg & test$lead.bg

      test$testcontact[test$bg] <- NA

      if(any(!test$testcontact, na.rm = T)){

        no.contact.pos <- which(!test$testcontact)

        suspect.lines$incorrect.contact <- no.contact.pos

        if(isTRUE(warn)){

          warning("Some measurements are not directly surrounded ",
                  "by blank measurements, at lines: ",
                  paste(no.contact.pos, collapse = ", "))

        }

      }

    } else if(method[1] == "previous"){

      test$testcontact <- test$lag.bg

      test$testcontact[test$bg] <- NA

      if(any(!test$testcontact, na.rm = T)){

        no.contact.pos <- which(!test$testcontact)

        suspect.lines$incorrect.contact <- no.contact.pos

        if(isTRUE(warn)){

          warning("Some measurements are not directly preceded ",
                  "by blank measurements, at lines: ",
                  paste(no.contact.pos, collapse = ", "))

        }

      }

    } else if(method[1] == "next"){

      test$testcontact <- test$lead.bg

      test$testcontact[test$bg] <- NA

      if(any(!test$testcontact, na.rm = T)){

        no.contact.pos <- which(!test$testcontact)

        suspect.lines$incorrect.contact <- no.contact.pos

        if(isTRUE(warn)){

          warning("Some measurements are not directly followed ",
                  "by blank measurements, at lines: ",
                  paste(no.contact.pos, collapse = ", "))

        }

      }

    }

  }

  if(isTRUE(id.contact)){

    if(method[1] == "approx" | method[1] == "mean"){

      test$testcontactid <- test$id == test$lag.id & test$id == test$lead.id

      test$testcontactid[test$bg] <- NA

      if(any(!test$testcontactid, na.rm = T)){

        no.contact.id.pos <- which(!test$testcontactid)

        suspect.lines$incorrect.id.contact <- no.contact.id.pos

        if(isTRUE(warn)){

          warning("Some measurements are not directly surrounded ",
                  "by blank measurements of similar id, at lines: ",
                  paste(no.contact.id.pos, collapse = ", "))

        }

      }

    } else if(method[1] == "previous"){

      test$testcontactid <- test$id == test$lag.id

      test$testcontactid[test$bg] <- NA

      if(any(!test$testcontactid, na.rm = T)){

        no.contact.id.pos <- which(!test$testcontactid)

        suspect.lines$incorrect.id.contact <- no.contact.id.pos

        if(isTRUE(warn)){

          warning("Some measurements are not directly preceded ",
                  "by blank measurements of similar id, at lines: ",
                  paste(no.contact.id.pos, collapse = ", "))

        }

      }

    } else if(method[1] == "next"){

      test$testcontactid <- test$id == test$lead.id

      test$testcontactid[test$bg] <- NA

      if(any(!test$testcontactid, na.rm = T)){

        no.contact.id.pos <- which(!test$testcontactid)

        suspect.lines$incorrect.id.contact <- no.contact.id.pos

        if(isTRUE(warn)){

          warning("Some measurements are not directly followed ",
                  "by blank measurements of similar id, at lines: ",
                  paste(no.contact.id.pos, collapse = ", "))

        }

      }

    }

  }

  # Other problems should be detectable;

  # Test range of background values

  if(!is.null(bglim)){

    if(length(bglim) == 1) bglim <- c(-bglim, bglim)

    db <- data.frame(id = id, bg = bg, x = x, t = t)

    test.bg <- db$bg & (db$x > max(bglim) | db$x < min(bglim))

    if(any(test.bg)){

      outlim.pos <- which(test.bg)

      suspect.lines$excessive.background <- outlim.pos

      if(isTRUE(warn)){

        warning("Some background measurements are outside the background limit",
                " range ('bglim'), at the lines ",
                paste(outlim.pos, collapse = ", "))

      }

    }

  }

  # Apply correction ----

  summary <- data.frame(id = id, x = x, bg = bg, t = t)

  summary$low  <- cumsum(as.integer(bg))
  summary$high <- rev(cumsum(rev(as.integer(bg))))

  background <- summary[summary$bg, ]

  summary$low.x <- background[summary$low,]$x
  summary$low.t <- background[summary$low,]$t

  summary$low.x[summary$bg] <- NA
  summary$low.t[summary$bg] <- NA

  rev.background <- background[rev(seq(nrow(background))), ]

  summary$high.x <- rev.background[summary$high, ]$x
  summary$high.t <- rev.background[summary$high, ]$t

  summary$high.x[summary$bg] <- NA
  summary$high.t[summary$bg] <- NA

  # Test duration within correction scheme ----

  if(!is.null(tlim)){

    if(method[1] == "approx" | method[1] == "mean"){

      summary$laps <- summary$high.t - summary$low.t

    } else if(method[1] == "previous"){

      summary$laps <- summary$t - summary$low.t

    } else if(method[1] == "next"){

      summary$laps <- summary$high.t - summary$t

    }

    out.time <- summary$laps > tlim

    if(any(out.time)){

      time.pos <- which(out.time)

      suspect.lines$excessive.time <- time.pos

      if(isTRUE(warn)){

        if(method[1] == "approx" | method[1] == "mean"){

          warning("Too much time (> ", tlim,
                  ") has elapsed between background measures at lines ",
                  paste(time.pos, collapse = ", "))

        } else {

          warning("Too much time (> ", tlim,
                  ") has elapsed between the background measurement and",
                  " the actual one, at lines ",
                  paste(time.pos, collapse = ", "))

        }

      }

    }

  }

  # Actual correction ----

  summary$initial <- seq(nrow(summary))

  if(method[1] == "approx"){

    xbg <- summary$low.x +
      ((summary$t - summary$low.t) *
         (( summary$high.x - summary$low.x )/( summary$high.t - summary$low.t )))

  } else if(method[1] == "mean"){

    xbg <- (summary$low.x + summary$high.x)/2

  } else if(method[1] == "previous"){

    xbg <-summary$ low.x

  } else if(method[1] == "next"){

    xbg <- summary$high.x

  }

  summary$xbg <- xbg

  summary$xcor <- summary$x - summary$xbg

  if(isTRUE(compress)) out <- summary[!summary$bg, ] else out <- summary

  outlist <- list(line = out$initial,
                  id = out$id,
                  t = out$t,
                  original.x = out$x,
                  corrected.x = out$xcor,
                  background.x = out$xbg,
                  previous.t = out$low.t,
                  previous.x = out$low.x,
                  next.t = out$high.t,
                  next.x = out$high.x,
                  suspect.lines = suspect.lines)

  return(outlist)

}
