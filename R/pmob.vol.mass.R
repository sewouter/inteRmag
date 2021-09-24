#' @title Check and fill a pmob for volume- and mass-related unit conversion.
#'
#' @description Check and fill a pmob (\strong{P}aleo\strong{M}agnetic
#' \strong{OB}ject) for volume- and mass-related unit conversion. THIS IS
#' CURRENTLY ONLY FOR DISCRETE SAMPLES, LONG CORES PARAMETERS ARE NOT TAKEN INTO
#' ACCOUNT.
#'
#' @param pmob The pmob to treat.
#' @param convert.mag Whether to convert the magnetic moment/magnetisation
#' units.
#' @param convert.magsus Whether to convert the magnetic susceptibility units.
#' @param digits The number of digits to consider to check the conversion
#' calculations to already provided values (when the same measurement is
#' provided in different units).
#' @param vol.main In the unique case where only density is available in the
#' pmob, this parameter details whether the mass will be calculated from the
#' default volume value (vol.main = TRUE), or on the contrary the volume will
#' be calculated from the default mass value (vol.main = FALSE).
#' @param default.vol The default value of volume. We mention a few standard
#' values;
#' 11.9 cm3 for the IODP standard (2.5 * 2.5 * 1.9 cm),
#' 10.1 \out{cm<sup>3</sup>} for the Bremen autosampler (2.3 * 2.3 * 1.9 cm),
#' 12.9 \out{cm<sup>3</sup>} for a 1 inch (2.54 cm) slice of a 1 inch-diameter drill core,
#' 10.1 \out{cm<sup>3</sup>} for a 2 cm slice of a 1 inch-diameter drill core.
#' 10.8 \out{cm<sup>3</sup>} for a 2.2 cm slice of a 2.5 cm-diameter drill core (CIMAN-ALP standard, or even more universal ??).
#' The default is here taken as 10.8 \out{cm<sup>3</sup>}, identical to the .ppl
#' files (on the arbitrary choice that this is the format that allows to link
#' inteRmag to PuffinPlot, a WYSIWYG paleomagnetic software that is a good fit
#' to be combined with our command-line scripting approach). Units are in
#' \out{m<sup>3</sup>}.
#' @param default.mass The default mass. This needs to be coherent with the
#' default volume and the default density. Units are in kg.
#' @param default.dens The default density of the sample, that allows to compute
#' either mass or volume based on the other. The default value is set as 2.7 g
#' per cm3, based on the assumed surface density of 2.7 \out{g cm<sup>-3</sup>},
#' mimicking the density of granite or granodiorite (Williamson & Adams, 1923).
#' Units are in \out{kg m<sup>-3</sup>}.
#' @param vol.lim The maximum volume expected for each sample; any volume value
#' (calculated or provided) exceeding it triggers a warning. This is set by
#' default to 13 \out{cm<sup>3</sup>}.
#' @param mass.lim The maximum mass expected for each sample; any mass value
#' (calculated or provided) exceeding it triggers a warning.
#' @param dens.lim The range of density expected for each sample. If only one
#' value is provided, this parameter sets the maximum density value expected for
#' the sample. Any discrepancy triggers a warning.
#'
#' @details If two values are provided for the same magnetisation measurement,
#' (e.g. if a measurement is provided with both the magnetic moment in Am2 and
#' in magnetisation by volume), the unnormalised values (xint, yint, zint, and
#' totmagsus) will be used in priority to compute the missing measurements, and
#' then the volume-normalised ones (xvol, yvol, zvol, and volmagsus) (in our
#' example, the mass normalised magnetisation will be obtained from the magnetic
#' moment). In theory, this should only be relevant in the rounding used on the
#' values.
#'
#' @references
#' \itemize{
#'   \item Williamson, E. D., and Adams, L. H., 1923. Density distribution in
#'   the Earth. Journal of the Washington Academy of Sciences 13, n°19, 413-428.
#' }
#'
#' @examples
#' pmob <- as.pmob(sampid = c("A1", "B2", "C3"),
#'                 xint = c(200, 450, NA),
#'                 yint = c(60,  30,  NA),
#'                 zint = c(10,  20,  NA),
#'                 xvol = c(NA, NA, 800000),
#'                 yvol = c(NA, NA,  60000),
#'                 zvol = c(NA, NA,  30000),
#'                 mass = c(0.020, 0.027, 0.010))
#'
#' pmob.vol.mass(pmob, convert.magsus = F)
#'
#' @export

pmob.vol.mass <- function(pmob,
                          convert.mag    = TRUE,
                          convert.magsus = TRUE,
                          digits = 5,
                          vol.main = TRUE,
                          default.vol  = 10.8  * 1e-6,
                          default.mass = 29.16 * 1e-3,
                          default.dens = 2.7   * 1e3,
                          vol.lim      = 13    * 1e-6,
                          mass.lim     = NULL,
                          dens.lim     = NULL)
{

  if(!is.pmob(pmob)) stop("Incorrect pmob object")

  check.default <- all.equal(default.mass/default.vol, default.dens)

  if(!check.default){
    stop("Dividing 'default.mass' by 'default.vol' ",
         "should give 'default.dens'.")
  }

  if(length(dens.lim) == 1) {
    dens.lim <- c(0, dens.lim)
  } else if(!is.null(dens.lim)) {
    dens.lim <- range(dens.lim)
  }

  spmob <- pmob

  if(is.null(spmob$isblank))  spmob<- pmob.add(pmob = spmob, isblank = F)
  if(is.null(spmob$isholder)) spmob <- pmob.add(pmob = spmob, isholder = F)

  is.in <- !spmob$isblank & !spmob$isholder

  select <- which(is.in)
  is.out <- which(!is.in)

  epmob  <- pmob[select,]

  if(is.null(pmob$vol))  epmob <- pmob.add(pmob = epmob, vol = NA)
  if(is.null(pmob$mass)) epmob <- pmob.add(pmob = epmob, mass = NA)
  if(is.null(pmob$dens)) epmob <- pmob.add(pmob = epmob, dens = NA)

  if(is.null(pmob$exactvol)) {
    exactvol                    <- rep(F, nrow(epmob))
    exactvol[!is.na(epmob$vol)] <- T

    epmob <- pmob.add(pmob = epmob, exactvol = exactvol)
  }

  if(is.null(pmob$exactmass)) {
    exactmass                     <- rep(F, nrow(epmob))
    exactmass[!is.na(epmob$mass)] <- T

    epmob <- pmob.add(pmob = epmob, exactmass = exactmass)
  }

  if(is.null(pmob$exactdens)) {
    exactdens                     <- rep(F, nrow(epmob))
    exactdens[!is.na(epmob$dens)] <- T

    epmob <- pmob.add(pmob = epmob, exactdens = exactdens)
  }

  npmob <- epmob

  # Identify lines ----

  imass <- is.na(epmob$mass)
  ivol  <- is.na(epmob$vol)
  idens <- is.na(epmob$dens)

  miss3    <- imass  & ivol  & idens  # all 3 missing

  onlymass <- !imass & ivol  & idens  # only mass known
  onlyvol  <- imass  & !ivol & idens  # only vol known
  onlydens <- imass  & ivol  & !idens # only dens known

  massvol  <- !imass & !ivol & idens  # mass and vol known
  massdens <- !imass & ivol  & !idens # mass and dens known
  voldens  <- imass  & !ivol & !idens # vol and dens known

  all3     <- !imass & !ivol & !idens # all 3 known

  # Add default values when no other choice ----

  npmob$vol[miss3 | onlydens] <- default.vol
  npmob$mass[miss3]           <- default.mass


  npmob$dens[miss3 | onlymass | onlyvol] <- default.dens

  # Calculate values when possible ----

  volcal  <- onlymass | onlydens | massdens
  masscal <- onlyvol  | onlydens | voldens

  npmob$vol[volcal]   <- npmob$mass[volcal] / npmob$dens[volcal]
  npmob$mass[masscal] <- npmob$vol[masscal] * npmob$dens[masscal]
  npmob$dens[massvol] <- npmob$mass[massvol] / npmob$vol[massvol]

  # Give incentive to good calculations ----

  npmob$exactvol[massdens] <- T
  npmob$exactmass[voldens] <- T
  npmob$exactdens[massvol] <- T

  if(!all.equal(round(npmob$mass/npmob$vol, 11), round(npmob$dens, 11))){

    prob.lines <- which(round(npmob$mass/npmob$vol, 11) != round(npmob$dens, 11))

    stop("$mass/$vol is not equal to $density at the following lines: ",
         paste(select[prob.lines], collapse = ", "))

  }

  if(!is.null(vol.lim)){

    excess.vol <- npmob$vol > vol.lim

    if(any(excess.vol)){

      max.vol  <- max(npmob$vol)
      prob.vol <- which(excess.vol)

      warning("Excessive volume values detected (max = ",
              sprintf("%.2E", max.vol),
              " m3, or ",
              round(max.vol * 10^6, 1),
              " cm3, when the limit is at ",
              sprintf("%.2E", vol.lim),
              " m3, or ",
              round(vol.lim * 10^6, 1),
              " cm3), at lines: \n",
              paste(select[prob.vol], collapse = ", "))

    }

  }

  if(!is.null(mass.lim)){

    excess.mass <- npmob$mass > mass.lim

    if(any(excess.mass)){

      max.mass  <- max(npmob$mass)
      prob.mass <- which(excess.mass)

      warning("Excessive mass values detected (max = ",
              round(max.mass, 4),
              " kg, or ",
              round(max.mass * 10^3, 1),
              " g, when the limit is at ",
              round(mass.lim, 4),
              " kg, or ",
              round(mass.lim * 10^3, 1),
              " g), at lines: \n",
              paste(select[prob.mass], collapse = ", "))

    }

  }

  if(!is.null(dens.lim)){

    out.dens <- npmob$dens > max(dens.lim) | npmob$dens < min(dens.lim)

    if(any(out.dens)){

      min.dens <- min(npmob$dens)
      max.dens <- max(npmob$dens)

      prob.dens <- which(out.dens)

      warning("Density values detected outside the expected range of ",
              round(min(dens.lim)), "-", round(max(dens.lim)), " kg/m3 or ",
              round(min(dens.lim) * 0.001, 2), "-",
              round(max(dens.lim) * 0.001, 2), " g/cm3 (actual range: ",
              round(min.dens), "-", round(max.dens), " kg/m3 or ",
              round(min.dens * 0.001, 2), "-",
              round(max.dens * 0.001, 2), " g/cm3), at lines: ",
              paste(select[prob.dens], collapse = ", "))
    }

  }

  if(isTRUE(convert.mag)){

    if(is.null(pmob$xint)) npmob <- pmob.add(pmob = npmob, xint = NA)
    if(is.null(pmob$yint)) npmob <- pmob.add(pmob = npmob, yint = NA)
    if(is.null(pmob$zint)) npmob <- pmob.add(pmob = npmob, zint = NA)

    if(is.null(pmob$xvol)) npmob <- pmob.add(pmob = npmob, xvol = NA)
    if(is.null(pmob$yvol)) npmob <- pmob.add(pmob = npmob, yvol = NA)
    if(is.null(pmob$zvol)) npmob <- pmob.add(pmob = npmob, zvol = NA)

    if(is.null(pmob$xmass)) npmob <- pmob.add(pmob = npmob, xmass = NA)
    if(is.null(pmob$ymass)) npmob <- pmob.add(pmob = npmob, ymass = NA)
    if(is.null(pmob$zmass)) npmob <- pmob.add(pmob = npmob, zmass = NA)

    na.xint <- is.na(npmob$xint)
    na.yint <- is.na(npmob$yint)
    na.zint <- is.na(npmob$zint)

    na.xvol <- is.na(npmob$xvol)
    na.yvol <- is.na(npmob$yvol)
    na.zvol <- is.na(npmob$zvol)

    na.xmass <- is.na(npmob$xmass)
    na.ymass <- is.na(npmob$ymass)
    na.zmass <- is.na(npmob$zmass)

    bloc.int  <- as.integer(na.xint)  + as.integer(na.yint)  + as.integer(na.zint)
    bloc.vol  <- as.integer(na.xvol)  + as.integer(na.yvol)  + as.integer(na.zvol)
    bloc.mass <- as.integer(na.xmass) + as.integer(na.ymass) + as.integer(na.zmass)

    if(any(bloc.int == 1 | bloc.int == 2)){
      prob.bloc.int.pos <- which(bloc.int == 1 | bloc.int == 2)

      warning("Missing xint, yint or zint values (one or two of the three)",
              " at lines: ",
              paste(select[prob.bloc.pos], collapse = ", "))
    }

    if(any(bloc.vol == 1 | bloc.vol == 2)){
      prob.bloc.vol.pos <- which(bloc.vol == 1 | bloc.vol == 2)

      warning("Missing xvol, yvol or zvol values (one or two of the three)",
              " at lines: ",
              paste(select[prob.bloc.vol.pos], collapse = ", "))
    }

    if(any(bloc.mass == 1 | bloc.mass == 2)){
      prob.bloc.mass.pos <- which(bloc.mass == 1 | bloc.mass == 2)

      warning("Missing xmass, ymass or zmass values (one or two of the three)",
              " at lines: ",
              paste(select[prob.bloc.mass.pos], collapse = ", "))
    }

    test <- data.frame(xint =  npmob$xint,  yint =  npmob$yint,  zint = npmob$zint,
                       xvol =  npmob$xvol,  yvol =  npmob$yvol,  zvol = npmob$zvol,
                       xmass = npmob$xmass, ymass = npmob$ymass, zmass = npmob$zmass,
                       xint.from.xvol = npmob$xvol * npmob$vol,
                       yint.from.yvol = npmob$yvol * npmob$vol,
                       zint.from.zvol = npmob$zvol * npmob$vol,
                       xint.from.xmass = npmob$xmass * npmob$mass,
                       yint.from.ymass = npmob$ymass * npmob$mass,
                       zint.from.zmass = npmob$zmass * npmob$mass,
                       xvol.from.xint = npmob$xint / npmob$vol,
                       yvol.from.yint = npmob$yint / npmob$vol,
                       zvol.from.zint = npmob$zint / npmob$vol,
                       xvol.from.xmass = (npmob$xmass * npmob$mass) / npmob$vol,
                       yvol.from.ymass = (npmob$ymass * npmob$mass) / npmob$vol,
                       zvol.from.zmass = (npmob$zmass * npmob$mass) / npmob$vol,
                       xmass.from.xint = npmob$xint / npmob$mass,
                       ymass.from.yint = npmob$yint / npmob$mass,
                       zmass.from.zint = npmob$zint / npmob$mass,
                       xmass.from.xvol = (npmob$xvol * npmob$vol) / npmob$mass,
                       ymass.from.yvol = (npmob$yvol * npmob$vol) / npmob$mass,
                       zmass.from.zvol = (npmob$zvol * npmob$vol) / npmob$mass)

    charfun <- function(x) sprintf(paste0("%.", digits,"E"), x)

    comp <- data.frame(xint = charfun(test$xint),
                       yint = charfun(test$yint),
                       zint = charfun(test$zint),
                       xvol = charfun(test$xvol),
                       yvol = charfun(test$yvol),
                       zvol = charfun(test$zvol),
                       xmass = charfun(test$xmass),
                       ymass = charfun(test$ymass),
                       zmass = charfun(test$zmass),
                       xint.from.xvol = charfun(test$xint.from.xvol),
                       yint.from.yvol = charfun(test$yint.from.yvol),
                       zint.from.zvol = charfun(test$zint.from.zvol),
                       xint.from.xmass = charfun(test$xint.from.xmass),
                       yint.from.ymass = charfun(test$yint.from.ymass),
                       zint.from.zmass = charfun(test$zint.from.zmass),
                       xvol.from.xint = charfun(test$xvol.from.xint),
                       yvol.from.yint = charfun(test$yvol.from.yint),
                       zvol.from.zint = charfun(test$zvol.from.zint),
                       xvol.from.xmass = charfun(test$xvol.from.xmass),
                       yvol.from.ymass = charfun(test$yvol.from.ymass),
                       zvol.from.zmass = charfun(test$zvol.from.zmass),
                       xmass.from.xint = charfun(test$xmass.from.xint),
                       ymass.from.yint = charfun(test$ymass.from.yint),
                       zmass.from.zint = charfun(test$zmass.from.zint),
                       xmass.from.xvol = charfun(test$xmass.from.xvol),
                       ymass.from.yvol = charfun(test$ymass.from.yvol),
                       zmass.from.zvol = charfun(test$zmass.from.zvol))

    # ----

    comp$int.vol.x  <- comp$xint == comp$xint.from.xvol
    comp$int.vol.y  <- comp$yint == comp$yint.from.yvol
    comp$int.vol.z  <- comp$zint == comp$zint.from.zvol

    comp$int.mass.x <- comp$xint == comp$xint.from.xmass
    comp$int.mass.y <- comp$yint == comp$yint.from.ymass
    comp$int.mass.z <- comp$zint == comp$zint.from.zmass

    comp$vol.mass.x <- comp$xvol == comp$xvol.from.xmass
    comp$vol.mass.y <- comp$yvol == comp$yvol.from.ymass
    comp$vol.mass.z <- comp$zvol == comp$zvol.from.zmass

    # ----

    comp$prob.int.vol.x <- comp$int.vol.x & !is.na(test$xint) & !is.na(test$xvol)
    comp$prob.int.vol.y <- comp$int.vol.y & !is.na(test$yint) & !is.na(test$yvol)
    comp$prob.int.vol.z <- comp$int.vol.z & !is.na(test$zint) & !is.na(test$zvol)

    comp$prob.int.mass.x <- comp$int.mass.x & !is.na(test$xint) & !is.na(test$xmass)
    comp$prob.int.mass.y <- comp$int.mass.y & !is.na(test$yint) & !is.na(test$ymass)
    comp$prob.int.mass.z <- comp$int.mass.z & !is.na(test$zint) & !is.na(test$zmass)

    comp$prob.vol.mass.x <- comp$vol.mass.x & !is.na(test$xvol) & !is.na(test$xmass)
    comp$prob.vol.mass.y <- comp$vol.mass.y & !is.na(test$yvol) & !is.na(test$ymass)
    comp$prob.vol.mass.z <- comp$vol.mass.z & !is.na(test$zvol) & !is.na(test$zmass)

    # ----

    if(any(comp$prob.int.vol.x)){
      prob.int.vol.x.pos <- which(comp$prob.int.vol.x)
      warning("xint cannot be obtained correctly from xvol at lines: ",
              paste(select[prob.int.vol.x.pos], collapse = ", "))
    }

    if(any(comp$prob.int.vol.y)){
      prob.int.vol.y.pos <- which(comp$prob.int.vol.y)
      warning("yint cannot be obtained correctly from yvol at lines: ",
              paste(select[prob.int.vol.y.pos], collapse = ", "))
    }

    if(any(comp$prob.int.vol.z)){
      prob.int.vol.z.pos <- which(comp$prob.int.vol.z)
      warning("zint cannot be obtained correctly from zvol at lines: ",
              paste(select[prob.int.vol.z.pos], collapse = ", "))
    }

    # ----

    if(any(comp$prob.int.mass.x)){
      prob.int.mass.x.pos <- which(comp$prob.int.mass.x)
      warning("xint cannot be obtained correctly from xmass at lines: ",
              paste(select[prob.int.mass.x.pos], collapse = ", "))
    }

    if(any(comp$prob.int.mass.y)){
      prob.int.mass.y.pos <- which(comp$prob.int.mass.y)
      warning("yint cannot be obtained correctly from ymass at lines: ",
              paste(select[prob.int.mass.y.pos], collapse = ", "))
    }

    if(any(comp$prob.int.mass.z)){
      prob.int.mass.z.pos <- which(comp$prob.int.mass.z)
      warning("zint cannot be obtained correctly from zmass at lines: ",
              paste(select[prob.int.mass.z.pos], collapse = ", "))
    }

    # ----

    if(any(comp$prob.vol.mass.x)){
      prob.vol.mass.x.pos <- which(comp$prob.vol.mass.x)
      warning("xvol cannot be obtained correctly from xmass at lines: ",
              paste(select[prob.vol.mass.x.pos], collapse = ", "))
    }

    if(any(comp$prob.vol.mass.y)){
      prob.vol.mass.y.pos <- which(comp$prob.vol.mass.y)
      warning("yvol cannot be obtained correctly from ymass at lines: ",
              paste(select[prob.vol.mass.y.pos], collapse = ", "))
    }

    if(any(comp$prob.vol.mass.z)){
      prob.vol.mass.z.pos <- which(comp$prob.vol.mass.z)
      warning("zvol cannot be obtained correctly from zmass at lines: ",
              paste(select[prob.vol.mass.z.pos], collapse = ", "))
    }

    # Convert
    # SAY IN DESCRIPTION THE ORDER OF PRIORITY -----------------------------------------------------

    npmob$xint[na.xint & !na.xvol] <- test$xint.from.xvol[na.xint & !na.xvol]
    npmob$yint[na.yint & !na.yvol] <- test$yint.from.yvol[na.yint & !na.yvol]
    npmob$zint[na.zint & !na.zvol] <- test$zint.from.zvol[na.zint & !na.zvol]

    npmob$xint[na.xint & na.xvol & !na.xmass] <- test$xint.from.xmass[na.xint & na.xvol & !na.xmass]
    npmob$yint[na.yint & na.yvol & !na.ymass] <- test$yint.from.ymass[na.yint & na.yvol & !na.ymass]
    npmob$zint[na.zint & na.zvol & !na.zmass] <- test$zint.from.zmass[na.zint & na.zvol & !na.zmass]

    # ----

    npmob$xvol[!na.xint & na.xvol] <- test$xvol.from.xint[!na.xint & na.xvol]
    npmob$yvol[!na.yint & na.yvol] <- test$yvol.from.yint[!na.yint & na.yvol]
    npmob$zvol[!na.zint & na.zvol] <- test$zvol.from.zint[!na.zint & na.zvol]

    npmob$xvol[na.xint & na.xvol & !na.xmass] <- test$xvol.from.xmass[na.xint & na.xvol & !na.xmass]
    npmob$yvol[na.yint & na.yvol & !na.ymass] <- test$yvol.from.ymass[na.yint & na.yvol & !na.ymass]
    npmob$zvol[na.zint & na.zvol & !na.zmass] <- test$zvol.from.zmass[na.zint & na.zvol & !na.zmass]

    # ----

    npmob$xmass[!na.xint & na.xmass] <- test$xmass.from.xint[!na.xint & na.xmass]
    npmob$ymass[!na.yint & na.ymass] <- test$ymass.from.yint[!na.yint & na.ymass]
    npmob$zmass[!na.zint & na.zmass] <- test$zmass.from.zint[!na.zint & na.zmass]

    npmob$xmass[na.xint & na.xmass & !na.xvol] <- test$xmass.from.xvol[na.xint & na.xmass & !na.xvol]
    npmob$ymass[na.yint & na.ymass & !na.yvol] <- test$ymass.from.yvol[na.yint & na.ymass & !na.yvol]
    npmob$zmass[na.zint & na.zmass & !na.zvol] <- test$zmass.from.zvol[na.zint & na.zmass & !na.zvol]

  }

  # ----

  if(isTRUE(convert.magsus)){

    if(is.null(pmob$totmagsus))  npmob <- pmob.add(pmob = npmob, totmagsus = NA)
    if(is.null(pmob$volmagsus))  npmob <- pmob.add(pmob = npmob, volmagsus = NA)
    if(is.null(pmob$massmagsus)) npmob <- pmob.add(pmob = npmob, massmagsus = NA)

    na.totms  <- is.na(npmob$totmagsus)
    na.volms  <- is.na(npmob$volmagsus)
    na.massms <- is.na(npmob$massmagsus)

    test.ms <- data.frame(totms = npmob$totmagsus,
                          volms = npmob$volmagsus,
                          massms = npmob$massmagsus,
                          totms.from.volms = npmob$volmagsus * npmob$vol,
                          totms.from.massms = npmob$massmagsus * npmob$mass,
                          volms.from.totms = npmob$totmagsus / npmob$vol,
                          volms.from.massms =
                            (npmob$massmagsus * npmob$mass) / npmob$vol,
                          massms.from.totms = npmob$totmagsus / npmob$mass,
                          massms.from.volms =
                            (npmob$volmagsus * npmob$vol) / npmob$mass)

    charfun <- function(x) sprintf(paste0("%.", digits,"E"), x)

    comp.ms <- data.frame(totms = charfun(test.ms$totms),
                          volms = charfun(test.ms$volms),
                          massms = charfun(test.ms$massms),
                          totms.from.volms = charfun(test.ms$totms.from.volms),
                          totms.from.massms =charfun(test.ms$totms.from.massms),
                          volms.from.totms =charfun(test.ms$volms.from.totms),
                          volms.from.massms = charfun(test.ms$volms.from.massms),
                          massms.from.totms = charfun(test.ms$massms.from.totms),
                          massms.from.volms = charfun(test.ms$massms.from.volms))

    comp.ms$tot.vol.ms  <- comp.ms$totms == comp.ms$totms.from.volms
    comp.ms$tot.mass.ms <- comp.ms$totms == comp.ms$totms.from.massms
    comp.ms$vol.mass.ms <- comp.ms$volms == comp.ms$volms.from.massms

    # ----

    comp.ms$prob.tot.vol.ms <- comp.ms$tot.vol.ms &
      !is.na(test.ms$totms) &
      !is.na(test.ms$volms)

    comp.ms$prob.tot.mass.ms <- comp.ms$tot.mass.ms &
      !is.na(test.ms$totms) &
      !is.na(test.ms$massms)

    comp.ms$prob.vol.mass.ms <- comp.ms$vol.mass.ms &
      !is.na(test.ms$volms) &
      !is.na(test.ms$massms)

    if(any(comp.ms$prob.tot.vol.ms)){
      prob.tot.vol.ms.pos <- which(comp.ms$prob.tot.vol.ms)
      warning("totmagsus cannot be obtained correctly from volmagsus at lines: ",
              paste(select[prob.tot.vol.ms.pos], collapse = ", "))
    }

    if(any(comp.ms$prob.tot.mass.ms)){
      prob.tot.mass.ms.pos <- which(comp.ms$prob.tot.mass.ms)
      warning("totmagsus cannot be obtained correctly from massmagsus at lines: ",
              paste(select[prob.tot.mass.ms.pos], collapse = ", "))
    }

    if(any(comp.ms$prob.vol.mass.ms)){
      prob.vol.mass.ms.pos <- which(comp.ms$prob.vol.mass.ms)
      warning("volmagsus cannot be obtained correctly from massmagsus at lines: ",
              paste(select[prob.vol.mass.ms.pos], collapse = ", "))
    }

    # Convert
    # SAY IN DESCRIPTION THE ORDER OF PRIORITY -------------------------------------------------

    npmob$totmagsus[na.totms & !na.volms] <- test.ms$totms.from.volms[na.totms & !na.volms]
    npmob$totmagsus[na.totms & na.volms & !na.massms] <- test.ms$totms.from.massms[na.totms & na.volms & !na.massms]

    # ----

    npmob$volmagsus[!na.totms & na.volms] <- test.ms$volms.from.totms[!na.totms & na.volms]
    npmob$volmagsus[na.totms & na.volms & !na.massms] <- test.ms$volms.from.massms[na.totms & na.volms & !na.massms]

    # ----

    npmob$massmagsus[!na.totms & na.massms] <- test.ms$massms.from.totms[!na.totms & na.massms]
    npmob$massmagsus[na.totms & !na.volms & na.massms] <- test.ms$massms.from.volms[na.totms & !na.volms & na.massms]

  }

  if(isTRUE(convert.mag) & isTRUE(convert.magsus)){

    remove <- colnames(pmob) %in% c("xint",  "yint",  "zint",
                                    "xvol",  "yvol",  "zvol",
                                    "xmass", "ymass", "zmass",
                                    "totmagsus", "volmagsus", "massmagsus",
                                    "vol",   "mass",  "dens")

    out <- pmob[, !remove, drop = F]
    out <- pmob.add(pmob = out,
                    xint = NA,  yint = NA,  zint = NA,
                    xvol = NA,  yvol = NA,  zvol = NA,
                    xmass = NA, ymass = NA, zmass = NA,
                    totmagsus = NA, volmagsus = NA, massmagsus = NA,
                    vol = NA,   mass = NA,  dens = NA)

    out$xint[select]       <- npmob$xint
    out$yint[select]       <- npmob$yint
    out$zint[select]       <- npmob$zint
    out$xvol[select]       <- npmob$xvol
    out$yvol[select]       <- npmob$yvol
    out$zvol[select]       <- npmob$zvol
    out$xmass[select]      <- npmob$xmass
    out$ymass[select]      <- npmob$ymass
    out$zmass[select]      <- npmob$zmass
    out$totmagsus[select]  <- npmob$totmagsus
    out$volmagsus[select]  <- npmob$volmagsus
    out$massmagsus[select] <- npmob$massmagsus
    out$vol[select]        <- npmob$vol
    out$mass[select]       <- npmob$mass
    out$dens[select]       <- npmob$dens

    if(!is.null(pmob$xint)) out$xint[is.out] <- pmob$xint[is.out]
    if(!is.null(pmob$yint)) out$yint[is.out] <- pmob$yint[is.out]
    if(!is.null(pmob$zint)) out$xint[is.out] <- pmob$zint[is.out]
    if(!is.null(pmob$xvol)) out$xvol[is.out] <- pmob$xvol[is.out]
    if(!is.null(pmob$yvol)) out$yvol[is.out] <- pmob$yvol[is.out]
    if(!is.null(pmob$zvol)) out$xvol[is.out] <- pmob$zvol[is.out]
    if(!is.null(pmob$xmass)) out$xmass[is.out] <- pmob$xmass[is.out]
    if(!is.null(pmob$ymass)) out$ymass[is.out] <- pmob$ymass[is.out]
    if(!is.null(pmob$zmass)) out$xmass[is.out] <- pmob$zmass[is.out]

    if(!is.null(pmob$totmagsus)){
      out$totmagsus[is.out]  <- pmob$totmagsus[is.out]
    }

    if(!is.null(pmob$volmagsus)){
      out$volmagsus[is.out]  <- pmob$volmagsus[is.out]
    }

    if(!is.null(pmob$massmagsus)){
      out$massmagsus[is.out] <- pmob$massmagsus[is.out]
    }

    if(!is.null(pmob$vol))  out$vol[is.out]  <- pmob$vol[is.out]
    if(!is.null(pmob$mass)) out$mass[is.out] <- pmob$mass[is.out]
    if(!is.null(pmob$dens)) out$dens[is.out] <- pmob$dens[is.out]

  } else if(isTRUE(convert.mag) & !isTRUE(convert.magsus)) {

    remove <- colnames(pmob) %in% c("xint",  "yint",  "zint",
                                    "xvol",  "yvol",  "zvol",
                                    "xmass", "ymass", "zmass",
                                    "vol",   "mass",  "dens")

    out <- pmob[, !remove, drop = F]
    out <- pmob.add(pmob = out,
                    xint = NA,  yint = NA,  zint = NA,
                    xvol = NA,  yvol = NA,  zvol = NA,
                    xmass = NA, ymass = NA, zmass = NA,
                    vol = NA,   mass = NA,  dens = NA)

    out$xint[select]  <- npmob$xint
    out$yint[select]  <- npmob$yint
    out$zint[select]  <- npmob$zint
    out$xvol[select]  <- npmob$xvol
    out$yvol[select]  <- npmob$yvol
    out$zvol[select]  <- npmob$zvol
    out$xmass[select] <- npmob$xmass
    out$ymass[select] <- npmob$ymass
    out$zmass[select] <- npmob$zmass
    out$vol[select]   <- npmob$vol
    out$mass[select]  <- npmob$mass
    out$dens[select]  <- npmob$dens

    if(!is.null(pmob$xint)) out$xint[is.out] <- pmob$xint[is.out]
    if(!is.null(pmob$yint)) out$yint[is.out] <- pmob$yint[is.out]
    if(!is.null(pmob$zint)) out$xint[is.out] <- pmob$zint[is.out]
    if(!is.null(pmob$xvol)) out$xvol[is.out] <- pmob$xvol[is.out]
    if(!is.null(pmob$yvol)) out$yvol[is.out] <- pmob$yvol[is.out]
    if(!is.null(pmob$zvol)) out$xvol[is.out] <- pmob$zvol[is.out]
    if(!is.null(pmob$xmass)) out$xmass[is.out] <- pmob$xmass[is.out]
    if(!is.null(pmob$ymass)) out$ymass[is.out] <- pmob$ymass[is.out]
    if(!is.null(pmob$zmass)) out$xmass[is.out] <- pmob$zmass[is.out]
    if(!is.null(pmob$vol))  out$vol[is.out]  <- pmob$vol[is.out]
    if(!is.null(pmob$mass)) out$mass[is.out] <- pmob$mass[is.out]
    if(!is.null(pmob$dens)) out$dens[is.out] <- pmob$dens[is.out]

  } else if(!isTRUE(convert.mag) & isTRUE(convert.magsus)) {

    remove <- colnames(pmob) %in% c("totmagsus", "volmagsus", "massmagsus",
                                    "vol",   "mass",  "dens")

    out <- pmob[, !remove, drop = F]
    out <- pmob.add(pmob = out,
                    totmagsus = NA, volmagsus = NA, massmagsus = NA,
                    vol = NA,   mass = NA,  dens = NA)

    out$totmagsus[select]  <- npmob$totmagsus
    out$volmagsus[select]  <- npmob$volmagsus
    out$massmagsus[select] <- npmob$massmagsus
    out$vol[select]        <- npmob$vol
    out$mass[select]       <- npmob$mass
    out$dens[select]       <- npmob$dens

    if(!is.null(pmob$totmagsus)){
      out$totmagsus[is.out]  <- pmob$totmagsus[is.out]
    }

    if(!is.null(pmob$volmagsus)){
      out$volmagsus[is.out]  <- pmob$volmagsus[is.out]
    }

    if(!is.null(pmob$massmagsus)){
      out$massmagsus[is.out] <- pmob$massmagsus[is.out]
    }

    if(!is.null(pmob$vol))  out$vol[is.out]  <- pmob$vol[is.out]
    if(!is.null(pmob$mass)) out$mass[is.out] <- pmob$mass[is.out]
    if(!is.null(pmob$dens)) out$dens[is.out] <- pmob$dens[is.out]


  } else if(!isTRUE(convert.mag) & !isTRUE(convert.magsus)) {

    remove <- colnames(pmob) %in% c("vol",   "mass",  "dens")

    out <- pmob[, !remove, drop = F]
    out <- pmob.add(pmob = out,
                    vol = NA,   mass = NA,  dens = NA)

    out$vol[select]  <- npmob$vol
    out$mass[select] <- npmob$mass
    out$dens[select] <- npmob$dens

    if(!is.null(pmob$vol))  out$vol[is.out]  <- pmob$vol[is.out]
    if(!is.null(pmob$mass)) out$mass[is.out] <- pmob$mass[is.out]
    if(!is.null(pmob$dens)) out$dens[is.out] <- pmob$dens[is.out]


  }

  return(out)

}





