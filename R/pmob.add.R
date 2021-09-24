#' @title Add new columns to a pmob
#'
#' @description Add new columns to a pmob (\strong{P}aleo\strong{M}agnetic
#' \strong{OB}jects).
#'
#' @param ... Personalized parameters (see \code{\link{as.pmob}}).
#' @param pmob pmob object to complement
#' @param sampid,specid,slotid,measid,siteid Identification parameters
#' (see \code{\link{as.pmob}}).
#' @param initialorder Initial order (see \code{\link{as.pmob}}).
#' @param depth Sampling depth (see \code{\link{as.pmob}}).
#' @param meassec,measmin,meashour,measday,measmonth,measyear
#' Time of the measurement (see \code{\link{as.pmob}}).
#' @param measdevice,magsusdevice Measurement devices
#' (see \code{\link{as.pmob}}).
#' @param xint,yint,zint Magnetic moments (see \code{\link{as.pmob}}).
#' @param blankcor,isblank,blankmethod Correction parameters (see \code{\link{as.pmob}}).
#' @param holdercor,isholder,holdermethod Correction parameters (see \code{\link{as.pmob}}).
#' @param xinterror,yinterror,zinterror Error (2σ) on the magnetic moment
#' (see \code{\link{as.pmob}}).
#' @param xvol,yvol,zvol Magnetization by volume (see \code{\link{as.pmob}}).
#' @param xmass,ymass,zmass Magnetization by mass (see \code{\link{as.pmob}}).
#' @param totmagsus,totmagsuserror,volmagsus,massmagsus Magnetic susceptibility
#' data (see \code{\link{as.pmob}}).
#' @param vol,mass,dens Sample properties (see \code{\link{as.pmob}}).
#' @param exactvol,exactmass,exactdens Sample properties reliability
#' (see \code{\link{as.pmob}}).
#' @param discrete Whether a single discrete sample is measured (TRUE), or a
#' continuous core (FALSE) (see \code{\link{as.pmob}}).
#' @param area Cross-sectional area of a long core sample (for continuous cores,
#' applies if the parameter \strong{discrete} is FALSE)
#' (see \code{\link{as.pmob}}).
#' @param defaultarea Long core sample properties reliability (see \code{\link{as.pmob}}).
#' @param xefflength,yefflength,zefflength Effective portion of length of a long
#' core that is measured by the sensor (see \code{\link{as.pmob}}).
#' @param xeffvol,yeffvol,zeffvol DEffective volume of a long core that is
#' measured by the sensors  (see \code{\link{as.pmob}}).
#' @param xeffmass,yeffmass,zeffmass Effective mass of a long core that is
#' measured by the sensors (see \code{\link{as.pmob}}).
#' @param magsuseffvol,magsuseffmass Magnetic susceptibilities of a long core
#' (see \code{\link{as.pmob}}).
#' @param deconvoluted Whether the continuous core measurements have been
#' deconvoluted (TRUE) or not (FALSE) (see \code{\link{as.pmob}}).
#' @param xintconv,yintconv,zintconv The original magnetic moment obtained
#' without deconvolution (see \code{\link{as.pmob}}).
#' @param sampaz,sampdip,samprot Sample orientation in the field
#' (see \code{\link{as.pmob}}).
#' @param coraz,cordip,corrot Sample orientation in the
#' measuring device (see \code{\link{as.pmob}}).
#' @param bedaz,beddip Orientation of beds (see \code{\link{as.pmob}}).
#' @param foldtrend,foldplunge Orientation of the folding
#' (see \code{\link{as.pmob}}).
#' @param magaz,usemagaz,solaraz Parameters determining the zimuth of the sample
#' (see \code{\link{as.pmob}}).
#' @param solarangle Solar angle (see \code{\link{as.pmob}}).
#' @param long,lat,height Coordinates (see \code{\link{as.pmob}}).
#' @param sampmin,samphour,sampday,sampmonth,sampyear
#' Time of sampling  (see \code{\link{as.pmob}}).
#' @param samptimezonemin,samptimezonehour Time zone of the sampling
#' (see \code{\link{as.pmob}}).
#' @param magvar,magazvarcor,bedazvarcor,foldtrendcarcor Parameters to
#' correct due to the magnetic variation (see \code{\link{as.pmob}}).
#' @param treatafx,treatafy,treatafz Treatment field by alternating field (AF)
#' demagnetization (see \code{\link{as.pmob}}).
#' @param treattemp Treatment temperature (see \code{\link{as.pmob}}).
#' @param treatirmx,treatirmy,treatirmz Treatment field by Isothermal Remanent
#' Magnetisation (IRM) (see \code{\link{as.pmob}}).
#' @param treatarmafx,treatarmafy,treatarmafz Treatment anhysteretic field
#'for Anhysteretic Remanent Magnetization (ARM) (see \code{\link{as.pmob}}).
#' @param treatarmbiasx,treatarmbiasy,treatarmbiasz Treatment bias field for
#' Anhysteretic Remanent Magnetization (ARM) (see \code{\link{as.pmob}}).
#' @param pcaanchor,pcacomponent,pcacomponentsingle Principal,Component Analysis
#' (PCA) parameters (see \code{\link{as.pmob}}).
#' @param circlecomponent,circlecomponentsingle circle computation parameters
#' (see \code{\link{as.pmob}}).
#' @param comments (see \code{\link{as.pmob}}).
#' @param originalfile,originalfileline information on the original date file
#' (see \code{\link{as.pmob}}).
#'
#' @examples
#' pmob <- as.pmob(sampid = c("Bosso1","Bosso1", "Bosso2", "Bosso3"),
#'                 specid = rep("b", 4))
#'
#' pmob.add(pmob = pmob, measyear = 2011, extra = 15,
#'          depth = c(20,20,30, 40))
#'
#' @export

pmob.add <- function(...,
                     pmob,
                     sampid       = NULL,
                     specid       = NULL,
                     slotid       = NULL,
                     measid       = NULL,
                     siteid       = NULL,
                     initialorder = NULL,
                     depth        = NULL,

                     meassec   = NULL,
                     measmin   = NULL,
                     meashour  = NULL,
                     measday   = NULL,
                     measmonth = NULL,
                     measyear  = NULL,

                     measdevice   = NULL,
                     magsusdevice = NULL,

                     xint = NULL,
                     yint = NULL,
                     zint = NULL,

                     blankcor    = NULL,
                     isblank     = NULL,
                     blankmethod = NULL,

                     holdercor = NULL,
                     isholder  = NULL,
                     holdermethod = NULL,

                     xinterror = NULL,
                     yinterror = NULL,
                     zinterror = NULL,

                     xvol = NULL,
                     yvol = NULL,
                     zvol = NULL,

                     xmass = NULL,
                     ymass = NULL,
                     zmass = NULL,

                     totmagsus = NULL,

                     totmagsuserror = NULL,

                     volmagsus  = NULL,
                     massmagsus = NULL,

                     vol  = NULL,
                     mass = NULL,
                     dens = NULL,

                     exactvol  = NULL,
                     exactmass = NULL,
                     exactdens = NULL,

                     discrete = NULL,

                     area = NULL,

                     defaultarea = NULL,

                     xefflength = NULL,
                     yefflength = NULL,
                     zefflength = NULL,

                     xeffvol = NULL,
                     yeffvol = NULL,
                     zeffvol = NULL,

                     xeffmass = NULL,
                     yeffmass = NULL,
                     zeffmass = NULL,

                     magsuseffvol  = NULL,
                     magsuseffmass = NULL,

                     deconvoluted = NULL,

                     xintconv = NULL,
                     yintconv = NULL,
                     zintconv = NULL,

                     sampaz  = NULL,
                     sampdip = NULL,
                     samprot = NULL,

                     coraz  = NULL,
                     cordip = NULL,
                     corrot = NULL,

                     bedaz  = NULL,
                     beddip = NULL,

                     foldtrend  = NULL,
                     foldplunge = NULL,

                     magaz = NULL,

                     usemagaz = NULL,

                     solaraz    = NULL,
                     solarangle = NULL,

                     long   = NULL,
                     lat    = NULL,
                     height = NULL,

                     sampmin   = NULL,
                     samphour  = NULL,
                     sampday   = NULL,
                     sampmonth = NULL,
                     sampyear  = NULL,

                     samptimezonemin  = NULL,
                     samptimezonehour = NULL,

                     magvar = NULL,

                     magazvarcor     = NULL,
                     bedazvarcor     = NULL,
                     foldtrendcarcor = NULL,

                     treatafx = NULL,
                     treatafy = NULL,
                     treatafz = NULL,

                     treattemp = NULL,

                     treatirmx = NULL,
                     treatirmy = NULL,
                     treatirmz = NULL,

                     treatarmafx = NULL,
                     treatarmafy = NULL,
                     treatarmafz = NULL,

                     treatarmbiasx = NULL,
                     treatarmbiasy = NULL,
                     treatarmbiasz = NULL,

                     pcaanchor             = NULL,
                     pcacomponent          = NULL,
                     pcacomponentsingle    = NULL,
                     circlecomponent       = NULL,
                     circlecomponentsingle = NULL,

                     comments = NULL,
                     originalfile = NULL,
                     originalfileline = NULL)
{
  if(!is.pmob(pmob)) stop("The 'pmob' parameter should be a pmob")

  initialname <-deparse(substitute(pmob))

  pmob.addon <- as.pmob(sampid        = sampid,
                        specid        = specid,
                        slotid        = slotid,
                        measid        = measid,
                        siteid        = siteid,
                        initialorder  = initialorder,
                        depth         = depth,

                        meassec   = meassec,
                        measmin   = measmin,
                        meashour  = meashour,
                        measday   = measday,
                        measmonth = measmonth,
                        measyear  = measyear,

                        measdevice   = measdevice,
                        magsusdevice = magsusdevice,

                        xint = xint,
                        yint = yint,
                        zint = zint,

                        blankcor    = blankcor,
                        isblank     = isblank,
                        blankmethod = blankmethod,

                        holdercor    = holdercor,
                        isholder     = isholder,
                        holdermethod = holdermethod,

                        xinterror = xinterror,
                        yinterror = yinterror,
                        zinterror = zinterror,

                        xvol = xvol,
                        yvol = yvol,
                        zvol = zvol,

                        xmass = xmass,
                        ymass = ymass,
                        zmass = zmass,

                        totmagsus = totmagsus,

                        totmagsuserror = totmagsuserror,

                        volmagsus  = volmagsus,
                        massmagsus = massmagsus,

                        vol  = vol,
                        mass = mass,
                        dens = dens,

                        exactvol  = exactvol,
                        exactmass = exactmass,
                        exactdens = exactdens,

                        discrete = discrete,

                        area = area,

                        defaultarea = defaultarea,

                        xefflength = xefflength,
                        yefflength = yefflength,
                        zefflength = zefflength,

                        xeffvol = xeffvol,
                        yeffvol = yeffvol,
                        zeffvol = zeffvol,

                        xeffmass = xeffmass,
                        yeffmass = yeffmass,
                        zeffmass = zeffmass,

                        magsuseffvol  = magsuseffvol,
                        magsuseffmass = magsuseffmass,

                        deconvoluted = deconvoluted,

                        xintconv = xintconv,
                        yintconv = yintconv,
                        zintconv = zintconv,

                        sampaz  = sampaz,
                        sampdip = sampdip,
                        samprot = samprot,

                        coraz  = coraz,
                        cordip = cordip,
                        corrot = corrot,

                        bedaz  = bedaz,
                        beddip = beddip,

                        foldtrend  = foldtrend,
                        foldplunge = foldplunge,

                        magaz = magaz,

                        usemagaz = usemagaz,

                        solaraz    = solaraz,
                        solarangle = solarangle,

                        long   = long,
                        lat    = lat,
                        height = height,

                        sampmin   = sampmin,
                        samphour  = samphour,
                        sampday   = sampday,
                        sampmonth = sampmonth,
                        sampyear  = sampyear,

                        samptimezonemin  = samptimezonemin,
                        samptimezonehour = samptimezonehour,

                        magvar = magvar,

                        magazvarcor     = magazvarcor,
                        bedazvarcor     = bedazvarcor,
                        foldtrendcarcor = foldtrendcarcor,

                        treatafx = treatafx,
                        treatafy = treatafy,
                        treatafz = treatafz,

                        treattemp = treattemp,

                        treatirmx = treatirmx,
                        treatirmy = treatirmy,
                        treatirmz = treatirmz,

                        treatarmafx = treatarmafx,
                        treatarmafy = treatarmafy,
                        treatarmafz = treatarmafz,

                        treatarmbiasx = treatarmbiasx,
                        treatarmbiasy = treatarmbiasy,
                        treatarmbiasz = treatarmbiasz,

                        pcaanchor             = pcaanchor,
                        pcacomponent          = pcacomponent,
                        pcacomponentsingle    = pcacomponentsingle,
                        circlecomponent       = circlecomponent,
                        circlecomponentsingle = circlecomponentsingle,

                        comments = comments,
                        originalfile = originalfile,
                        originalfileline = originalfileline,
                        ...
  )

  nr1 <- nrow(pmob)
  nr2 <- nrow(pmob.addon)

  if(nr2 == 0){
    return(pmob)
  }

  if(nr2 == 1){
    pmob.addon <- pmob.addon[rep(1, nr1),]
  } else if(nr2 != nr1){
    stop("The added parameters should have length 1, ",
         "or have as many elements than pmob has columns (", nr1, ").")
  }

  remove <- !(names(pmob.addon) %in% c("initialorder", "version"))

  pmob.addon <- pmob.addon[, remove, drop = F]

  dupl <- names(pmob.addon) %in% names(pmob)

  if(any(dupl)){
    warning("The following column(s) is (are) already found in ", initialname,
            ":\n - ",
            paste(names(pmob.addon)[dupl], collapse = "\n - "),
            "\nThey will not be overwriten in the resulting pmob",
            "\nTo overwrite the pmob, you can remove them\nbeforehand by running",
            " the following lines of code:\n\n",
            paste(initialname, "$", names(pmob.addon)[dupl],
                  " <- NULL\n", sep = ""))

    pmob.addon <- pmob.addon[,!dupl, drop = F]

    if(ncol(pmob.addon) == 0) return(pmob)
  }

  res <- cbind(pmob, pmob.addon)

  res <- pmob.col.sort(res)

  row.names(res) <- NULL

  return(res)

}
