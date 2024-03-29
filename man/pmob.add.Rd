% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmob.add.R
\name{pmob.add}
\alias{pmob.add}
\title{Add new columns to a pmob}
\usage{
pmob.add(
  ...,
  pmob,
  sampid = NULL,
  specid = NULL,
  slotid = NULL,
  measid = NULL,
  siteid = NULL,
  initialorder = NULL,
  depth = NULL,
  meassec = NULL,
  measmin = NULL,
  meashour = NULL,
  measday = NULL,
  measmonth = NULL,
  measyear = NULL,
  measdevice = NULL,
  magsusdevice = NULL,
  xint = NULL,
  yint = NULL,
  zint = NULL,
  blankcor = NULL,
  isblank = NULL,
  blankmethod = NULL,
  holdercor = NULL,
  isholder = NULL,
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
  volmagsus = NULL,
  massmagsus = NULL,
  vol = NULL,
  mass = NULL,
  dens = NULL,
  exactvol = NULL,
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
  magsuseffvol = NULL,
  magsuseffmass = NULL,
  deconvoluted = NULL,
  xintconv = NULL,
  yintconv = NULL,
  zintconv = NULL,
  sampaz = NULL,
  sampdip = NULL,
  samprot = NULL,
  coraz = NULL,
  cordip = NULL,
  corrot = NULL,
  bedaz = NULL,
  beddip = NULL,
  foldtrend = NULL,
  foldplunge = NULL,
  magaz = NULL,
  usemagaz = NULL,
  solaraz = NULL,
  solarangle = NULL,
  long = NULL,
  lat = NULL,
  height = NULL,
  sampmin = NULL,
  samphour = NULL,
  sampday = NULL,
  sampmonth = NULL,
  sampyear = NULL,
  samptimezonemin = NULL,
  samptimezonehour = NULL,
  magvar = NULL,
  magazvarcor = NULL,
  bedazvarcor = NULL,
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
  pcaanchor = NULL,
  pcacomponent = NULL,
  pcacomponentsingle = NULL,
  circlecomponent = NULL,
  circlecomponentsingle = NULL,
  comments = NULL,
  originalfile = NULL,
  originalfileline = NULL
)
}
\arguments{
\item{...}{Personalized parameters (see \code{\link{as.pmob}}).}

\item{pmob}{pmob object to complement}

\item{sampid, specid, slotid, measid, siteid}{Identification parameters
(see \code{\link{as.pmob}}).}

\item{initialorder}{Initial order (see \code{\link{as.pmob}}).}

\item{depth}{Sampling depth (see \code{\link{as.pmob}}).}

\item{meassec, measmin, meashour, measday, measmonth, measyear}{Time of the measurement (see \code{\link{as.pmob}}).}

\item{measdevice, magsusdevice}{Measurement devices
(see \code{\link{as.pmob}}).}

\item{xint, yint, zint}{Magnetic moments (see \code{\link{as.pmob}}).}

\item{blankcor, isblank, blankmethod}{Correction parameters (see \code{\link{as.pmob}}).}

\item{holdercor, isholder, holdermethod}{Correction parameters (see \code{\link{as.pmob}}).}

\item{xinterror, yinterror, zinterror}{Error (2σ) on the magnetic moment
(see \code{\link{as.pmob}}).}

\item{xvol, yvol, zvol}{Magnetization by volume (see \code{\link{as.pmob}}).}

\item{xmass, ymass, zmass}{Magnetization by mass (see \code{\link{as.pmob}}).}

\item{totmagsus, totmagsuserror, volmagsus, massmagsus}{Magnetic susceptibility
data (see \code{\link{as.pmob}}).}

\item{vol, mass, dens}{Sample properties (see \code{\link{as.pmob}}).}

\item{exactvol, exactmass, exactdens}{Sample properties reliability
(see \code{\link{as.pmob}}).}

\item{discrete}{Whether a single discrete sample is measured (TRUE), or a
continuous core (FALSE) (see \code{\link{as.pmob}}).}

\item{area}{Cross-sectional area of a long core sample (for continuous cores,
applies if the parameter \strong{discrete} is FALSE)
(see \code{\link{as.pmob}}).}

\item{defaultarea}{Long core sample properties reliability (see \code{\link{as.pmob}}).}

\item{xefflength, yefflength, zefflength}{Effective portion of length of a long
core that is measured by the sensor (see \code{\link{as.pmob}}).}

\item{xeffvol, yeffvol, zeffvol}{DEffective volume of a long core that is
measured by the sensors  (see \code{\link{as.pmob}}).}

\item{xeffmass, yeffmass, zeffmass}{Effective mass of a long core that is
measured by the sensors (see \code{\link{as.pmob}}).}

\item{magsuseffvol, magsuseffmass}{Magnetic susceptibilities of a long core
(see \code{\link{as.pmob}}).}

\item{deconvoluted}{Whether the continuous core measurements have been
deconvoluted (TRUE) or not (FALSE) (see \code{\link{as.pmob}}).}

\item{xintconv, yintconv, zintconv}{The original magnetic moment obtained
without deconvolution (see \code{\link{as.pmob}}).}

\item{sampaz, sampdip, samprot}{Sample orientation in the field
(see \code{\link{as.pmob}}).}

\item{coraz, cordip, corrot}{Sample orientation in the
measuring device (see \code{\link{as.pmob}}).}

\item{bedaz, beddip}{Orientation of beds (see \code{\link{as.pmob}}).}

\item{foldtrend, foldplunge}{Orientation of the folding
(see \code{\link{as.pmob}}).}

\item{magaz, usemagaz, solaraz}{Parameters determining the zimuth of the sample
(see \code{\link{as.pmob}}).}

\item{solarangle}{Solar angle (see \code{\link{as.pmob}}).}

\item{long, lat, height}{Coordinates (see \code{\link{as.pmob}}).}

\item{sampmin, samphour, sampday, sampmonth, sampyear}{Time of sampling  (see \code{\link{as.pmob}}).}

\item{samptimezonemin, samptimezonehour}{Time zone of the sampling
(see \code{\link{as.pmob}}).}

\item{magvar, magazvarcor, bedazvarcor, foldtrendcarcor}{Parameters to
correct due to the magnetic variation (see \code{\link{as.pmob}}).}

\item{treatafx, treatafy, treatafz}{Treatment field by alternating field (AF)
demagnetization (see \code{\link{as.pmob}}).}

\item{treattemp}{Treatment temperature (see \code{\link{as.pmob}}).}

\item{treatirmx, treatirmy, treatirmz}{Treatment field by Isothermal Remanent
Magnetisation (IRM) (see \code{\link{as.pmob}}).}

\item{treatarmafx, treatarmafy, treatarmafz}{Treatment anhysteretic field
for Anhysteretic Remanent Magnetization (ARM) (see \code{\link{as.pmob}}).}

\item{treatarmbiasx, treatarmbiasy, treatarmbiasz}{Treatment bias field for
Anhysteretic Remanent Magnetization (ARM) (see \code{\link{as.pmob}}).}

\item{pcaanchor, pcacomponent, pcacomponentsingle}{Principal,Component Analysis
(PCA) parameters (see \code{\link{as.pmob}}).}

\item{circlecomponent, circlecomponentsingle}{circle computation parameters
(see \code{\link{as.pmob}}).}

\item{comments}{(see \code{\link{as.pmob}}).}

\item{originalfile, originalfileline}{information on the original date file
(see \code{\link{as.pmob}}).}
}
\description{
Add new columns to a pmob (\strong{P}aleo\strong{M}agnetic
\strong{OB}jects).
}
\examples{
pmob <- as.pmob(sampid = c("Bosso1","Bosso1", "Bosso2", "Bosso3"),
                specid = rep("b", 4))

pmob.add(pmob = pmob, measyear = 2011, extra = 15,
         depth = c(20,20,30, 40))

}
