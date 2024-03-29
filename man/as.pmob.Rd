% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.pmob.R
\name{as.pmob}
\alias{as.pmob}
\title{Generate a pmob}
\usage{
as.pmob(
  ...,
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
\item{...}{personalized parameters, which can be of length 1 or n (length of
the other parameters), and of class 'numeric', 'integer', 'logical',
'character' or 'factor'. Header name is the name of the argument (e.g.
\code{as.pmob(sampleid = "Bosso1", nameofcolumn = "information")}). These can
also be provided as a list.}

\item{sampid}{Identification of the sample [alphanumeric].}

\item{specid}{Identification of the specimen  [alphanumeric].}

\item{slotid}{Identification of the tray slot containing the sample for
discrete measurements on a multi-sample tray [alphanumeric].}

\item{measid}{Identification of the measurement [alphanumeric].}

\item{siteid}{Identification of the sampling site [alphanumeric].}

\item{initialorder}{Position or order of the measurements in the initial file
[positive non null integer].}

\item{depth}{Sampling depth [normal numeric; in m].}

\item{meassec}{Time of the measurement [normal numeric in [0,60[; in
seconds].}

\item{measmin}{Time of the measurement [integer in [0,59]; in minutes].}

\item{meashour}{Time of the measurement [integer in [0,23]; in hours].}

\item{measday}{Time of the measurement [integer in [1,31]; in days].}

\item{measmonth}{Time of the measurement [integer in [0,12]; in months].}

\item{measyear}{Time of the measurement [integer; in years].}

\item{measdevice}{Name of the measurement device [alphanumeric].}

\item{magsusdevice}{Name of the measurement device for magnetic
susceptibility [alphanumeric].}

\item{xint}{Magnetic moment in the x direction of the measuring device
[scientific numeric; in \out{A m<sup>2</sup>}].}

\item{yint}{Magnetic moment in the y direction of the measuring device
[scientific numeric; in \out{A m<sup>2</sup>}].}

\item{zint}{Magnetic moment in the z direction of the measuring device
[scientific numeric; in \out{A m<sup>2</sup>}].}

\item{blankcor}{Whether the magnetic moment has been corrected for the
background level [TRUE/FALSE].}

\item{isblank}{Whether this is a measurement of the background level of
magnetic moment (serves for automated correction) [TRUE/FALSE].}

\item{blankmethod}{The methodology used for computing the blank value at the
time of the sample measurement for correction. Supported values are "approx"
for a linear approximation between the previous and the next blank values,
using the time of the measurement, "mean" stands for a mean of the previous
and the next blank values, "previous" takes the previous blank value, and
"next", the next one [alphanumeric].}

\item{holdercor}{Whether the magnetic moment has been corrected for the
holder [TRUE/FALSE].}

\item{isholder}{Whether this is a measurement of the magnetic moment of the
empty holder (serves for automated correction) [TRUE/FALSE].}

\item{holdermethod}{The methodology used for computing the holder value at
the time of the sample measurement for correction. Supported values are
"approx" for a linear approximation between the previous and the next holder
values, using the time of the measurement, "mean" stands for a mean of the
previous and the next holder values, "previous" takes the previous holder
value, and "next", the next one [alphanumeric].}

\item{xinterror}{Error (2σ) on the magnetic moment in the
x direction of the measuring device
[scientific numeric; in \out{A m<sup>2</sup>}].}

\item{yinterror}{Error (2σ) on the magnetic moment in the
y direction of the measuring device
[scientific numeric; in \out{A m<sup>2</sup>}].}

\item{zinterror}{Error (2σ) on the magnetic moment in the
z direction of the measuring device
[scientific numeric; in \out{A m<sup>2</sup>}].}

\item{xvol}{Magnetization by volume in the x direction of the measuring
device [scientific numeric; in \out{A m<sup>-1</sup>}].}

\item{yvol}{Magnetization by volume in the y direction of the measuring
device [scientific numeric; in \out{A m<sup>-1</sup>}].}

\item{zvol}{Magnetization by volume in the z direction of the measuring
device [scientific numeric; in \out{A m<sup>-1</sup>}].}

\item{xmass}{Magnetization by mass in the x direction of the measuring device
[scientific numeric; in \out{A m<sup>2</sup> kg<sup>-1</sup>}].}

\item{ymass}{Magnetization by mass in the y direction of the measuring device
[scientific numeric; in \out{A m<sup>2</sup> kg<sup>-1</sup>}].}

\item{zmass}{Magnetization by mass in the z direction of the measuring device
[scientific numeric; in \out{A m<sup>2</sup> kg<sup>-1</sup>}].}

\item{totmagsus}{Total magnetic susceptibility, which is defined by the
magnetic moment m of a sample (expressed in \out{A m<sup>2</sup>}) induced
by an external field H (expressed in \out{A m<sup>-1</sup>}). The total
magnetic susceptibility is m/H (expressed in \out{m<sup>3</sup>})
 [scientific numeric; in \out{m<sup>3</sup>}]. See also Tauxe, 2010.}

\item{totmagsuserror}{Error (2σ) on the total magnetic susceptibility
[scientific numeric; in \out{m<sup>3</sup>}].}

\item{volmagsus}{Magnetic susceptibility by volume [scientific numeric; in SI
(dimensionless)].}

\item{massmagsus}{Magnetic susceptibility by mass [scientific numeric; in
\out{m<sup>3</sup> kg<sup>-1</sup>}].}

\item{vol}{Volume of the sample. For long cores, it is used for the entire
core sample to compute density [positive scientific numeric; in
\out{m<sup>3</sup>}].}

\item{mass}{Mass of the sample. For long cores, it is used for the entire
sample to compute density [positive scientific numeric; in kg].}

\item{dens}{Density of the sample. In discrete samples it is used to
compute volume from mass of vice-versa. In continuous samples it is used to
compute effective mass from the effective volume [positive scientific
number; in \out{kg m<sup>-3</sup>}].}

\item{exactvol}{Whether the volume is only an approximation (FALSE) or is
truly measured (TRUE) [TRUE/FALSE].}

\item{exactmass}{Whether the mass is only an approximation (FALSE) or is
truly measured (TRUE) [TRUE/FALSE].}

\item{exactdens}{Whether the density is only an approximation (FALSE) or
is truly measured (TRUE) [TRUE/FALSE].}

\item{discrete}{Whether a single discrete sample is measured (TRUE), or a
continuous core (FALSE) [TRUE/FALSE].}

\item{area}{Cross-sectional area of a long core sample (for continuous cores,
applies if the parameter \strong{discrete} is FALSE) [positive scientific
numeric; in \out{m<sup>2</sup>}].}

\item{defaultarea}{Whether the area is only a default value (TRUE) or the
actual standard cross-sectional area of the core (FALSE). We mention that
the actual cross-sectional area can vary throughout the core, e.g., if the
core is deteriorated: this is not taken into account by this parameter
[TRUE/FALSE].}

\item{xefflength}{Effective portion of length of a long core that is measured
by the x sensor. This is multiplied by \strong{area} to determine the
effective volume of the sample measured for x [positive scientific numeric;
in m].}

\item{yefflength}{Effective portion of length of a long core that is measured
by the y sensor. This is multiplied by \strong{area} to determine the
effective volume of the sample measured for y [positive scientific numeric;
in m].}

\item{zefflength}{Effective portion of length of a long core that is measured
by the z sensor. This is multiplied by \strong{area} to determine the
effective volume of the sample measured for z [positive scientific numeric;
in m].}

\item{xeffvol}{Effective volume of a long core that is measured by the x
sensor [positive scientific numeric; in \out{m<sup>3</sup>}]}

\item{yeffvol}{Effective volume of a long core that is measured by the y
sensor [positive scientific numeric; in \out{m<sup>3</sup>}]}

\item{zeffvol}{DEffective volume of a long core that is measured by the z
sensor [positive scientific numeric; in \out{m<sup>3</sup>}]}

\item{xeffmass}{Effective mass of a long core that is measured by the x
sensor [positive scientific numeric; in kg].}

\item{yeffmass}{Effective mass of a long core that is measured by the y
sensor [positive scientific numeric; in kg].}

\item{zeffmass}{Effective mass of a long core that is measured by the z
sensor [positive scientific numeric; in kg].}

\item{magsuseffvol}{Effective volume of a long core that is measured for
magnetic susceptibility [positive scientific number; in \out{m<sup>3</sup>}].}

\item{magsuseffmass}{Effective mass of a long core that is measured for
magnetic susceptibility [positive scientific number; in kg].}

\item{deconvoluted}{Whether the continuous core measurements have been
deconvoluted (TRUE) or not (FALSE) [TRUE/FALSE].}

\item{xintconv}{The original magnetic moment in the x direction of the
measuring device, obtained without deconvolution [scientific number; in
\out{A m<sup>2</sup>}].}

\item{yintconv}{The original magnetic moment in the y direction of the
measuring device, obtained without deconvolution [scientific number; in
\out{A m<sup>2</sup>}].}

\item{zintconv}{The original magnetic moment in the z direction of the
measuring device, obtained without deconvolution [scientific number; in
\out{A m<sup>2</sup>}].}

\item{sampaz}{Sample azimuth in the field (see Fig. 2 in supporting
material): it is the angle measured clockwise from the north of the
horizontal projection of the field arrow (Tauxe 2010) [normal numeric in
[0,360[; in arc degree].}

\item{sampdip}{Sample dip in the field (see Fig. 2 in supporting material):
it is the angle of the  field arrow from the horizontal. It is positive
downward, and ranges from +90° for straight down to -90° for straight up
(Tauxe, 2010) [normal numeric in [-90,90]; in arc degree].}

\item{samprot}{Rotation of the sample on its axis, taken in the field (see
Fig. 2 in supporting material). It is measured clockwise from the 12 o’clock
summit (or X sample coordinate, see Fig. 2B in supporting material), on the
sample part opposite of the field arrow. \strong{SPECIAL CASE: IN DOWNWARD
VERTICAL SAMPLES} (dip = 90), the rotation is the angle between the azimuth
and the field arrow.
\strong{SPECIAL CASE: IN UPWARD VERTICAL SAMPLES} (dip = -90), the rotation
is the angle between the azimuth and the field arrow + 180° (as the upward
dip brings the 12 o’clock position (or X sample coordinate, see Fig. 2B in
supporting material) on the sample part opposite of the field arrow, to face
the opposite direction of azimuth). [normal numeric in [0,360[; in arc
degree].}

\item{coraz}{Azimuth of the sample in the measuring device (see
Fig. 3, 4 & 5 in supplementary material): this is a correction for the
difference between the magnetometer coordinates and the sample coordinates.
This is equivalent to sample azimuth, considering that the x, y and z
coordinates of the magnetometer are attributed to the North, East and
downward directions respectively [normal numeric in [0,360[; in arc degree].}

\item{cordip}{Dip of the sample in the measuring device (see
Fig. 3, 4 & 5 in supplementary material):  this is a correction for the
difference between the magnetometer coordinates and the sample coordinates.
This is equivalent to sample dip, considering that the x, y and z coordinates
of the magnetometer are attributed to the North, East and downward directions
respectively [normal numeric in [-90,90]; in arc degree].}

\item{corrot}{Rotation of the sample in the measuring device (see
Fig. 3, 4 & 5 in supplementary material):  this is a correction for the
difference between the magnetometer coordinates and the sample coordinates.
This is equivalent to sample rotation, considering that the x, y and z
coordinates of the magnetometer are attributed to the North, East and
downward directions respectively [normal numeric in [0,360[; in arc degree].}

\item{bedaz}{Bedding azimuth, or dip direction (see Fig. 1 in supplementary
material): it is the azimuth (the angle taken eastward or clockwise from the
north) of the line perpendicular to the plane oriented in the stratigraphic
upward direction, projected on a horizontal plane. This accounts for
stratigraphically overturned beds (see Fig. 1d-f in supplementary
material) [normal numeric in [0,360[; in arc degree].}

\item{beddip}{Bedding dip (see Fig. 1 in supplementary material): it is the
plane's maximum angular deviation from the horizontal. It is positive
downward, and ranges from +90° for straight down to 0° for horizontal.
Stratigraphically overturned beds are indicated with dip values in the
]90°,180°] interval [normal numeric in [0,180]; in arc degree].}

\item{foldtrend}{Fold trend (see Fig. 1 in supplementary material): it is the
trend of the axis of a fold (or its azimuth, i.e. the angle taken eastward or
clockwise from the north) [normal numeric in [0,360[; in arc degree].}

\item{foldplunge}{Fold plunge (see Fig. 1 in supplementary material): it is
the angle that the axis of a fold makes with the horizontal. It is positive
downwards, and ranges from +90° form straight down to 0° for horizontal.
Stratigraphically overturned folds are indicated with plunge values in the
]90°,180°] interval [normal numeric in [0,180]; in arc degree].}

\item{magaz}{Magnetic azimuth of the sample in the field [numeric in [0,360[;
in arc degree].}

\item{usemagaz}{If TRUE, the magnetic azimuth (\strong{magaz}) will be used
as the sample azimuth (\strong{sampaz}), if FALSE it is the solar azimuth
(\strong{solaraz}) that will be used as the sample azimuth
(\strong{sampaz}).}

\item{solaraz}{Solar azimuth of the sample in the field [numeric in [0,360[;
in arc degree].}

\item{solarangle}{Angle that the shadow of a vertical needle does with the
core’s field arrow at the moment of sampling (see the \strong{sampmin},
\strong{samphour}, \strong{sampday}, \strong{sampmonth} and
\strong{sampyear} parameters): this is used to compute the solar azimuth
(\strong{solaraz}) [numeric in [0,360[; in arc degree].}

\item{long}{Longitude [decimal degree in ]-180,180]; in decimal degree].}

\item{lat}{Latitude [decimal degree in [-90,90]; in decimal degree].}

\item{height}{Height [numeric; in m].}

\item{sampmin}{Time of sampling [normal numeric in [0,60[; in minutes].}

\item{samphour}{Time of sampling [integer in [0,23]; in hours].}

\item{sampday}{Time of sampling [integer in [1,31]; in days].}

\item{sampmonth}{Time of sampling [integer in [0,12]; in months].}

\item{sampyear}{Time of sampling [integer; in years].}

\item{samptimezonemin}{Time zone of the sampling, expressed relative to
the UTC (Coordinated Universal Time), e.g. UTC +12:45 for New Zealand. This
has to take into account the possible daylight saving time at the moment of
the sampling, which changes the effective time zone [Integer in [0,59]; in
minutes].}

\item{samptimezonehour}{Time zone of the sampling, expressed relative to
the UTC (Coordinated Universal Time), e.g. UTC +12:45 for New Zealand. This
has to take into account the possible daylight saving time at the moment of
the sampling, which changes the effective time zone [integer, usually in
[-12,14]; in hours].}

\item{magvar}{Magnetic variation: it is the angle on the horizontal plane
between the magnetic north and the geographic north. It is measured as the
angle made by the magnetic North eastward (clockwise) from the geographic
North [numeric in [0,360[; in arc degrees].}

\item{magazvarcor}{Whether the magnetic azimuth parameter (\strong{magaz})
has to be corrected for the magnetic variation [TRUE/FALSE].}

\item{bedazvarcor}{Whether the bedding azimuth parameter (\strong{bedaz})
has to be corrected for the magnetic variation. [TRUE/FALSE].}

\item{foldtrendcarcor}{Whether the fold trend parameter (\strong{foldtrend})
has to be corrected for the magnetic variation [TRUE/FALSE].}

\item{treatafx}{Treatment field by alternating field (AF) demagnetization in
the x direction of the measuring device [normal numeric; in T].}

\item{treatafy}{Treatment field by alternating field (AF) demagnetization in
the y direction of the measuring device [normal numeric; in T].}

\item{treatafz}{Treatment field by alternating field (AF) demagnetization in
the z direction of the measuring device [normal numeric; in T].}

\item{treattemp}{Treatment temperature in Kelvin [positive normal numeric;
in K].}

\item{treatirmx}{Treatment field by Isothermal Remanent Magnetisation (IRM)
in the x direction of the measuring device [normal numeric; in T].}

\item{treatirmy}{Treatment field by Isothermal Remanent Magnetisation (IRM)
in the y direction of the measuring device [normal numeric; in T].}

\item{treatirmz}{Treatment field by Isothermal Remanent Magnetisation (IRM)
in the z direction of the measuring device [normal numeric; in T].}

\item{treatarmafx}{Treatment anhysteretic field in the x direction of the
measuring device for Anhysteretic Remanent Magnetization (ARM) [normal
numeric; in T].}

\item{treatarmafy}{Treatment anhysteretic field in the y direction of the
measuring device for Anhysteretic Remanent Magnetization (ARM) [normal
numeric; in T].}

\item{treatarmafz}{Treatment anhysteretic field in the z direction of the
measuring device for Anhysteretic Remanent Magnetization (ARM) [normal
numeric; in T].}

\item{treatarmbiasx}{Treatment bias field in x direction of the measuring
device for Anhysteretic Remanent Magnetization (ARM). In most settings the
bias field comes from the z axis of the measuring device, this parameter is
only set to allow for unconventional settings [normal numeric; in T].}

\item{treatarmbiasy}{Treatment bias field in y direction of the measuring
device for Anhysteretic Remanent Magnetization (ARM). In most settings the
bias field comes from the z axis of the measuring device, this parameter is
only set to allow for unconventional settings [normal numeric; in T].}

\item{treatarmbiasz}{Treatment bias field in z direction for Anhysteretic
Remanent Magnetization (ARM) [normal numeric; in T].}

\item{pcaanchor}{Whether to anchor the axis of a Principal Component Analysis
(PCA) at the origin (x = 0, y = 0, z = 0) [TRUE/FALSE].}

\item{pcacomponent}{Identifies groups of measurements considered in one
component (e.g. “C1”, “C2”, “Outlier”), to be used for Principal Component
Analysis (PCA) [alphanumeric].}

\item{pcacomponentsingle}{Which group of measurements to use for Principal
Component Analysis (PCA) in files only considering one single group
[alphanumeric].}

\item{circlecomponent}{Identifies groups of measurements considered in one
component (e.g. “C1”, “C2”, “Outlier”), to be used for circle computation
[alphanumeric].}

\item{circlecomponentsingle}{Which group of measurements to use for circle
computation in files only considering one single group. [alphanumeric].}

\item{comments}{Additional comments [alphanumeric].}

\item{originalfile}{The name of the original file from which the data is from
[alphanumeric].}

\item{originalfileline}{The line, in the original file, of the sample
identification information (used for \strong{sampid}) related to the x
magnetic moment data. If not available, the line of the sample identification
related to the x magnetisation by volume or by mass (in this order). This is
to improve the manual or automated localisation of mistakes in the original
file, based on automated or manual checks of the pmob.}
}
\description{
Function to generate a pmob (\strong{P}aleo\strong{M}agnetic
\strong{OB}ject); you could do it using the \code{data.frame} function, but
this function formats everything neatly and includes a few useful checks.
}
\examples{
as.pmob(exc = 1, sampid = c("Bosso1", "MonteAcuto1"))

}
\references{
\itemize{
  \item Tauxe, L., 2010. Essentials of Paleomagnetism. University of
  California Press.
}
}
