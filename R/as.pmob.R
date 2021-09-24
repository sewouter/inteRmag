#' @title Generate a pmob
#'
#' @description Function to generate a pmob (\strong{P}aleo\strong{M}agnetic
#' \strong{OB}ject); you could do it using the \code{data.frame} function, but
#' this function formats everything neatly and includes a few useful checks.
#'
#' @param ... personalized parameters, which can be of length 1 or n (length of
#' the other parameters), and of class 'numeric', 'integer', 'logical',
#' 'character' or 'factor'. Header name is the name of the argument (e.g.
#' \code{as.pmob(sampleid = "Bosso1", nameofcolumn = "information")}). These can
#' also be provided as a list.
#' @param sampid Identification of the sample [alphanumeric].
#' @param specid  Identification of the specimen  [alphanumeric].
#' @param slotid Identification of the tray slot containing the sample for
#' discrete measurements on a multi-sample tray [alphanumeric].
#' @param measid Identification of the measurement [alphanumeric].
#' @param siteid Identification of the sampling site [alphanumeric].
#' @param initialorder Position or order of the measurements in the initial file
#' [positive non null integer].
#' @param depth Sampling depth [normal numeric; in m].
#' @param meassec Time of the measurement [normal numeric in [0,60[; in
#' seconds].
#' @param measmin Time of the measurement [integer in [0,59]; in minutes].
#' @param meashour Time of the measurement [integer in [0,23]; in hours].
#' @param measday Time of the measurement [integer in [1,31]; in days].
#' @param measmonth Time of the measurement [integer in [0,12]; in months].
#' @param measyear Time of the measurement [integer; in years].
#' @param measdevice Name of the measurement device [alphanumeric].
#' @param magsusdevice Name of the measurement device for magnetic
#' susceptibility [alphanumeric].
#' @param xint Magnetic moment in the x direction of the measuring device
#' [scientific numeric; in \out{A m<sup>2</sup>}].
#' @param yint Magnetic moment in the y direction of the measuring device
#' [scientific numeric; in \out{A m<sup>2</sup>}].
#' @param zint Magnetic moment in the z direction of the measuring device
#' [scientific numeric; in \out{A m<sup>2</sup>}].
#' @param blankcor Whether the magnetic moment has been corrected for the
#' background level [TRUE/FALSE].
#' @param isblank Whether this is a measurement of the background level of
#' magnetic moment (serves for automated correction) [TRUE/FALSE].
#' @param blankmethod The methodology used for computing the blank value at the
#' time of the sample measurement for correction. Supported values are "approx"
#' for a linear approximation between the previous and the next blank values,
#' using the time of the measurement, "mean" stands for a mean of the previous
#' and the next blank values, "previous" takes the previous blank value, and
#' "next", the next one [alphanumeric].
#' @param holdercor Whether the magnetic moment has been corrected for the
#' holder [TRUE/FALSE].
#' @param isholder Whether this is a measurement of the magnetic moment of the
#' empty holder (serves for automated correction) [TRUE/FALSE].
#' @param holdermethod The methodology used for computing the holder value at
#' the time of the sample measurement for correction. Supported values are
#' "approx" for a linear approximation between the previous and the next holder
#' values, using the time of the measurement, "mean" stands for a mean of the
#' previous and the next holder values, "previous" takes the previous holder
#' value, and "next", the next one [alphanumeric].
#' @param xinterror Error (2σ) on the magnetic moment in the
#' x direction of the measuring device
#' [scientific numeric; in \out{A m<sup>2</sup>}].
#' @param yinterror Error (2σ) on the magnetic moment in the
#' y direction of the measuring device
#' [scientific numeric; in \out{A m<sup>2</sup>}].
#' @param zinterror Error (2σ) on the magnetic moment in the
#' z direction of the measuring device
#' [scientific numeric; in \out{A m<sup>2</sup>}].
#' @param xvol Magnetization by volume in the x direction of the measuring
#' device [scientific numeric; in \out{A m<sup>-1</sup>}].
#' @param yvol Magnetization by volume in the y direction of the measuring
#' device [scientific numeric; in \out{A m<sup>-1</sup>}].
#' @param zvol Magnetization by volume in the z direction of the measuring
#' device [scientific numeric; in \out{A m<sup>-1</sup>}].
#' @param xmass Magnetization by mass in the x direction of the measuring device
#' [scientific numeric; in \out{A m<sup>2</sup> kg<sup>-1</sup>}].
#' @param ymass Magnetization by mass in the y direction of the measuring device
#' [scientific numeric; in \out{A m<sup>2</sup> kg<sup>-1</sup>}].
#' @param zmass Magnetization by mass in the z direction of the measuring device
#' [scientific numeric; in \out{A m<sup>2</sup> kg<sup>-1</sup>}].
#' @param totmagsus Total magnetic susceptibility, which is defined by the
#' magnetic moment m of a sample (expressed in \out{A m<sup>2</sup>}) induced
#' by an external field H (expressed in \out{A m<sup>-1</sup>}). The total
#' magnetic susceptibility is m/H (expressed in \out{m<sup>3</sup>})
#'  [scientific numeric; in \out{m<sup>3</sup>}]. See also Tauxe, 2010.
#' @param totmagsuserror Error (2σ) on the total magnetic susceptibility
#' [scientific numeric; in \out{m<sup>3</sup>}].
#' @param volmagsus Magnetic susceptibility by volume [scientific numeric; in SI
#' (dimensionless)].
#' @param massmagsus Magnetic susceptibility by mass [scientific numeric; in
#' \out{m<sup>3</sup> kg<sup>-1</sup>}].
#' @param vol Volume of the sample. For long cores, it is used for the entire
#' core sample to compute density [positive scientific numeric; in
#' \out{m<sup>3</sup>}].
#' @param mass Mass of the sample. For long cores, it is used for the entire
#' sample to compute density [positive scientific numeric; in kg].
#' @param dens Density of the sample. In discrete samples it is used to
#' compute volume from mass of vice-versa. In continuous samples it is used to
#' compute effective mass from the effective volume [positive scientific
#' number; in \out{kg m<sup>-3</sup>}].
#' @param exactvol Whether the volume is only an approximation (FALSE) or is
#' truly measured (TRUE) [TRUE/FALSE].
#' @param exactmass Whether the mass is only an approximation (FALSE) or is
#' truly measured (TRUE) [TRUE/FALSE].
#' @param exactdens Whether the density is only an approximation (FALSE) or
#' is truly measured (TRUE) [TRUE/FALSE].
#' @param discrete Whether a single discrete sample is measured (TRUE), or a
#' continuous core (FALSE) [TRUE/FALSE].
#' @param area Cross-sectional area of a long core sample (for continuous cores,
#' applies if the parameter \strong{discrete} is FALSE) [positive scientific
#' numeric; in \out{m<sup>2</sup>}].
#' @param defaultarea Whether the area is only a default value (TRUE) or the
#' actual standard cross-sectional area of the core (FALSE). We mention that
#' the actual cross-sectional area can vary throughout the core, e.g., if the
#' core is deteriorated: this is not taken into account by this parameter
#' [TRUE/FALSE].
#' @param xefflength Effective portion of length of a long core that is measured
#' by the x sensor. This is multiplied by \strong{area} to determine the
#' effective volume of the sample measured for x [positive scientific numeric;
#' in m].
#' @param yefflength Effective portion of length of a long core that is measured
#' by the y sensor. This is multiplied by \strong{area} to determine the
#' effective volume of the sample measured for y [positive scientific numeric;
#' in m].
#' @param zefflength Effective portion of length of a long core that is measured
#' by the z sensor. This is multiplied by \strong{area} to determine the
#' effective volume of the sample measured for z [positive scientific numeric;
#' in m].
#' @param xeffvol Effective volume of a long core that is measured by the x
#' sensor [positive scientific numeric; in \out{m<sup>3</sup>}]
#' @param yeffvol Effective volume of a long core that is measured by the y
#' sensor [positive scientific numeric; in \out{m<sup>3</sup>}]
#' @param zeffvol DEffective volume of a long core that is measured by the z
#' sensor [positive scientific numeric; in \out{m<sup>3</sup>}]
#' @param xeffmass Effective mass of a long core that is measured by the x
#' sensor [positive scientific numeric; in kg].
#' @param yeffmass Effective mass of a long core that is measured by the y
#' sensor [positive scientific numeric; in kg].
#' @param zeffmass Effective mass of a long core that is measured by the z
#' sensor [positive scientific numeric; in kg].
#' @param magsuseffvol Effective volume of a long core that is measured for
#' magnetic susceptibility [positive scientific number; in \out{m<sup>3</sup>}].
#' @param magsuseffmass Effective mass of a long core that is measured for
#' magnetic susceptibility [positive scientific number; in kg].
#' @param deconvoluted Whether the continuous core measurements have been
#' deconvoluted (TRUE) or not (FALSE) [TRUE/FALSE].
#' @param xintconv The original magnetic moment in the x direction of the
#' measuring device, obtained without deconvolution [scientific number; in
#' \out{A m<sup>2</sup>}].
#' @param yintconv The original magnetic moment in the y direction of the
#' measuring device, obtained without deconvolution [scientific number; in
#' \out{A m<sup>2</sup>}].
#' @param zintconv The original magnetic moment in the z direction of the
#' measuring device, obtained without deconvolution [scientific number; in
#' \out{A m<sup>2</sup>}].
#' @param sampaz Sample azimuth in the field (see Fig. 2 in supporting
#' material): it is the angle measured clockwise from the north of the
#' horizontal projection of the field arrow (Tauxe 2010) [normal numeric in
#' [0,360[; in arc degree].
#' @param sampdip Sample dip in the field (see Fig. 2 in supporting material):
#' it is the angle of the  field arrow from the horizontal. It is positive
#' downward, and ranges from +90° for straight down to -90° for straight up
#' (Tauxe, 2010) [normal numeric in [-90,90]; in arc degree].
#' @param samprot Rotation of the sample on its axis, taken in the field (see
#' Fig. 2 in supporting material). It is measured clockwise from the 12 o’clock
#' summit (or X sample coordinate, see Fig. 2B in supporting material), on the
#' sample part opposite of the field arrow. \strong{SPECIAL CASE: IN DOWNWARD
#' VERTICAL SAMPLES} (dip = 90), the rotation is the angle between the azimuth
#' and the field arrow.
#' \strong{SPECIAL CASE: IN UPWARD VERTICAL SAMPLES} (dip = -90), the rotation
#' is the angle between the azimuth and the field arrow + 180° (as the upward
#' dip brings the 12 o’clock position (or X sample coordinate, see Fig. 2B in
#' supporting material) on the sample part opposite of the field arrow, to face
#' the opposite direction of azimuth). [normal numeric in [0,360[; in arc
#' degree].
#' @param coraz Azimuth of the sample in the measuring device (see
#' Fig. 3, 4 & 5 in supplementary material): this is a correction for the
#' difference between the magnetometer coordinates and the sample coordinates.
#' This is equivalent to sample azimuth, considering that the x, y and z
#' coordinates of the magnetometer are attributed to the North, East and
#' downward directions respectively [normal numeric in [0,360[; in arc degree].
#' @param cordip Dip of the sample in the measuring device (see
#' Fig. 3, 4 & 5 in supplementary material):  this is a correction for the
#' difference between the magnetometer coordinates and the sample coordinates.
#' This is equivalent to sample dip, considering that the x, y and z coordinates
#' of the magnetometer are attributed to the North, East and downward directions
#' respectively [normal numeric in [-90,90]; in arc degree].
#' @param corrot Rotation of the sample in the measuring device (see
#' Fig. 3, 4 & 5 in supplementary material):  this is a correction for the
#' difference between the magnetometer coordinates and the sample coordinates.
#' This is equivalent to sample rotation, considering that the x, y and z
#' coordinates of the magnetometer are attributed to the North, East and
#' downward directions respectively [normal numeric in [0,360[; in arc degree].
#' @param bedaz Bedding azimuth, or dip direction (see Fig. 1 in supplementary
#' material): it is the azimuth (the angle taken eastward or clockwise from the
#' north) of the line perpendicular to the plane oriented in the stratigraphic
#' upward direction, projected on a horizontal plane. This accounts for
#' stratigraphically overturned beds (see Fig. 1d-f in supplementary
#' material) [normal numeric in [0,360[; in arc degree].
#' @param beddip Bedding dip (see Fig. 1 in supplementary material): it is the
#' plane's maximum angular deviation from the horizontal. It is positive
#' downward, and ranges from +90° for straight down to 0° for horizontal.
#' Stratigraphically overturned beds are indicated with dip values in the
#' ]90°,180°] interval [normal numeric in [0,180]; in arc degree].
#' @param foldtrend Fold trend (see Fig. 1 in supplementary material): it is the
#' trend of the axis of a fold (or its azimuth, i.e. the angle taken eastward or
#' clockwise from the north) [normal numeric in [0,360[; in arc degree].
#' @param foldplunge Fold plunge (see Fig. 1 in supplementary material): it is
#' the angle that the axis of a fold makes with the horizontal. It is positive
#' downwards, and ranges from +90° form straight down to 0° for horizontal.
#' Stratigraphically overturned folds are indicated with plunge values in the
#' ]90°,180°] interval [normal numeric in [0,180]; in arc degree].
#' @param magaz Magnetic azimuth of the sample in the field [numeric in [0,360[;
#' in arc degree].
#' @param usemagaz If TRUE, the magnetic azimuth (\strong{magaz}) will be used
#' as the sample azimuth (\strong{sampaz}), if FALSE it is the solar azimuth
#' (\strong{solaraz}) that will be used as the sample azimuth
#' (\strong{sampaz}).
#' @param solaraz Solar azimuth of the sample in the field [numeric in [0,360[;
#' in arc degree].
#' @param solarangle Angle that the shadow of a vertical needle does with the
#' core’s field arrow at the moment of sampling (see the \strong{sampmin},
#' \strong{samphour}, \strong{sampday}, \strong{sampmonth} and
#' \strong{sampyear} parameters): this is used to compute the solar azimuth
#' (\strong{solaraz}) [numeric in [0,360[; in arc degree].
#' @param long Longitude [decimal degree in ]-180,180]; in decimal degree].
#' @param lat Latitude [decimal degree in [-90,90]; in decimal degree].
#' @param height Height [numeric; in m].
#' @param sampmin Time of sampling [normal numeric in [0,60[; in minutes].
#' @param samphour Time of sampling [integer in [0,23]; in hours].
#' @param sampday Time of sampling [integer in [1,31]; in days].
#' @param sampmonth Time of sampling [integer in [0,12]; in months].
#' @param sampyear Time of sampling [integer; in years].
#' @param samptimezonemin Time zone of the sampling, expressed relative to
#' the UTC (Coordinated Universal Time), e.g. UTC +12:45 for New Zealand. This
#' has to take into account the possible daylight saving time at the moment of
#' the sampling, which changes the effective time zone [Integer in [0,59]; in
#' minutes].
#' @param samptimezonehour Time zone of the sampling, expressed relative to
#' the UTC (Coordinated Universal Time), e.g. UTC +12:45 for New Zealand. This
#' has to take into account the possible daylight saving time at the moment of
#' the sampling, which changes the effective time zone [integer, usually in
#' [-12,14]; in hours].
#' @param magvar Magnetic variation: it is the angle on the horizontal plane
#' between the magnetic north and the geographic north. It is measured as the
#' angle made by the magnetic North eastward (clockwise) from the geographic
#' North [numeric in [0,360[; in arc degrees].
#' @param magazvarcor Whether the magnetic azimuth parameter (\strong{magaz})
#' has to be corrected for the magnetic variation [TRUE/FALSE].
#' @param bedazvarcor Whether the bedding azimuth parameter (\strong{bedaz})
#' has to be corrected for the magnetic variation. [TRUE/FALSE].
#' @param foldtrendcarcor Whether the fold trend parameter (\strong{foldtrend})
#' has to be corrected for the magnetic variation [TRUE/FALSE].
#' @param treatafx Treatment field by alternating field (AF) demagnetization in
#' the x direction of the measuring device [normal numeric; in T].
#' @param treatafy Treatment field by alternating field (AF) demagnetization in
#' the y direction of the measuring device [normal numeric; in T].
#' @param treatafz Treatment field by alternating field (AF) demagnetization in
#' the z direction of the measuring device [normal numeric; in T].
#' @param treattemp Treatment temperature in Kelvin [positive normal numeric;
#' in K].
#' @param treatirmx Treatment field by Isothermal Remanent Magnetisation (IRM)
#' in the x direction of the measuring device [normal numeric; in T].
#' @param treatirmy Treatment field by Isothermal Remanent Magnetisation (IRM)
#' in the y direction of the measuring device [normal numeric; in T].
#' @param treatirmz Treatment field by Isothermal Remanent Magnetisation (IRM)
#' in the z direction of the measuring device [normal numeric; in T].
#' @param treatarmafx Treatment anhysteretic field in the x direction of the
#' measuring device for Anhysteretic Remanent Magnetization (ARM) [normal
#' numeric; in T].
#' @param treatarmafy Treatment anhysteretic field in the y direction of the
#' measuring device for Anhysteretic Remanent Magnetization (ARM) [normal
#' numeric; in T].
#' @param treatarmafz Treatment anhysteretic field in the z direction of the
#' measuring device for Anhysteretic Remanent Magnetization (ARM) [normal
#' numeric; in T].
#' @param treatarmbiasx Treatment bias field in x direction of the measuring
#' device for Anhysteretic Remanent Magnetization (ARM). In most settings the
#' bias field comes from the z axis of the measuring device, this parameter is
#' only set to allow for unconventional settings [normal numeric; in T].
#' @param treatarmbiasy Treatment bias field in y direction of the measuring
#' device for Anhysteretic Remanent Magnetization (ARM). In most settings the
#' bias field comes from the z axis of the measuring device, this parameter is
#' only set to allow for unconventional settings [normal numeric; in T].
#' @param treatarmbiasz Treatment bias field in z direction for Anhysteretic
#' Remanent Magnetization (ARM) [normal numeric; in T].
#' @param pcaanchor Whether to anchor the axis of a Principal Component Analysis
#' (PCA) at the origin (x = 0, y = 0, z = 0) [TRUE/FALSE].
#' @param pcacomponent Identifies groups of measurements considered in one
#' component (e.g. “C1”, “C2”, “Outlier”), to be used for Principal Component
#' Analysis (PCA) [alphanumeric].
#' @param pcacomponentsingle Which group of measurements to use for Principal
#' Component Analysis (PCA) in files only considering one single group
#' [alphanumeric].
#' @param circlecomponent Identifies groups of measurements considered in one
#' component (e.g. “C1”, “C2”, “Outlier”), to be used for circle computation
#' [alphanumeric].
#' @param circlecomponentsingle Which group of measurements to use for circle
#' computation in files only considering one single group. [alphanumeric].
#' @param comments Additional comments [alphanumeric].
#' @param originalfile The name of the original file from which the data is from
#' [alphanumeric].
#' @param originalfileline The line, in the original file, of the sample
#' identification information (used for \strong{sampid}) related to the x
#' magnetic moment data. If not available, the line of the sample identification
#' related to the x magnetisation by volume or by mass (in this order). This is
#' to improve the manual or automated localisation of mistakes in the original
#' file, based on automated or manual checks of the pmob.
#'
#' @references
#' \itemize{
#'   \item Tauxe, L., 2010. Essentials of Paleomagnetism. University of
#'   California Press.
#' }
#'
#' @examples
#' as.pmob(exc = 1, sampid = c("Bosso1", "MonteAcuto1"))
#'
#' @export
#' @importFrom StratigrapheR homogenise merge_list
#' @import lubridate

as.pmob <- function(...,
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

                    blankcor = NULL,
                    isblank  = NULL,
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

  pm.list <- list(sampid        = sampid,
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

                  blankcor = blankcor,
                  isblank  = isblank,
                  blankmethod = blankmethod,

                  holdercor = holdercor,
                  isholder  = isholder,
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

                  originalfile     = originalfile,
                  originalfileline = originalfileline)

  pm.length <- sapply(pm.list, length)

  null.length.pos <- which(pm.length == 0L)
  pres.length.pos <- which(pm.length != 0L)

  n <- max(pm.length)
  l <- pm.list[-null.length.pos]

  ext.pm <- homogenise(n = n, l = l, cycle = FALSE)

  proto <- pm.list

  proto[pres.length.pos] <- ext.pm

  # Check functions ----

  pm.alphanum <- function(arg) if(!is.null(arg)) return(as.character(arg))

  pm.integer <- function(arg, header,
                         xmin = NULL, min.in = TRUE,
                         xmax = NULL, max.in = TRUE,
                         unique = FALSE)
  {
    if(!is.null(arg)) {

      if(!(inherits(arg, c("numeric", "integer")) | all(is.na(arg)))){
        stop("The '", header,
             "' parameter should be of class 'numeric' or 'integer'",
             call. = FALSE)
      }

      if(!is.null(xmin)) {

        if(isTRUE(any(arg < xmin))) {
          stop("Values in '", header, "' should be higher than ", xmin,
               call. = FALSE)
        }

        if(!min.in & isTRUE(any(xmin == arg))){
          stop("Values in '", header, "' should not be equal to ", xmin,
               call. = FALSE)
        }

      }

      if(!is.null(xmax)) {

        if(isTRUE(any(arg > xmax))) {
          stop("Values in '", header, "' should be lower than ", xmax,
               call. = FALSE)
        }

        if(!max.in & isTRUE(any(xmax == arg))){
          stop("Values in '", header, "' should not be equal to ", xmax,
               call. = FALSE)
        }

      }

      if(!all.equal(round(arg,0), round(arg,16))) {
        warning("Values in '", header, "' should be strict integers",
                call. = FALSE)
      }

      if(unique) {
        if(any(duplicated(arg))& !is.na(arg)) {
          stop("The values in '", header,
               "' should not be repeated in any way",
               call. = FALSE)
        }
      }

      return(as.integer(arg))
    }
  }

  pm.norm.num <- function(arg, header,
                          xmin = NULL, min.in = TRUE,
                          xmax = NULL, max.in = TRUE)
  {

    if(!is.null(arg)) {
      if(!(inherits(arg, c("numeric", "integer")) | all(is.na(arg)))){
        stop("The '", header,
             "' parameter should be of class 'numeric' or 'integer'",
             call. = FALSE)
      }


      if(!is.null(xmin)) {

        if(isTRUE(any(arg < xmin))) {
          stop("Values in '", header, "' should be higher than ", xmin,
               call. = FALSE)
        }

        if(!min.in & isTRUE(any(arg == xmin))){
          stop("Values in '", header, "' should not be equal to ", xmin,
               call. = FALSE)
        }

      }

      if(!is.null(xmax)) {

        if(isTRUE(any(arg > xmax))) {
          stop("Values in '", header, "' should be lower than ", xmax,
               call. = FALSE)
        }

        if(!max.in & isTRUE(any(xmax == arg))){
          stop("Values in '", header, "' should not be equal to ", xmax,
               call. = FALSE)
        }

      }

      return(as.numeric(arg))

    }
  }

  pm.sci.num <- function(arg, header, pos = FALSE)
  {
    if(!is.null(arg)) {

      if(!(inherits(arg, c("numeric", "integer")) | all(is.na(arg)))){
        stop("The '", header,
             "' parameter should be of class 'numeric' or 'integer'",
             call. = FALSE)
      }

      if(pos & isTRUE(any(arg < 0))) {
        stop("The values in '", header, "' should be positive", call. = FALSE)
      }

      return(as.numeric(arg))
    }
  }

  pm.boolean <- function(arg, header)
  {
    if(!is.null(arg)) {

      if(!(inherits(arg, "logical") | all(is.na(arg)))){
        stop("The '", header,
             "' parameter should be of class 'logical'",
             call. = FALSE)
      }

      return(as.logical(arg))

    }
  }

  # Accumulation into a pmob object ----

  pmob <- list()

  # ----

  pmob$sampid <- pm.alphanum(proto$sampid)
  pmob$specid <- pm.alphanum(proto$specid)
  pmob$slotid <- pm.alphanum(proto$slotid)
  pmob$measid <- pm.alphanum(proto$measid)
  pmob$siteid <- pm.alphanum(proto$siteid)

  if(is.null(proto$initialorder)){
    pmob$initialorder <- as.integer(seq_len(n))
  } else {
    pmob$initialorder  <- pm.integer(proto$initialorder, "initialorder",
                                     xmin = 1L,
                                     unique = T)
  }

  pmob$depth <- pm.norm.num(proto$depth, "depth")

  pmob$meassec <- pm.norm.num(proto$meassec, "meassec",
                              xmin = 0,
                              xmax = 60, max.in = FALSE)

  pmob$measmin <- pm.integer(proto$measmin, "measmin",
                             xmin = 0L,
                             xmax = 60L, max.in = FALSE)

  pmob$meashour <- pm.integer(proto$meashour, "meashour",
                              xmin = 0L,
                              xmax = 24L, max.in = FALSE)

  pmob$measday <- pm.integer(proto$measday, "measday",
                             xmin = 0L, xmax = 31L)

  pmob$measmonth <- pm.integer(proto$measmonth, "measmonth",
                               xmin = 0L, xmax = 12L)

  pmob$measyear <- pm.integer(proto$measyear, "measyear")

  pmob$measdevice   <- pm.alphanum(proto$measdevice)
  pmob$magsusdevice <- pm.alphanum(proto$magsusdevice)

  pmob$xint <- pm.sci.num(proto$xint, "xint")
  pmob$yint <- pm.sci.num(proto$yint, "yint")
  pmob$zint <- pm.sci.num(proto$zint, "zint")

  pmob$blankcor    <- pm.boolean(proto$blankcor, "blankcor")
  pmob$isblank     <- pm.boolean(proto$isblank,   "isblank")
  pmob$blankmethod <- pm.alphanum(proto$blankmethod)

  pmob$holdercor    <- pm.boolean(proto$holdercor, "holdercor")
  pmob$isholder     <- pm.boolean(proto$isholder,   "isholder")
  pmob$holdermethod <- pm.alphanum(proto$holdermethod)

  pmob$xinterror <- pm.sci.num(proto$xinterror, "xinterror")
  pmob$yinterror <- pm.sci.num(proto$yinterror, "yinterror")
  pmob$zinterror <- pm.sci.num(proto$zinterror, "zinterror")

  pmob$xvol <- pm.sci.num(proto$xvol, "xvol")
  pmob$yvol <- pm.sci.num(proto$yvol, "yvol")
  pmob$zvol <- pm.sci.num(proto$zvol, "zvol")

  pmob$xmass <- pm.sci.num(proto$xmass, "xmass")
  pmob$ymass <- pm.sci.num(proto$ymass, "ymass")
  pmob$zmass <- pm.sci.num(proto$zmass, "zmass")

  pmob$totmagsus <- pm.sci.num(proto$totmagsus, "totmagsus")

  pmob$totmagsuserror <- pm.sci.num(proto$totmagsuserror, "totmagsuserror")

  pmob$volmagsus  <- pm.sci.num(proto$volmagsus,  "volmagsus")
  pmob$massmagsus <- pm.sci.num(proto$massmagsus, "massmagsus")

  pmob$vol  <- pm.sci.num(proto$vol,  "vol",      pos = T)
  pmob$mass <- pm.sci.num(proto$mass, "mass",    pos = T)
  pmob$dens <- pm.sci.num(proto$dens, "dens", pos = T)

  pmob$exactvol     <- pm.boolean(proto$exactvol, "exactvol")
  pmob$exactmass    <- pm.boolean(proto$exactmass, "exactmass")
  pmob$exactdens <- pm.boolean(proto$exactdens, "exactdens")

  pmob$discrete <- pm.boolean(proto$discrete, "discrete")

  pmob$area <- pm.sci.num(proto$area, "area", pos = T)

  pmob$defaultarea <- pm.boolean(proto$defaultarea, "defaultarea")

  pmob$xefflength <- pm.sci.num(proto$xefflength, "xefflength", pos = T)
  pmob$yefflength <- pm.sci.num(proto$yefflength, "yefflength", pos = T)
  pmob$zefflength <- pm.sci.num(proto$zefflength, "zefflength", pos = T)

  pmob$xeffvol <- pm.sci.num(proto$xeffvol, "xeffvol", pos = T)
  pmob$yeffvol <- pm.sci.num(proto$yeffvol, "yeffvol", pos = T)
  pmob$zeffvol <- pm.sci.num(proto$zeffvol, "zeffvol", pos = T)

  pmob$xeffmass <- pm.sci.num(proto$xeffmass, "xeffmass", pos = T)
  pmob$yeffmass <- pm.sci.num(proto$yeffmass, "yeffmass", pos = T)
  pmob$zeffmass <- pm.sci.num(proto$zeffmass, "zeffmass", pos = T)

  pmob$magsuseffvol  <- pm.sci.num(proto$magsuseffvol,  "magsuseffvol",  pos = T)
  pmob$magsuseffmass <- pm.sci.num(proto$magsuseffmass, "magsuseffmass", pos = T)

  pmob$deconvoluted <- pm.boolean(proto$deconvoluted, "deconvoluted")

  pmob$xintconv <- pm.sci.num(proto$xintconv, "xintconv")
  pmob$yintconv <- pm.sci.num(proto$yintconv, "yintconv")
  pmob$zintconv <- pm.sci.num(proto$zintconv, "zintconv")

  pmob$sampaz <- pm.norm.num(proto$sampaz, "sampaz",
                             xmin = 0,
                             xmax = 360, max.in = F)

  pmob$sampdip <- pm.norm.num(proto$sampdip, "sampdip",
                              xmin = -90,
                              xmax = 90)

  pmob$samprot <- pm.norm.num(proto$samprot, "samprot",
                              xmin = 0,
                              xmax = 360, max.in = F)

  pmob$coraz <- pm.norm.num(proto$coraz, "coraz",
                            xmin = 0,
                            xmax = 360, max.in = F)

  pmob$cordip <- pm.norm.num(proto$cordip, "cordip",
                             xmin = -90,
                             xmax = 90)

  pmob$corrot <- pm.norm.num(proto$corrot, "corrot",
                             xmin = 0,
                             xmax = 360, max.in = F)

  pmob$bedaz <- pm.norm.num(proto$bedaz, "bedaz",
                            xmin = 0,
                            xmax = 360, max.in = F)

  pmob$beddip <- pm.norm.num(proto$beddip, "beddip",
                             xmin = 0,
                             xmax = 180)

  pmob$foldtrend <- pm.norm.num(proto$foldtrend, "foldtrend",
                                xmin = 0,
                                xmax = 360, max.in = F)

  pmob$foldplunge <- pm.norm.num(proto$foldplunge, "foldplunge",
                                 xmin = 0,
                                 xmax = 180)

  pmob$magaz <- pm.norm.num(proto$magaz, "magaz",
                            xmin = 0,
                            xmax = 360, max.in = F)

  pmob$usemagaz <- pm.boolean(proto$usemagaz, "usemagaz")

  pmob$solaraz <- pm.norm.num(proto$solaraz, "solaraz",
                              xmin = 0,
                              xmax = 360, max.in = F)

  pmob$solarangle <- pm.norm.num(proto$solarangle, "solarangle",
                                 xmin = 0,
                                 xmax = 360, max.in = F)

  pmob$long <- pm.norm.num(proto$long, "long",
                           xmin = -180, min.in = F,
                           xmax = 180)

  pmob$lat <- pm.norm.num(proto$lat, "lat",
                          xmin = -90,
                          xmax = 90)

  pmob$height <- pm.norm.num(proto$height, "height")

  pmob$sampmin <- pm.norm.num(proto$sampmin, "sampmin",
                              xmin = 0,
                              xmax = 60, max.in = FALSE)

  pmob$samphour <- pm.integer(proto$samphour, "samphour",
                              xmin = 0L,
                              xmax = 24L, max.in = FALSE)

  pmob$sampday <- pm.integer(proto$sampday, "sampday",
                             xmin = 0L, xmax = 31L)

  pmob$sampmonth <- pm.integer(proto$sampmonth, "sampmonth",
                               xmin = 0L, xmax = 12L)

  pmob$sampyear <- pm.integer(proto$sampyear, "sampyear")

  pmob$samptimezonemin <- pm.integer(proto$sampimezonemin,
                                     "samptimezonemin",
                                     xmin = 0L,
                                     xmax = 60L, max.in = FALSE)

  pmob$samptimezonehour <- pm.integer(proto$samptimezonehour,
                                      "samptimezonehour")

  pmob$magvar <- pm.norm.num(proto$magvar, "magvar",
                             xmin = 0,
                             xmax = 360, max.in = F)

  pmob$magazvarcor     <- pm.boolean(proto$magazvarcor,     "magazvarcor")
  pmob$bedazvarcor     <- pm.boolean(proto$bedazvarcor,     "bedazvarcor")
  pmob$foldtrendvarcor <- pm.boolean(proto$foldtrendvarcor, "foldtrendvarcor")

  pmob$treatafx <- pm.norm.num(proto$treatafx, "treatafx")
  pmob$treatafy <- pm.norm.num(proto$treatafy, "treatafy")
  pmob$treatafz <- pm.norm.num(proto$treatafz, "treatafz")

  pmob$treattemp <- pm.norm.num(proto$treattemp, "treattemp", xmin = 0)

  pmob$treatirmx <- pm.norm.num(proto$treatirmx, "treatirmx")
  pmob$treatirmy <- pm.norm.num(proto$treatirmy, "treatirmy")
  pmob$treatirmz <- pm.norm.num(proto$treatirmz, "treatirmz")

  pmob$treatarmafx <- pm.norm.num(proto$treatarmafx, "treatarmxafx")
  pmob$treatarmafy <- pm.norm.num(proto$treatarmafy, "treatarmxafy")
  pmob$treatarmafz <- pm.norm.num(proto$treatarmafz, "treatarmxafz")

  pmob$treatarmbiasx <- pm.norm.num(proto$treatarmbiasx, "treatarmbiasx")
  pmob$treatarmbiasy <- pm.norm.num(proto$treatarmbiasy, "treatarmbiasy")
  pmob$treatarmbiasz <- pm.norm.num(proto$treatarmbiasz, "treatarmbiasz")

  pmob$pcaanchor <- pm.boolean(proto$pcaanchor, "pcaanchor")

  pmob$pcacomponent       <- pm.alphanum(proto$pcacomponent)
  pmob$pcacomponentsingle <- pm.alphanum(proto$pcacomponentsingle)

  pmob$circlecomponent       <- pm.alphanum(proto$circlecomponent)
  pmob$circlecomponentsingle <- pm.alphanum(proto$circlecomponentsingle)

  pmob$comments <- pm.alphanum(proto$comments)

  pmob$originalfile <- pm.alphanum(proto$originalfile)
  pmob$originalfileline <- pm.integer(proto$originalfileline,
                                      "originalfileline", xmin = 1)

  pmob$version <- rep("pmob 0.0.0.9020", n)

  # Validity of measurement date ----

  if(any(pmob$measyear + 0.1 < 1833L, na.rm = T)) {
    warning("The first usable magnetometer was made in 1833",
            " by Carl Friedrich Gauss: \nAre you sure that your",
            " measurements are",
            " older than that ?")
  }

  if(!is.null(pmob$meaday) &
     !is.null(pmob$measmonth) &
     !is.null(pmob$measyear)){

    meas.rawdate <- paste(pmob$measday,
                          pmob$measmonth,
                          pmob$measyear, sep = "-")

    meas.redate <- as.Date(meas.rawdate,
                           format("%d-%m-%Y"))

    meas.na.date <- is.na(meas.redate)

    if(any(meas.na.date &
           !is.na(pmob$measday) &
           ! is.na(pmob$measmonth) &
           !is.na(pmob$measyear))){

      meas.invalid.pos <- which(meas.na.date &
                                  !is.na(pmob$measday) &
                                  ! is.na(pmob$measonth) &
                                  !is.na(pmob$measyear))

      stop("The following measurement dates (day-month-year) are not valid:\n",
           paste(unique(meas.rawdate[meas.invalid.pos]), collapse = "\n"))
    }
  }

  # Validity of sampling date ----

  if(any(pmob$sampyear + 0.1 < -1500L, na.rm = T)) {
    warning("The first known shadow clocks were made in 1500 BC,",
            " in Mesopotamia.",
            " \nAlso: what the hell ??? No way you have samples collected",
            " in the year ",
            abs(min(pmob$measyear)), " BC")
  }

  if(!is.null(pmob$sampday) &
     !is.null(pmob$sampmonth) &
     !is.null(pmob$sampyear)){

    meas.rawdate <- paste(pmob$sampday,
                          pmob$sampmonth,
                          pmob$sampyear, sep = "-")

    meas.redate <- as.Date(meas.rawdate,
                           format("%d-%m-%Y"))

    meas.na.date <- is.na(meas.redate)

    if(any(meas.na.date &
           !is.na(pmob$sampday) &
           ! is.na(pmob$sampmonth) &
           !is.na(pmob$sampyear))){

      meas.invalid.pos <- which(meas.na.date &
                                  !is.na(pmob$sampday) &
                                  ! is.na(pmob$sampmonth) &
                                  !is.na(pmob$sampyear))

      stop("The following sampling dates (day-month-year) are not valid:\n",
           paste(unique(meas.rawdate[meas.invalid.pos]), collapse = "\n"))
    }
  }

  # ----

  res <- data.frame(pmob, stringsAsFactors = F)

  glob <- list(...)

  l.glob <- length(glob)

  is.list.glob <- lapply(glob, function(x) inherits(x, "list"))

  if(l.glob == 1){

    if(is.list.glob[[1]]){

      list.glob <- glob[[1]]

    } else {

      list.glob <- list(...)

    }

  } else {

    list.glob <- list(...)

  }

  add <- as.data.frame(homogenise(n = n, l = list.glob, cycle = F),
                       stringsAsFactors = F)

  if(nrow(add) == n) res <- cbind(res, add)

  return(res)


}
