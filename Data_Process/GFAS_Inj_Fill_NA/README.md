README -- describes contents of HEMCO/GFAS/v2021-09 directory
10 September 2021
Lixu Jin, University of Montana
lixu.jin@umontana.edu

Overview:
================================================================================

This directory contains injection heights derived from
fire observations and meteorological information from the operational weather
forecasts of ECMWF.

Injection heights are weighted by grids with emissions to tackle with grid-dependent issue. 
Grids where emissions equals 0 will be masked. 

The processing script is provided in ./GFAS_injs.ipynb
The script could be customized for users' purposes (e.g., different injections/resolutions).

FRP observations currently assimilated in this data set are the NASA Terra MODIS
and Aqua MODIS active fire products - 

              (http://modis-fire.umd.edu/pages/ActiveFire.php)

This data set includes: Fire Radiative Power (FRP), dry matter burnt and biomass
burning emissions.

Citation:
--------------------------------------------------------------------------------
Contains modified Copernicus Atmosphere Monitoring Service Information 2018

Please acknowledge the use of this data set according to the terms of the
Copernicus CAMS License agreement:

"Where the Licensee makes or contributes to a publication or distribution
containing adapted or modified CAMS Information, the Licensee shall provide the
following or any similar notice:"

'Contains modified Copernicus Atmosphere Monitoring Service Information [Year]';

Reference:
--------------------------------------------------------------------------------
Francesca Di Giuseppe, Samuel RÃ©my, Florian Pappenberger, and Fredrik
Wetterhall: Combining fire radiative power observations with the fire weather
index improves the estimation of fire emissions, Atmos. Chem. Phys. Discuss.,
https://doi.org/10.5194/acp-2017-790

RÃ©my, S., Veira, A., Paugam, R., Sofiev, M., Kaiser, J. W., Marenco, F., Burton,
S. P., Benedetti, A., Engelen, R. J., Ferrare, R., and Hair, J. W.: Two global
data sets of daily fire emission injection heights since 2003, Atmos. Chem.
Phys., 17, 2921-2942, https://doi.org/10.5194/acp-17-2921-2017, 2017

N. Andela (VUA), J.W. Kaiser (ECMWF, KCL), A. Heil (FZ JÃ¼lich), T.T. van Leeuwen
(VUA), G.R. van der Werf (VUA), M.J. Wooster (KCL), S. Remy (ECMWF) and
M.G. Schultz (FZ JÃ¼lich), Assessment of the Global Fire Assimilation System
(GFASv1)

Kaiser, J. W., Heil, A., Andreae, M. O., Benedetti, A., Chubarova, N., Jones,
L., Morcrette, J.-J., Razinger, M., Schultz, M. G., Suttie, M., and van der
Werf, G. R. (2012). Biomass burning emissions estimated with a global fire
assimilation system based on observed fire radiative power. BG, 9:527-554

Xu et al. (2010) New GOES imager algorithms for cloud and active fire detection
and fire radiative power assessment across North, South and Central America.
RSE Vol. 114

Heil et al. (2010) Assessment of the Real-Time Fire Emissions (GFASv0) by MACC,
ECMWF Tech. Memo No. 628

Di Giuseppe, F, Remy, S, Pappenberger, F, Wetterhall, F (2016): Improving GFAS
and CAMS biomass burning estimations by means of the Global ECMWF Fire Forecast
system (GEFF), ECMWF Tech. Memo No. 790

Files:
================================================================================

YYYY/GFAS_YYYYMM_hgts_RES.nc

 -- Gridded plume rise model parameters, GFAS analysis surface parameters,
    gridded satellite parameters

    Resolution  : 4x5/0.25x0.3125
    Units       : m (plume rise model parameters)
    Timestamps  : Monthly, 2003/01 through present month
    Compression : Level 4 (i.e. nccopy -d4)

    Exhaustive data set information can be found on the CAMS GFAS page in the
    ECMWF wiki:

                          https://confluence.ecmwf.int/

How the NetCDF files were created:
================================================================================

The data set files were created as follows.

1) Retrieve NetCDF format data in half month chunks from ECMWF.

2) Create new NetCDF format data file to contain entire month's worth of data,
   setting 'title', 'conventions', and 'history' attributes as per COARDS.

3) Create output dimensions, converting 'time' from 

                      hours since 1900-01-01 00:00:0.0

   to

                      hours since 1970-01-01 00:00:0.0

   and reversing the sequence of latitude values.

4) For each data set variable:
       - set 'units', 'long_name', 'missing_value'
       - set the modal value in each input half-monthly variable to zero
         (if it is an emission variable), or the output missing value
         otherwise. There does not appear to be a consistent missing value
         in the output from the ECMWF data API for this data set in NetCDF
         format; as this is sparse data, the modal value will be missing data
       - concatenate input half-monthly data, reversing the latitude dimension

4) For the mean altitude of maximum injection variable:
       - where there is no CO emission, i.e. no fire, mask mean altitude of
         maximum injection value
       - where is CO emission, i.e. fire, average mean altitude of maximum
         injection value at GEOS-Chem resolution
       - the above two steps were taken to make data more HEMCO-friendly.
================================================================================
