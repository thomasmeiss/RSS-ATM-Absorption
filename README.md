# RSS-ATM-Absorption
atmospheric absorption models for water vapor (H2O) and dry air

12/21/2022
by Thomas Meissner
Remote Sensing Systems

RSS atmospheric absorption models for water vapor (H2O) and dry air at microwave frequencies

based on: 

[1] F. Wentz and T. Meissner, 2016, Atmospheric absorption model for dry air and water vapor 
      at microwave frequencies below 100GHz derived from spaceborne radiometer observations, 
      Radio Science 51, doi:10.1002/2015RS005858. 

[2]  K. Wentz, T. Meissner, F. Wentz, and A. Manaster, 2024, Absolute Intercalibration of 
     Spaceborne Microwave Radiometers. J. Atmos. Oceanic Technol., 41, 1121–1138, 
     https://doi.org/10.1175/JTECH-D-24-0024.1.     

[3] P. Rosenkranz, 2017, Line-by-line microwave radiative transfer (non-scattering) Fortran programs, 
       doi: 10.21982/M81013. Available online at https://rscl-grss.org/ 2017-05-15.   

major changes from [1] to [2]
1. water vapor continuum absoprtion
   use Rosenkranz 2019 FB [2]
   use Rosenkranz 2019 SB [2] , reduced by 7%
2. non-resonant oxygen absoprtion [1]
   the width parameter gamma_c has temperature dependence:
   T^-1.5 for C-band and above
   T^-0.8 for L-band
   smooth transition (sigmoid function) between L-band and C-band
