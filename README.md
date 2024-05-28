# RSS-ATM-Absorption-
atmospheric absorption models for water vapor (H2O) and dry air

12/21/2022
by Thomas Meissner
Remote Sensing Systems

RSS atmospheric absorption models for water vapor (H2O) and dry air at microwave frequencies

based on: 

[1] F. Wentz and T. Meissner, 2016, Atmospheric absorption model for dry air and water vapor 
      at microwave frequencies below 100GHz derived from spaceborne radiometer observations, 
      Radio Science 51, doi:10.1002/2015RS005858. 

[2] P. Rosenkranz, 2017, Line-by-line microwave radiative transfer (non-scattering) Fortran programs, 
       doi: 10.21982/M81013. Available online at https://rscl-grss.org/ 2017-05-15.   

major changes from [1]
1. water vapor continuum absoprtion
   use Rosenkranz 2019 FB [2]
   use Rosenkranz 2019 SB [2] , reduced by 7%
2. non-resonant oxygen absoprtion [1]
   the width parameter gamma_c has temperature dependence:
   T^-1.5 for C-band and above
   T^-0.8 for L-band
   smooth transition (sigmoid function) between L-band and C-band
