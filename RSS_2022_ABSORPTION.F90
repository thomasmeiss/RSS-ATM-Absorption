! atmospheric absorption models for water vapor (H2O) and dry air
! 
! updated 12/21/2022
! by Thomas Meissner
! see 
! [1]  \\THOMAS\C\atmosphere\Atmospheric Absorption Models (H2O, O2).pptx - Shortcut
! [2]  Frank's memo \\THOMAS\C\atmosphere\memo\memo4.txt
! [3]  F. Wentz and T. Meissner, 2016, Atmospheric absorption model for dry air and water vapor 
!      at microwave frequencies below 100GHz derived from spaceborne radiometer observations, 
!      Radio Science 51, doi:10.1002/2015RS005858. 
! [4]  P. Rosenkranz, 2017, Line-by-line microwave radiative transfer (non-scattering) Fortran programs, 
!      doi: 10.21982/M81013. Available online at https://rscl-grss.org/ 2017-05-15.   

! based on Frank's absorption models for the RSS V8 processor
! major changes
! 1. water vapor continuum absoprtion
!    use Rosenkranz 2019 FB
!    use Rosenkranz 2019 SB, reduced by 7%
! 2. non-resonant oxygen absoprtion
!    the width parameter gamma_c has temperature dependence:
!    T^-1.5 for C-band and above
!    T^-0.8 for L-band
!    smooth transition (sigmoid function) between L-band and C-band

module RSS_2022_ABSORPTION

public
save

contains


subroutine absh2o_RSS_2022(p,t,pv,freq,  absh2o, absh2o_comp)
! model [3] with modifications [2] and [1].
implicit none

real(4), intent(in)            :: p,t,pv,freq
! units mb Kelvin mb GHz

real(4), intent(out)           :: absh2o
real(4), optional, intent(out) :: absh2o_comp(3) 
! units neper/km

real(4)                        :: alphac(3)
  
integer(4), parameter          :: nlines=15

integer(4), save               :: istart=1
integer(4)                     :: i
real(4)                        :: s,base
real(4)                        :: tht,pwet,pdry,ga,gasq,xterm,rnuneg,rnupos
real(8)                        :: sum
real(4)                        :: chi,chisq,freqsq,f0sq,u

real(4), parameter             :: db_to_neper = 0.2302585094 !convert db/km to neper/km

! SB adjustment
real(4), parameter             :: zeta_SB=0.93   !0.83 ! strength of SB continuum


!     line frequencies:
real(4), dimension(nlines), parameter :: &
f0 = (/22.2351, 183.3101, 321.2256, 325.1529, 380.1974, 439.1508, &
443.0183, 448.0011, 470.8890, 474.6891, 488.4911, 556.9360,       &
620.7008, 752.0332, 916.1712/)

!     line intensities at 300k:
real(4), dimension(nlines), save :: b1


!     t coeff. of intensities:
real(4), dimension(nlines), save :: b2

!     air-broadened width parameters at 300k:
real(4), dimension(nlines), save :: b3

!     self-broadened width parameters at 300k:
real(4), dimension(nlines), save :: b5

!     t-exponent of air-broadening:
real(4), dimension(nlines), save :: b4

!     t-exponent of self-broadening:
real(4), dimension(nlines), save :: b6


if(istart.eq.1) then
	  istart=0

      b1 = (/.1310e-13, .2273e-11, .8036e-13, .2694e-11, .2438e-10,             & 
              .2179e-11, .4624e-12, .2562e-10, .8369e-12, .3263e-11, .6659e-12, &
              .1531e-08, .1707e-10, .1011e-08, .4227e-10/)	 
       
      b2 = (/2.144, .668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405, 3.597, 2.379, 2.852, .159, 2.391, .396, 1.441/)

      b3 = (/.0281, .0281, .023, .0278, .0287, .021, .0186, .0263, .0215, .0236, .026, .0321, .0244, .0306, .0267/)
 
      b5 = (/.1349, .1491, .108, .135, .1541, .090, .0788, .1275, .0983, .1095, .1313, .1320, .1140, .1253, .1275/)
      
      b4  = (/.69, .64, .67, .68, .54, .63, .60, .66, .66, .65, .69, .69, .71, .68, .70/)
 
      b6 = (/.61, .85, .54, .74, .89, .52, .50, .67, .65, .64, .72, 1.0, .68, .84, .78/)
        
      b1=1.8281089E+14*b1/f0**2
	  b5=b5/b3  !convert b5 to Liebe notation
      b3(1)=b3(1)/1.040 
      !modification 1: reduce air broadened line width of 22.24 GHz by 4% (Payne et al. 2008, RK 2019)
      b3(2)=b3(2)*1.03  
      !increase air broadend line width of 183.81 GHz by 3% (Payne et al. 2008), bug fixed 05/09/2023 (Katherine)
      
endif
 
if(pv.le.0.) then
      absh2o=0.0
      if (present(absh2o_comp)) absh2o_comp=0.0
	  return
endif
 

pwet=0.1*pv
pdry=0.1*p-pwet
tht = 300./t
xterm=1-tht
freqsq=freq*freq

sum = 0.

do i=1,nlines ! sum over lines

f0sq=f0(i)*f0(i)
ga=b3(i)*(pdry*tht**b4(i) + b5(i)*pwet*tht**b6(i))
gasq = ga*ga
s = b1(i)*exp(b2(i)*xterm)
rnuneg = f0(i)-freq
rnupos = f0(i)+freq
base = ga/(562500. + gasq)  !use clough's definition of local line contribution

if(i.ne.1) then
      if(abs(rnuneg).lt.750) sum = sum + s*(ga/(gasq + rnuneg**2) - base)
      if(abs(rnupos).lt.750) sum = sum + s*(ga/(gasq + rnupos**2) - base)

else
      chi=0.07*ga   !modification 2: change of wings of 22 GHz line. see Frank's memo4.txt [2]

	  if(freq.lt.19) then
	  u=abs(freq-19.)/16.5
	  if(u.lt.0) u=0
	  if(u.gt.1) u=1
	  chi=0.07*ga + 0.93*ga*u*u*(3-2*u)  !see Frank's memo4.txt [2] 
	  endif	
	  chisq=chi*chi
      sum=sum +     s*2*((ga-chi)*freqsq + (ga+chi)*(f0sq+gasq-chisq))/((freqsq-f0sq-gasq+chisq)**2 + 4*freqsq*gasq)
endif

enddo ! sum over lines

if(sum.lt.0) sum=0


alphac(1) = pwet*freq*(tht**3.5)*sum                                                  ! line
alphac(2) = (10*pwet)*(10*pdry)*freq*(tht**3)*5.94E-10/(db_to_neper*0.1820)           ! RK 2019 FB [4]
alphac(3) = zeta_SB*(10*pwet)*(10*pwet)*freq*(tht**7.5)*1.42E-8/(db_to_neper*0.1820)  ! RK 2019 SB * zeta [4],[1]
alphac    = (db_to_neper*0.1820)*freq*alphac

absh2o = alphac(1) + alphac(2) + alphac(3)

if (present(absh2o_comp)) absh2o_comp = alphac

return
end subroutine absh2o_RSS_2022


subroutine absoxy_RSS_2022(p,t,pv,freq, absoxy)
! total dry air abs
! has N2 continuum added
! model [3] with modifications [2] and [1].


!     Revision History 

!     TM 12/21/2022: non-resonant O2 width temperatur dependenc changes smoothly frpm theta^0.8 at L-band to tht^1.5 at C-band.

!     original name: fdabsoxy_1992_modified
!     july 17 2015  changed oct 11 2016.  changes made to vapor model. see 'memo4.txt'
!     july 17 2015 veresion save in O:\skytemp3\version_07172015
!
!     june 11 2015 changed july 17 2015.  oxyopc adjustment above 37 ghz changed slightly. see 'O:\gmi\abs_cal\memo20.txt'

!     modified version of Liebe 1992 oxygen model 
!     from Atmospheric 60-GHz Oxygen Spectrum:.. Liebe, Rosenkranz, Hufford, 1992
!     coded: June 2009 1992 by f.wentz
!           inputs: t, temperature (k)
!                   p, total pressure (mb)
!                   pv, water vapor pressure (mb)
!                   freq, frequency (ghz)
!           output: absoxy, oxygen absorption coefficient (neper/km)
!
!
!     It is the same as fdabsoxy_1989 except for the a5 and a6 coefs for finding delta have different values.
!     Also this 1992 versions says delta is proprotional to total pressure p rather than pdry. 
!     Also note in this routine, the 1.e-3 scaling is done in the startup block.
!
!     compared to abo2_rk (the code Rosenkranz sent us), this routine gives very similar results if you set the apterm to 0.
!     for my freqs 6-85 ghz, the largest dif was 0.0003 at the coldest vapor at 85.5 GHz.
!     the apterm adds 0.003 at 85 ghz
!     Apart from the apterm, you can essentially say fdabsoxy_1992 and abo2_rk are the same for my purposes.
!
!     this routine has been modified in the following ways (see 'memo8.txt')
!     1.  non-resonance continuum temperature coef changed from 0.8 to 1.5
!     2.  a p*p continuum was added
!     these modifications were done June 22 2009 

implicit none

real(4), intent(in)   ::   p,t,pv,freq
real(4), intent(out)  ::   absoxy


integer(4), parameter        ::   nlines=44
integer(4), save             ::   istart=1
integer(4)                   ::   i
real(4)                      ::   gamoxy
real(4)                      ::   tht,pwet,pdry,ga,gasq,delta,rnuneg,rnupos,ff,zterm,apterm,sftot,xterm
real(4)                      ::   ga1, ga2
real(4), dimension(6,nlines), save ::   h
real(4), dimension(nlines), save   ::   f0, a1, a2, a3, a4, a5, a6
real(8)                      ::   sum
      
real(4), parameter            ::  db_to_neper = 0.2302585094 !convert db/km to neper/km
     
real(4), parameter            ::  btrans=2.5, ctrans=3.5
real(4)                       ::  xtrans

if (istart == 1) then

    istart = 0
    
    a4( 1:38) = 0.0
    a4(39:44) = 0.6

!          freq          a1      a2       a3       a5          a6
    h = reshape ((/                                             &  
       50.474238,    0.94e-6,  9.694,  8.60e-3,  0.210,  0.685, &
       50.987749,    2.46e-6,  8.694,  8.70e-3,  0.190,  0.680, &
       51.503350,    6.08e-6,  7.744,  8.90e-3,  0.171,  0.673, &
       52.021410,   14.14e-6,  6.844,  9.20e-3,  0.144,  0.664, &
       52.542394,   31.02e-6,  6.004,  9.40e-3,  0.118,  0.653, &
       53.066907,   64.10e-6,  5.224,  9.70e-3,  0.114,  0.621, &
       53.595749,  124.70e-6,  4.484, 10.00e-3,  0.200,  0.508, &
       54.130000,  228.00e-6,  3.814, 10.20e-3,  0.291,  0.375, &
       54.671159,  391.80e-6,  3.194, 10.50e-3,  0.325,  0.265, &
       55.221367,  631.60e-6,  2.624, 10.79e-3,  0.224,  0.295, &
       55.783802,  953.50e-6,  2.119, 11.10e-3, -0.144,  0.613, &
       56.264775,  548.90e-6,  0.015, 16.46e-3,  0.339, -0.098, &
       56.363389, 1344.00e-6,  1.660, 11.44e-3, -0.258,  0.655, &
       56.968206, 1763.00e-6,  1.260, 11.81e-3, -0.362,  0.645, &
       57.612484, 2141.00e-6,  0.915, 12.21e-3, -0.533,  0.606, &
       58.323877, 2386.00e-6,  0.626, 12.66e-3, -0.178,  0.044, &
       58.446590, 1457.00e-6,  0.084, 14.49e-3,  0.650, -0.127, &
       59.164207, 2404.00e-6,  0.391, 13.19e-3, -0.628,  0.231, &
       59.590983, 2112.00e-6,  0.212, 13.60e-3,  0.665, -0.078, &
       60.306061, 2124.00e-6,  0.212, 13.82e-3, -0.613,  0.070, &
       60.434776, 2461.00e-6,  0.391, 12.97e-3,  0.606, -0.282, &
       61.150560, 2504.00e-6,  0.626, 12.48e-3,  0.090, -0.058, &
       61.800154, 2298.00e-6,  0.915, 12.07e-3,  0.496, -0.662, &
       62.411215, 1933.00e-6,  1.260, 11.71e-3,  0.313, -0.676, &
       62.486260, 1517.00e-6,  0.083, 14.68e-3, -0.433,  0.084, &
       62.997977, 1503.00e-6,  1.665, 11.39e-3,  0.208, -0.668, &
       63.568518, 1087.00e-6,  2.115, 11.08e-3,  0.094, -0.614, &
       64.127767,  733.50e-6,  2.620, 10.78e-3, -0.270, -0.289, &
       64.678903,  463.50e-6,  3.195, 10.50e-3, -0.366, -0.259, &
       65.224071,  274.80e-6,  3.815, 10.20e-3, -0.326, -0.368, &
       65.764772,  153.00e-6,  4.485, 10.00e-3, -0.232, -0.500, &
       66.302091,   80.09e-6,  5.225,  9.70e-3, -0.146, -0.609, &
       66.836830,   39.46e-6,  6.005,  9.40e-3, -0.147, -0.639, &
       67.369598,   18.32e-6,  6.845,  9.20e-3, -0.174, -0.647, &
       67.900867,    8.01e-6,  7.745,  8.90e-3, -0.198, -0.655, &
       68.431005,    3.30e-6,  8.695,  8.70e-3, -0.210, -0.660, &
       68.960311,    1.28e-6,  9.695,  8.60e-3, -0.220, -0.665, &
      118.750343,  945.00e-6,  0.009, 16.30e-3, -0.031,  0.008, &
      368.498350,   67.90e-6,  0.049, 19.20e-3,  0.0,    0.0,   &
      424.763124,  638.00e-6,  0.044, 19.16e-3,  0.0,    0.0,   & 
      487.249370,  235.00e-6,  0.049, 19.20e-3,  0.0,    0.0,   &
      715.393150,   99.60e-6,  0.145, 18.10e-3,  0.0,    0.0,   &
      773.839675,  671.00e-6,  0.130, 18.10e-3,  0.0,    0.0,   &
      834.145330,  180.00e-6,  0.147, 18.10e-3,  0.0,    0.0    /), (/6,nlines/))
      
      f0(:)=h(1,:)
      a1(:)=h(2,:)/h(1,:)
      a2(:)=h(3,:)
      a3(:)=h(4,:)
      a5(:)=0.001*h(5,:)
      a6(:)=0.001*h(6,:)
 
endif ! istart=1

 
tht = 300.0/t
pwet=0.1*pv
pdry=0.1*p-pwet
xterm=1.0 - tht

sum = 0.
do i=1,nlines
      ga = a3(i)*(pdry*tht**(0.8-a4(i)) + 1.1*tht*pwet)
      gasq=ga*ga
	  delta=(a5(i) + a6(i)*tht)*p*tht**0.8
      rnuneg = f0(i)-freq
      rnupos = f0(i)+freq
      ff = (ga-rnuneg*delta)/(gasq+rnuneg**2) +  (ga-rnupos*delta)/(gasq+rnupos**2)
      sum = sum + ff*a1(i)*exp(a2(i)*xterm)
enddo
if(sum.lt.0) sum=0.0

!     add nonresonant contribution

ga1=5.6e-3*(pdry+1.1*pwet)*1.097687*(tht**0.8)  
! Liebe temperature dependence at L-band and below   
! 1.097687=(300/267)**0.8: overall normalization 

ga2=5.6e-3*(pdry+1.1*pwet)*(tht**1.5)  
! Wentz Meissner 2016 at C-band and above
! transition from L to C bands

xtrans = 1.0 - 1.0/(1.0 + exp(-(btrans*(freq-ctrans))))
ga = xtrans*ga1 + (1.0-xtrans)*ga2
      
zterm=ga*(1.+(freq/ga)**2)
apterm=1.4e-10*(1-1.2e-5*freq**1.5)*pdry*tht**1.5
if(apterm.lt.0.0) apterm=0.0

sftot=pdry*freq*(tht**2) * (tht*sum + 6.14e-4/zterm + apterm)

gamoxy=0.1820*freq*sftot
!x    if(freq.gt.37) gamoxy=gamoxy + 0.1820*43.e-10 *pdry**2*tht**3*(freq-37.)**1.7  !prior to 7/17/2015
if(freq.gt.37) gamoxy=gamoxy + 0.1820*26.e-10 *pdry**2*tht**3*(freq-37.)**1.8  !implemented 7/17/2015.

absoxy = db_to_neper * gamoxy

return
end subroutine absoxy_RSS_2022


end module RSS_2022_ABSORPTION