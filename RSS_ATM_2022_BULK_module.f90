module RSS_ATM_2022_BULK_module
public
save

public :: read_freq_table
public :: read_ATM_model_table
public :: get_RSS_ATM_BULK
private:: find_permittivity_meissner_wentz, dielectric_meissner_wentz

integer(4), parameter           :: nfreq_max = 50 ! max number of bands to be processed
character(len=250), parameter   :: ftab_file = 'bulk_ATM_frequency_band_table.txt'
character(len=250), parameter   :: tab_path =  'bulk_tables\'

real(4), parameter              :: dskin=0.3  ! difference between skin (ERA5) and subskin (MW sensor) SST
real(4), parameter              :: c=29.979   ! velocity of light. Fixed bug. See Richard's email 
real(4), parameter              :: pi=3.14159

real(4), parameter              :: rtol=0.001, eps=0.0001

integer(4), parameter           :: iunit=3

integer(4)                      :: nfreq_tab
real(4), dimension(nfreq_max)   :: freq_arr

integer(4)                      :: nt_bin  
integer(4)                      :: nv_bin
integer(4)                      :: ntv_bin
integer(4)                      :: nto_bin

real(4)                         :: dt_bin,  t_min  
real(4)                         :: dv_bin,  v_min
real(4)                         :: dtv_bin, tv_min
real(4)                         :: dto_bin, to_min

real(8), allocatable, dimension(:,:)    :: arr1_TV, arr1_TO, arr1_TL, arr1_TU, arr1_TD
real(8), allocatable, dimension(:,:,:)  :: arr2_TV, arr2_TO, arr2_TL, arr2_TU, arr2_TD

real(8), allocatable, dimension(:,:,:)  :: arr_AV
real(8), allocatable, dimension(:,:)    :: arr_AO
real(8), allocatable, dimension(:,:)    :: arr_AL

real(8), allocatable, dimension(:,:)    :: acoef_AV

contains


subroutine get_RSS_ATM_BULK (IBAND, SURTEP, V, L, THT,               &
                             F,                                      &                                                           
                             TEFF_O, TEFF_V, TEFF_L, TEFF_U, TEFF_D, &
                             AO, AV, AL, OPACITY, TRAN, TBUP, TBDW)

implicit none

integer(4), intent(in)       :: IBAND 
!frequency band index see bulk_ATM_frequency_band_table.txt'

real(4), intent(in)          :: SURTEP
! subskin sea surface temperature. measured by microwave satellite (Kelvin)

real(4), intent(in)          :: V
! columnar water vapor (mm)

real(4), intent(in)          :: L
! columnar liquid cloud water (mm)

real(4), intent(in), optional :: THT
! Earth Incidence Angle (deg)

real(4), intent(out), optional :: F
! center frequency of band (GHz), for consistency check

real(4), intent(out), optional :: TEFF_O
! effective average atmospheric temperature of air (Kelvin)

real(4), intent(out), optional :: TEFF_V
! effective average atmospheric temperature of water vapor (Kelvin) 

real(4), intent(out), optional :: TEFF_L
! effective average atmospheric temperature of liquid cloud water vapor (Kelvin) 

real(4), intent(out), optional :: TEFF_U
! effective atmospheric up-welling temperature (Kelvin) 

real(4), intent(out), optional :: TEFF_D
! effective atmospheric down-welling temperature (Kelvin) 

real(4), intent(out), optional :: AO
! total columnar oxygen absorption (neper) 

real(4), intent(out), optional :: AV
! total columnar water vapor absorption (neper) 

real(4), intent(out), optional :: AL
! total columnar liquid cloud water absorption (neper) 

real(4), intent(out), optional :: OPACITY
! total columnar atmospheric absorption = opacity (neper)

real(4), intent(out), optional :: TRAN
! atmospheric transmittance (dimensionless)  

real(4), intent(out), optional :: TBUP
! up-welling brightness temperature (Kelvin)

real(4), intent(out), optional :: TBDW
! down-welling brightness temperature (Kelvin)


real(4)    :: xtht, ts, xf, wavlen, costht, xpath
real(4)    :: xteff_o, xteff_v, xteff_l, xteff_u, xteff_d, xav, xao, xal, xopacity, xtran, xtbup, xtbdw
real(4)    :: xts, xv, xtv, xto
real(4)    :: t0, t1, t_low, brief_t
real(4)    :: v0, v1, v_low, brief_v
real(4)    :: tv0, tv1, tv_low, brief_tv
real(4)    :: to0, to1, to_low, brief_to
integer(4) :: it1, it2, iv1, iv2, itv1, itv2, ito1, ito2
real(4)    :: xarr11, xarr12, xarr21, xarr22, xarr1, xarr2, xarr
real(4)    :: zarr1, zarr2, zarr

real(4)    :: ctemp
complex(4) :: permit


integer(4), save :: istart=1, iwarning=1
integer(4)       :: istop

if (iwarning==1 .and. L>0.25) then
    iwarning=0
    write(*,*) 'L= ',L
    write(*,*) 'RSS ATM BULK only set up for non-raining atmosphere and requires L<0.25.'
    write(*,*) 'high cloud water in BULK ATM. PGM STOP.'
    write(*,*) 'enter 0 to continue.'
    write(*,*) 'enter 1 to stop.'
    read(*,*)  istop
    if (istop==1) stop
endif

if (istart==1) then
    istart=0
    call read_freq_table
    call read_ATM_model_table
endif 

if (iband > nfreq_max) then
    write(*,*) iband, nfreq_max
    write(*,*) 'invalid frequency band index in get_RSS_ATM_BULK'
    stop
endif

XF = FREQ_ARR(IBAND)

TS = SURTEP - DSKIN
! The tables were derived from ERA5. The ERA5 skin TS is about 0.3 K colder than satellite SST

! TS bin
t0 = t_min+dt_bin/2.0
t1 = t0 + (nt_bin-1)*dt_bin 
xts = ts
if (xts .le. t0+eps)  xts=t0+eps
if (xts .ge. t1-eps)  xts=t1-eps
it1 = 1 + floor((xts - t0)/dt_bin)
it2 = it1+1
if (it1 .lt. 1  .or. it2 .gt. nt_bin) then   
    write(*,*) ts,xts,t0,t1,it1,it2
    write(*,*) 'sstbin error in get_RSS_ATM_BULK '
    stop
endif 
t_low = t0 + (it1-1)*dt_bin
brief_t = (xts - t_low) / dt_bin


! V bin
v0 = v_min+dv_bin/2.0
v1 = v0 + (nv_bin-1)*dv_bin 
xv = v
if (xv .le. v0+eps)  xv=v0+eps
if (xv .ge. v1-eps)  xv=v1-eps
iv1 = 1 + floor((xv - v0)/dv_bin)
iv2 = iv1+1
if (iv1 .lt. 1  .or. iv2 .gt. nv_bin) then   
    write(*,*) v,xv,v0,v1,iv1,iv2
    write(*,*) 'vbin error in get_RSS_ATM_BULK '
    stop
endif 
v_low = v0 + (iv1-1)*dv_bin
brief_v = (xv - v_low) / dv_bin

! TO
ZARR1  = arr1_TO(iband,iv1)
ZARR2  = arr1_TO(iband,iv2)
ZARR   = (1.0 - brief_v)*ZARR1 + brief_v*ZARR2
XARR11 = arr2_TO(iband,it1,iv1)
XARR12 = arr2_TO(iband,it1,iv2)
XARR21 = arr2_TO(iband,it2,iv1)
XARR22 = arr2_TO(iband,it2,iv2)
XARR1 = (1.0 - brief_v)*XARR11 + brief_v*XARR12
XARR2 = (1.0 - brief_v)*XARR21 + brief_v*XARR22
XARR  = (1.0 - brief_t)*XARR1  + brief_t*XARR2
XTEFF_O = ZARR + XARR

! TV
ZARR1  = arr1_TV(iband,iv1)
ZARR2  = arr1_TV(iband,iv2)
ZARR   = (1.0 - brief_v)*ZARR1 + brief_v*ZARR2
XARR11 = arr2_TV(iband,it1,iv1)
XARR12 = arr2_TV(iband,it1,iv2)
XARR21 = arr2_TV(iband,it2,iv1)
XARR22 = arr2_TV(iband,it2,iv2)
XARR1 = (1.0 - brief_v)*XARR11 + brief_v*XARR12
XARR2 = (1.0 - brief_v)*XARR21 + brief_v*XARR22
XARR  = (1.0 - brief_t)*XARR1  + brief_t*XARR2
XTEFF_V = ZARR + XARR

! TL
ZARR1  = arr1_TL(iband,iv1)
ZARR2  = arr1_TL(iband,iv2)
ZARR   = (1.0 - brief_v)*ZARR1 + brief_v*ZARR2
XARR11 = arr2_TL(iband,it1,iv1)
XARR12 = arr2_TL(iband,it1,iv2)
XARR21 = arr2_TL(iband,it2,iv1)
XARR22 = arr2_TL(iband,it2,iv2)
XARR1 = (1.0 - brief_v)*XARR11 + brief_v*XARR12
XARR2 = (1.0 - brief_v)*XARR21 + brief_v*XARR22
XARR  = (1.0 - brief_t)*XARR1  + brief_t*XARR2
XTEFF_L = ZARR + XARR

! TU
ZARR1  = arr1_TU(iband,iv1)
ZARR2  = arr1_TU(iband,iv2)
ZARR   = (1.0 - brief_v)*ZARR1 + brief_v*ZARR2
XARR11 = arr2_TU(iband,it1,iv1)
XARR12 = arr2_TU(iband,it1,iv2)
XARR21 = arr2_TU(iband,it2,iv1)
XARR22 = arr2_TU(iband,it2,iv2)
XARR1 = (1.0 - brief_v)*XARR11 + brief_v*XARR12
XARR2 = (1.0 - brief_v)*XARR21 + brief_v*XARR22
XARR  = (1.0 - brief_t)*XARR1  + brief_t*XARR2
XTEFF_U = ZARR + XARR

! TD
ZARR1  = arr1_TD(iband,iv1)
ZARR2  = arr1_TD(iband,iv2)
ZARR   = (1.0 - brief_v)*ZARR1 + brief_v*ZARR2
XARR11 = arr2_TD(iband,it1,iv1)
XARR12 = arr2_TD(iband,it1,iv2)
XARR21 = arr2_TD(iband,it2,iv1)
XARR22 = arr2_TD(iband,it2,iv2)
XARR1 = (1.0 - brief_v)*XARR11 + brief_v*XARR12
XARR2 = (1.0 - brief_v)*XARR21 + brief_v*XARR22
XARR  = (1.0 - brief_t)*XARR1  + brief_t*XARR2
XTEFF_D = ZARR + XARR

! TV bin
tv0 = tv_min+dtv_bin/2.0
tv1 = tv0 + (ntv_bin-1)*dtv_bin
xtv = xteff_v
if (xtv .le. tv0+eps) xtv=tv0+eps
if (xtv .ge. tv1-eps) xtv=tv1-eps
itv1 = 1+ floor((xtv - tv0)/dtv_bin)
itv2 = itv1+1
if (itv1 .lt. 1 .or. itv2 .gt. ntv_bin) then 
  write(*,*) xteff_v,xtv,tv0,tv1,itv1,itv2
  write(*,*) 'tvbin error in get_RSS_ATM_BULK '
  stop
endif
tv_low = tv0 + (itv1-1)*dtv_bin
brief_tv = (xtv - tv_low) / dtv_bin

! AV
XARR11 = arr_AV(iband,itv1,iv1)
XARR12 = arr_AV(iband,itv1,iv2)
XARR21 = arr_AV(iband,itv2,iv1)
XARR22 = arr_AV(iband,itv2,iv2)
XARR1 = (1.0 - brief_v)*XARR11  + brief_v*XARR12
XARR2 = (1.0 - brief_v)*XARR21  + brief_v*XARR22
XARR  = (1.0 - brief_tv)*XARR1  + brief_tv*XARR2
XAV   = (acoef_AV(iband,1)*V + acoef_AV(iband,2)*V*V) * XARR

! TO bin
! extrapolate if outside interval boundary
to0 = to_min+dto_bin/2.0
to1 = to0 + (nto_bin-1)*dto_bin
xto = xteff_o
ito1 = 1 + floor((xto - to0)/dto_bin)
if (ito1 .lt. 1)         ito1 = 1
if (ito1 .ge. nto_bin-1) ito1 = nto_bin-1
ito2 = ito1+1
to_low = to0 + (ito1-1)*dto_bin
brief_to = (xto - to_low) / dto_bin

! AO
XARR1 = arr_AO(iband,ito1)
XARR2 = arr_AO(iband,ito2)
XAO = (1.0 - brief_to)*XARR1 + brief_to*XARR2

! AL
wavlen = c/xf
ctemp = xteff_l ! the FORTRAN MW routine requires temperature input in Kelvin
call find_permittivity_meissner_wentz(xf,ctemp,0.0,   permit) 
ZARR1  = arr_AL(iband,it1)
ZARR2  = arr_AL(iband,it2)
ZARR   = (1.0 - brief_t)*ZARR1 + brief_t*ZARR2
xal = 0.1*(6.0*pi/wavlen)*L*aimag((1.0-permit)/(2.0+permit)) + ZARR*L 


xopacity = xao+xav+xal
if (xopacity .lt. 0.0) xopacity=0.0

if (present(f))         f=xf

if (present(teff_o))    teff_o=xteff_o  
if (present(teff_v))    teff_v=xteff_v
if (present(teff_l))    teff_l=xteff_l
if (present(teff_u))    teff_u=xteff_u
if (present(teff_d))    teff_d=xteff_d

if (present(ao))        ao=xao  
if (present(av))        av=xav  
if (present(al))        al=xal  
if (present(opacity))   opacity=xopacity


if (present(tran) .or. present(tbup) .or. present(tbdw)) then

    if (.not.(present(tht))) then
        write(*,*) 'input error in get_RSS_ATM_BULK'
        write(*,*) 'computation of tran, tbup, tbdw requires tht input.'
        stop
    endif

    xtht=tht
    ! AMSR ATBD equation (15)
    ! (1+hratio)/sqrt(costht**2+hratio*(2+hratio)), hratio=.00035
    costht = COSD(XTHT)
    xpath=1.00035/sqrt(costht*costht+7.001225e-4)   

    XTRAN = exp(-XOPACITY*XPATH)
    XTBUP = (1.0-XTRAN)*XTEFF_U
    XTBDW = (1.0-XTRAN)*XTEFF_D 
    
    if (present(tran)) tran=xtran   
    if (present(tbup)) tbup=xtbup   
    if (present(tbdw)) tbdw=xtbdw   

endif  
  
return
end subroutine get_RSS_ATM_BULK


subroutine read_freq_table
implicit none

character(len=150) :: string_dummy
integer(4)         :: ifreq, jfreq
logical(4)         :: lexist

inquire(file=ftab_file, exist=lexist)
if (.not.(lexist)) then
    write(*,*) ftab_file
    write(*,*) 'file does not exist'
    stop
endif

open(unit=iunit,form='formatted',file=ftab_file,status='old',action='read')
read(iunit,*) nfreq_tab
read(iunit,*) string_dummy
do ifreq=1,nfreq_tab
    read(iunit,*,end=199) jfreq, freq_arr(jfreq)
    if (ifreq /= jfreq) then
        write(*,*) ifreq,jfreq,'index mismatch in read_freq_table. pgm stop.'
        stop 
    endif
enddo
close(iunit)
return

199 continue
write(*,*) ifreq, 'error in read_freq_table. pgm stop.'
close(3)
stop

end subroutine read_freq_table



subroutine read_ATM_model_table
implicit none

character(len=250) :: tab_file
character(len=100) :: str_file
logical(4)         :: lexist

integer(4)         :: mfreq_max, mfreq_tab
real(4)            :: dvv, vvmin 
real(4)            :: dtt, ttmin
integer(4)         :: mtt, mvv 

! TO
str_file = 'table_bulk_TO.dat'
tab_file = trim(tab_path)//trim(str_file)
inquire(file=tab_file, exist=lexist)
if (.not.(lexist)) then
    write(*,*) tab_file
    write(*,*) 'file does not exist'
    stop
endif
open(unit=iunit,file=tab_file,action='read',form='unformatted', access='stream',status='old')  
read(iunit) mfreq_max, mfreq_tab
if (mfreq_max/=nfreq_max .or. mfreq_tab/=nfreq_tab) then
    write(*,*) mfreq_max,nfreq_max,mfreq_tab,nfreq_tab
    write(*,*) tab_file
    write(*,*) 'freq mismatch'
    stop
endif
read(iunit) mtt, dtt, ttmin
read(iunit) mvv, dvv, vvmin
nt_bin = mtt
nv_bin = mvv
dt_bin = dtt
dv_bin = dvv
t_min  = ttmin
v_min  = vvmin
if (.not.allocated(arr1_TO))   allocate(arr1_TO(nfreq_max,nv_bin))
if (.not.allocated(arr2_TO))   allocate(arr2_TO(nfreq_max,nt_bin,nv_bin))
read(iunit) arr1_TO
read(iunit) arr2_TO
close(iunit)

! TV
str_file = 'table_bulk_TV.dat'
tab_file = trim(tab_path)//trim(str_file)
inquire(file=tab_file, exist=lexist)
if (.not.(lexist)) then
    write(*,*) tab_file
    write(*,*) 'file does not exist'
    stop
endif
open(unit=iunit,file=tab_file,action='read',form='unformatted', access='stream',status='old')  
read(iunit) mfreq_max, mfreq_tab
if (mfreq_max/=nfreq_max .or. mfreq_tab/=nfreq_tab) then
    write(*,*) mfreq_max,nfreq_max,mfreq_tab,nfreq_tab
    write(*,*) tab_file
    write(*,*) 'freq mismatch'
    stop
endif
read(iunit) mtt, dtt, ttmin
read(iunit) mvv, dvv, vvmin
if (nt_bin /= mtt .or. nv_bin /= mvv .or. &
    abs(dtt-dt_bin)>rtol .or. abs(dvv-dv_bin)>rtol .or. abs(t_min-ttmin)>rtol .or. abs(v_min-vvmin)>rtol) then
    write(*,*) tab_file
    write(*,*) mtt, dtt, ttmin, nt_bin, dt_bin, t_min
    write(*,*) mvv, dvv, vvmin, nv_bin, dv_bin, v_min
    write(*,*) 'input mismatch'
    stop
endif
dt_bin = dtt
dv_bin = dvv
t_min  = ttmin
v_min  = vvmin
if (.not.allocated(arr1_TV))   allocate(arr1_TV(nfreq_max,nv_bin))
if (.not.allocated(arr2_TV))   allocate(arr2_TV(nfreq_max,nt_bin,nv_bin))
read(iunit) arr1_TV
read(iunit) arr2_TV
close(iunit)

! TL
str_file = 'table_bulk_TL.dat'
tab_file = trim(tab_path)//trim(str_file)
inquire(file=tab_file, exist=lexist)
if (.not.(lexist)) then
    write(*,*) tab_file
    write(*,*) 'file does not exist'
    stop
endif
open(unit=iunit,file=tab_file,action='read',form='unformatted', access='stream',status='old')  
read(iunit) mfreq_max, mfreq_tab
if (mfreq_max/=nfreq_max .or. mfreq_tab/=nfreq_tab) then
    write(*,*) mfreq_max,nfreq_max,mfreq_tab,nfreq_tab
    write(*,*) tab_file
    write(*,*) 'freq mismatch'
    stop
endif
read(iunit) mtt, dtt, ttmin
read(iunit) mvv, dvv, vvmin
if (nt_bin /= mtt .or. nv_bin /= mvv .or. &
    abs(dtt-dt_bin)>rtol .or. abs(dvv-dv_bin)>rtol .or. abs(t_min-ttmin)>rtol .or. abs(v_min-vvmin)>rtol) then
    write(*,*) tab_file
    write(*,*) mtt, dtt, ttmin, nt_bin, dt_bin, t_min
    write(*,*) mvv, dvv, vvmin, nv_bin, dv_bin, v_min
    write(*,*) 'input mismatch'
    stop
endif
dt_bin = dtt
dv_bin = dvv
t_min  = ttmin
v_min  = vvmin
if (.not.allocated(arr1_TL))   allocate(arr1_TL(nfreq_max,nv_bin))
if (.not.allocated(arr2_TL))   allocate(arr2_TL(nfreq_max,nt_bin,nv_bin))
read(iunit) arr1_TL
read(iunit) arr2_TL
close(iunit)

! TU
str_file = 'table_bulk_TU.dat'
tab_file = trim(tab_path)//trim(str_file)
inquire(file=tab_file, exist=lexist)
if (.not.(lexist)) then
    write(*,*) tab_file
    write(*,*) 'file does not exist'
    stop
endif
open(unit=iunit,file=tab_file,action='read',form='unformatted', access='stream',status='old')  
read(iunit) mfreq_max, mfreq_tab
if (mfreq_max/=nfreq_max .or. mfreq_tab/=nfreq_tab) then
    write(*,*) mfreq_max,nfreq_max,mfreq_tab,nfreq_tab
    write(*,*) tab_file
    write(*,*) 'freq mismatch'
    stop
endif
read(iunit) mtt, dtt, ttmin
read(iunit) mvv, dvv, vvmin
if (nt_bin /= mtt .or. nv_bin /= mvv .or. &
    abs(dtt-dt_bin)>rtol .or. abs(dvv-dv_bin)>rtol .or. abs(t_min-ttmin)>rtol .or. abs(v_min-vvmin)>rtol) then
    write(*,*) tab_file
    write(*,*) mtt, dtt, ttmin, nt_bin, dt_bin, t_min 
    write(*,*) mvv, dvv, vvmin, nv_bin, dv_bin, v_min 
    write(*,*) 'input mismatch'
    stop
endif
dt_bin = dtt
dv_bin = dvv
t_min  = ttmin
v_min  = vvmin
if (.not.allocated(arr1_TU))   allocate(arr1_TU(nfreq_max,nv_bin))
if (.not.allocated(arr2_TU))   allocate(arr2_TU(nfreq_max,nt_bin,nv_bin))
read(iunit) arr1_TU
read(iunit) arr2_TU
close(iunit)

! TD
str_file = 'table_bulk_TD.dat'
tab_file = trim(tab_path)//trim(str_file)
inquire(file=tab_file, exist=lexist)
if (.not.(lexist)) then
    write(*,*) tab_file
    write(*,*) 'file does not exist'
    stop
endif
open(unit=iunit,file=tab_file,action='read',form='unformatted', access='stream',status='old')  
read(iunit) mfreq_max, mfreq_tab
if (mfreq_max/=nfreq_max .or. mfreq_tab/=nfreq_tab) then
    write(*,*) mfreq_max,nfreq_max,mfreq_tab,nfreq_tab
    write(*,*) tab_file
    write(*,*) 'freq mismatch'
    stop
endif
read(iunit) mtt, dtt, ttmin
read(iunit) mvv, dvv, vvmin
if (nt_bin /= mtt .or. nv_bin /= mvv .or. &
    abs(dtt-dt_bin)>rtol .or. abs(dvv-dv_bin)>rtol .or. abs(t_min-ttmin)>rtol .or. abs(v_min-vvmin)>rtol) then
    write(*,*) tab_file
    write(*,*) mtt, dtt, ttmin, nt_bin, dt_bin, t_min 
    write(*,*) mvv, dvv, vvmin, nv_bin, dv_bin, v_min 
    write(*,*) 'input mismatch'
    stop
endif
dt_bin = dtt
dv_bin = dvv
t_min  = ttmin
v_min  = vvmin
if (.not.allocated(arr1_TD))   allocate(arr1_TD(nfreq_max,nv_bin))
if (.not.allocated(arr2_TD))   allocate(arr2_TD(nfreq_max,nt_bin,nv_bin))
read(iunit) arr1_TD
read(iunit) arr2_TD
close(iunit)

! AO
str_file = 'table_bulk_AO.dat'
tab_file = trim(tab_path)//trim(str_file)
inquire(file=tab_file, exist=lexist)
if (.not.(lexist)) then
    write(*,*) tab_file
    write(*,*) 'file does not exist'
    stop
endif
open(unit=iunit,file=tab_file,action='read',form='unformatted', access='stream',status='old')  
read(iunit) mfreq_max, mfreq_tab
if (mfreq_max/=nfreq_max .or. mfreq_tab/=nfreq_tab) then
    write(*,*) mfreq_max,nfreq_max,mfreq_tab,nfreq_tab
    write(*,*) tab_file
    write(*,*) 'freq mismatch'
    stop
endif
read(iunit) mtt, dtt, ttmin
nto_bin = mtt
dto_bin = dtt
to_min  = ttmin
if (.not.allocated(arr_AO))   allocate(arr_AO(nfreq_max,nto_bin))
read(iunit) arr_AO
close(iunit)

! AV
str_file = 'table_bulk_AV.dat'
tab_file = trim(tab_path)//trim(str_file)
inquire(file=tab_file, exist=lexist)
if (.not.(lexist)) then
    write(*,*) tab_file
    write(*,*) 'file does not exist'
    stop
endif
open(unit=iunit,file=tab_file,action='read',form='unformatted', access='stream',status='old')  
read(iunit) mfreq_max, mfreq_tab
if (mfreq_max/=nfreq_max .or. mfreq_tab/=nfreq_tab) then
    write(*,*) mfreq_max,nfreq_max,mfreq_tab,nfreq_tab
    write(*,*) tab_file
    write(*,*) 'freq mismatch'
    stop
endif
read(iunit) mtt, dtt, ttmin
read(iunit) mvv, dvv, vvmin
if (nv_bin /= mvv .or.  abs(dvv-dv_bin)>rtol .or. abs(v_min-vvmin)>rtol  ) then
    write(*,*) tab_file
    write(*,*) mvv, dvv, vvmin, nv_bin, dv_bin, v_min 
    write(*,*) 'input mismatch'
    stop
endif
ntv_bin = mtt
dtv_bin = dtt
dv_bin  = dvv
tv_min  = ttmin
v_min   = vvmin
if (.not.allocated(arr_AV))   allocate(arr_AV(nfreq_max,ntv_bin,nv_bin))
if (.not.allocated(acoef_AV)) allocate(acoef_AV(nfreq_max,2))
read(iunit) acoef_AV
read(iunit) arr_AV
close(iunit)

! AL
str_file = 'table_bulk_AL.dat'
tab_file = trim(tab_path)//trim(str_file)
inquire(file=tab_file, exist=lexist)
if (.not.(lexist)) then
    write(*,*) tab_file
    write(*,*) 'file does not exist'
    stop
endif
open(unit=iunit,file=tab_file,action='read',form='unformatted', access='stream',status='old')  
read(iunit) mfreq_max, mfreq_tab
if (mfreq_max/=nfreq_max .or. mfreq_tab/=nfreq_tab) then
    write(*,*) mfreq_max,nfreq_max,mfreq_tab,nfreq_tab
    write(*,*) tab_file
    write(*,*) 'freq mismatch'
    stop
endif
read(iunit) mtt, dtt, ttmin
if (nt_bin /= mtt .or. abs(dtt-dt_bin)>rtol .or. abs(t_min-ttmin)>rtol ) then
    write(*,*) tab_file
    write(*,*) mtt, dtt, ttmin, nt_bin, dt_bin, t_min 
    write(*,*) 'input mismatch'
    stop
endif
nt_bin = mtt
dt_bin = dtt
t_min  = ttmin
if (.not.allocated(arr_AL))   allocate(arr_AL(nfreq_max,nt_bin))
read(iunit) arr_AL
close(iunit)


return
end subroutine read_ATM_model_table

!   input:
!   name   parameter  unit  range
!   freq   frequency  [ghz] 1 to 400
!   t      sst        [k]   248.16 k (-25 c) to 313.16 k (40 c) for pure water
!                           271.16 k (-2  c) to 307.16 k (34 c) for saline water
!   s      salinity   [ppt]  0 to 40
!
!   output:
!   eps    complex dielectric constant 
!          negative imaginary part to be consistent with wentz1 conventionc
!
!
subroutine	find_permittivity_meissner_wentz(freq,t,s,   eps)
implicit none
    
	real(4), parameter :: f0=17.97510
    real(4) freq,t,sst,s
	real(4) e0s,e1s,e2s,n1s,n2s,sig
    complex(4) :: j = (0.0,1.0), eps

	sst = t - 273.15 ! [celsius]
    call dielectric_meissner_wentz(sst,s,  e0s,e1s,e2s,n1s,n2s,sig)

!     debye law (2 relaxation wavelengths)
    eps = (e0s - e1s)/(1.0 - j*(freq/n1s)) + (e1s - e2s)/(1.0 - j*(freq/n2s)) + e2s +  j*sig*f0/freq
	eps = conjg(eps)

return 
end	subroutine find_permittivity_meissner_wentz    


subroutine dielectric_meissner_wentz(sst_in,s,   e0s,e1s,e2s,n1s,n2s,sig)
!
!     complex dielectric constant: eps
!     [MW 2004, MW 2012].
!     
!     Changes from [MW 2012]:
!     1. Typo (sign) in the printed version of coefficient d3 in Table 7. Its value should be -0.35594E-06.
!     2. Changed SST behavior of coefficient b2 from:
!     b2 = 1.0 + s*(z(10) + z(11)*sst) to
!     b2 = 1.0 + s*(z(10) + 0.5*z(11)*(sst + 30)) 
!
!!
!     input:
!     name   parameter  unit  range
!     sst      sst        [c]   -25 c to 40 c for pure water
!                               -2  c to 34 c for saline water
!     s      salinity   [ppt]  0 to 40
!
!     output:
!     eps    complex dielectric constant
!            negative imaginary part to be consistent with wentz1 convention
!

implicit none


    real(4), intent(in)  :: sst_in,s
    real(4), intent(out) :: e0s,e1s,e2s,n1s,n2s,sig
 
    real(4), dimension(11), parameter :: &
      x=(/ 5.7230e+00, 2.2379e-02, -7.1237e-04, 5.0478e+00, -7.0315e-02, 6.0059e-04, 3.6143e+00, &
           2.8841e-02, 1.3652e-01,  1.4825e-03, 2.4166e-04 /)
    
    real(4), dimension(13), parameter :: &
      z=(/ -3.56417e-03,  4.74868e-06,  1.15574e-05,  2.39357e-03, -3.13530e-05, &
            2.52477e-07, -6.28908e-03,  1.76032e-04, -9.22144e-05, -1.99723e-02, &
            1.81176e-04, -2.04265e-03,  1.57883e-04  /)  ! 2004

    real(4), dimension(3), parameter :: a0coef=(/ -0.33330E-02,  4.74868e-06,  0.0e+00/)
    real(4), dimension(5), parameter :: b1coef=(/0.23232E-02, -0.79208E-04, 0.36764E-05, -0.35594E-06, 0.89795E-08/)
 
    real(4) :: e0,e1,e2,n1,n2
    real(4) :: a0,a1,a2,b1,b2
    real(4) :: sig35,r15,rtr15,alpha0,alpha1

    real(4) :: sst,sst2,sst3,sst4,s2
    
    sst=sst_in
    if(sst.lt.-30.16) sst=-30.16  !protects against n1 and n2 going zero for very cold water
    
    sst2=sst*sst
    sst3=sst2*sst
    sst4=sst3*sst

    s2=s*s
 
    !     pure water
    e0    = (3.70886e4 - 8.2168e1*sst)/(4.21854e2 + sst) ! stogryn et al.
    e1    = x(1) + x(2)*sst + x(3)*sst2
    n1    = (45.00 + sst)/(x(4) + x(5)*sst + x(6)*sst2)
    e2    = x(7) + x(8)*sst
    n2    = (45.00 + sst)/(x(9) + x(10)*sst + x(11)*sst2)
    
    !     saline water
    !     conductivity [s/m] taken from stogryn et al.
    sig35 = 2.903602 + 8.60700e-2*sst + 4.738817e-4*sst2 - 2.9910e-6*sst3 + 4.3047e-9*sst4
    r15   = s*(37.5109+5.45216*s+1.4409e-2*s2)/(1004.75+182.283*s+s2)
    alpha0 = (6.9431+3.2841*s-9.9486e-2*s2)/(84.850+69.024*s+s2)
    alpha1 = 49.843 - 0.2276*s + 0.198e-2*s2
    rtr15 = 1.0 + (sst-15.0)*alpha0/(alpha1+sst)
    
    sig = sig35*r15*rtr15
    
    !    permittivity
    a0 = exp(a0coef(1)*s + a0coef(2)*s2 + a0coef(3)*s*sst)  
    e0s = a0*e0
    
    if(sst.le.30) then
        b1 = 1.0 + s*(b1coef(1) + b1coef(2)*sst + b1coef(3)*sst2 + b1coef(4)*sst3 + b1coef(5)*sst4)
    else
        b1 = 1.0 + s*(9.1873715e-04 + 1.5012396e-04*(sst-30))
    endif
      
    n1s = n1*b1
    
    a1  = exp(z(7)*s + z(8)*s2 + z(9)*s*sst)
    e1s = e1*a1

    b2 = 1.0 + s*(z(10) + 0.5*z(11)*(sst + 30))
    n2s = n2*b2
    
    a2 = 1.0  + s*(z(12) + z(13)*sst)
    e2s = e2*a2
    
return
end subroutine  dielectric_meissner_wentz



end module RSS_ATM_2022_BULK_module