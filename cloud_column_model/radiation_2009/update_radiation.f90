 subroutine update_radiation (rmonth,rday,rlat,cosz,tb,qb,rho,p0,p00, &
                              pi,dz0,tland,qcl,qrn,qci,qcs,qcg,dpt,dqv,&
                              rcp,rc2p,rrp,rpp,rsp,rap,rgp,rhp,&
                              rsw,rlw,rflux,rams_micro)
 use module_ra_goddard_gce
! implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
!   This routine i) one-dimensionalizes GCE input, ii) add strat layers above GCE layers, and
!   iii) drive goddard radiation code, and iv) output heating rate and energy budgets. 
! 
! Methods: 
!  * F90 freeformat
!  * Put all GCE common block here. So there is no more common block below this point. 
!  * Developed as non-implicit codiing, but be careful it is implicite codeing now for GCE.
!  * Put horizontal loop here in order to 1-dimensionalize radiation driver. 
! 
! History:
! 03/2009  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!-----------------------------------------------------------------------------------------------------
!  D. Posselt 12/22/2009
!  Modified this to be called as a column model from an external driver
!  Stripped out unneccessary common blocks and statistics

      include 'dimensions.h'
!       parameter (nxf=260,nyf=260,nzf=43,nt=28640,itt=2) 
!      parameter (nxf=260,nyf=260,nzf=43,nt=28640,itt=182) 
!       parameter (lb=2,kb=1)
! define decomposition from rmp_switch.h
!       parameter (npes=8,ncol=16,nrow=1)
! define partial dimension for computation by decomposition
!       parameter (nx=(nxf-lb*2-1)/ncol+1+lb*2)
!       parameter (ny=(nyf-lb*2-1)/nrow+1+lb*2)
!       parameter (nz=nzf)
! define partial dimension for fft by transpose decomposition
!       parameter (nyc=(nyf-lb*2-1)/ncol+1+lb*2)
!       parameter (nyr= ny)
!       parameter (nxc= nx)
!       parameter (nxr=(nxf-lb*2-1)/nrow+1+lb*2)
!       parameter (nzc=(nzf-kb*2-1)/ncol+1+kb*2)
!       parameter (nzr=(nzf-kb*2-1)/nrow+1+kb*2)
! integer,parameter :: nxf=20,nyf=20,nzf=43,nt=28640,itt=2
! integer,parameter ::lb=2,kb=1
! define decomposition from rmp_switch.h
! integer,parameter ::npes=1,ncol=1,nrow=1
! define partial dimension for computation by decomposition
! integer,parameter ::nx=(nxf-lb*2-1)/ncol+1+lb*2
! integer,parameter ::ny=(nyf-lb*2-1)/nrow+1+lb*2
! integer,parameter ::nz=nzf
! define partial dimension for fft by transpose decomposition
! integer,parameter ::nyc=(nyf-lb*2-1)/ncol+1+lb*2
! integer,parameter ::nyr= ny
! integer,parameter ::nxc= nx
! integer,parameter ::nxr=(nxf-lb*2-1)/nrow+1+lb*2
! integer,parameter ::nzc=(nzf-kb*2-1)/ncol+1+kb*2
! integer,parameter ::nzr=(nzf-kb*2-1)/nrow+1+kb*2
 integer,parameter :: nx1=1
 integer,parameter :: ny1=1
!
! input parameter
!
 real,intent(in) :: rmonth  !month of simulation start
 real,intent(in) :: rday    !day of simulation start
 real,intent(in) :: rlat    !latitude of simulation
 real,intent(in) :: cosz !cosine of solar zenith angle (0~1) [-]
 real,intent(in) :: tland(nx,ny) !skin temperature [K]
 real,dimension(nz),intent(in) :: &
      tb , & ! mean temperature [K]
      qb , & ! mean water vapor [g/g]
      rho, & ! density [g/cm3]
      p0 , & ! layer pressure [g*cm/s2/cm2] = [dyn/cm2] = [1e-3 mb] 
      pi , & ! exner function (= (p0*1e-1/1e-5)**(287./1004.) [-]
      p00, & ! level pressure (staggered) [g*cm/s2/cm2] = [dyn/cm2] = [1e-3 mb]
      dz0    ! layer depth profile [cm]
 integer, intent(in) :: rams_micro !flag for whether to use RAMS micro fields (1=yes)
! (Both SW and LW radiation option)
! This integer parameter determines the numberr of stratosphere layer above the 
! top of actual simulation layer. If Zero, radiative transfer calculation takes place
! within the actual simulated model layers. If greater than zero, additional stratospheric
! layer can improve SW and LW radiation budget by computing molecular scattereing & absorption
! in the strasphere. 
! BUT REMEMBER, THIS ADDITION CAN SLOW DOWN THE COMPUTATION.
  integer,parameter :: alev_strat = 5    !recommended about 5 layers
  real,parameter    :: frac_p = 0.25        !fractional reductin rate of pressure

! D. Posselt -- Variable declarations
 integer, parameter :: kles=nz-1
 real,dimension(nx,ny,nz) :: qcl,qrn,qci,qcg,qcs,dpt,dqv  !Goddard Microphysics 
 real,dimension(nx,ny,nz) ::  rcp,rrp,rpp,rsp,rap,rgp,rhp,rc2p !RAMS microphysics
  real :: gce_reff(nx,ny,nz,3)   !effective radius [micron]??
  real :: rams_reff(nx,ny,nz,8)  !effective radius [micron]??
 real,dimension(nx,ny,nz) :: rsw,rlw
 real rflux(nx1,ny1,8)  !radiation budget at TOA and surface [W/m2]
 real,dimension(nx,ny) :: rldown,swdown,albedo

!
! input common block 
!
!  integer :: imax,iles1,iles,il2,jmax,jles1,jles,jl2, kmax,kles,kl2  !MPI parameter
!  common/mpi_parameter/imax,iles1,iles,il2,jmax,jles1,jles,jl2,kmax,kles,kl2 
!  real,dimension(nx,ny,nz) :: qcl,qrn,qci,qcg,qcs,dpt,dqv  !Goddard Microphysics 
!  common/b1cr/ qcl,qrn  ! cloud liquid, rain mixing ratio [g/g]
!  common/b1ig/ qci,qcg  ! cloud ice, graupel mixing ratio [g/g]
!  common/b1s/  qcs      ! snow mixing ratio [g/g]
!  common/b1tq/ dpt,dqv  ! purturbation of potential temperature and water vapor [K], [g/g]
!  real,dimension(nx,ny,nz) ::  rcp,rrp,rpp,rsp,rap,rgp,rhp,rc2p !RAMS microphysics
!  common/rams_mixr1/ rcp,rrp,rpp,rsp,rap,rgp,rhp,rc2p !mixing ratio [g/g]
! Variables from RAMS micro used in radiation
!   real :: gce_reff(nx,ny,nz,3)   !effective radius [micron]??
!   real :: rams_reff(nx,ny,nz,8)  !effective radius [micron]??
!   common/rams_reff/ gce_reff, rams_reff
!
! output common block
!
!  real,dimension(nx,ny,nz) :: rsw,rlw
!  common/slwave/   rsw, rlw  !-> output (short/longwave heating rate) for updating theta (gcempi)
!  real rflux(nx1,ny1,8)  !radiation budget at TOA and surface [W/m2]
!  common/radflux/rflux
!  real,dimension(nx,ny) :: rldown,swdown,albedo
!  common/var6/     rldown,swdown,albedo  ! rldown&swdown output to PLACE, albedo input from PLACE
!
! constant parameters 
!
 real,parameter :: qv_min     = 1.e-6   !minimum water vapor mixing ratio [g/g]
 real,parameter :: temp_min   = 190.    !minimum temperature [K]
 real,parameter :: q_thresh   = 1e-6    !minimum condensate mixing ratio to account cloud optical depth
 real,parameter :: alb_const  = 0.07    !constant albedo  (Toshi: Should be input from gce.config)
 real,parameter :: emiss_const= 0.98    !surface emissivity for IR band (Toshi: Should be input from gce.config)
!
! local parmeters
!
 integer :: i,j,k,k_strat,k_rad
 real :: julian !real of julian date
 integer,parameter :: mxlyr = nz-2+alev_strat  !actual maximum layer number for radiation routine
 real :: &
    qv     ,& !water vapor mixing ratio [g/g] 
    tskin  ,& !surface skin temperature [K]
    tsurf  ,& !near-surface air temperature [K]
    emiss  ,& !surface emissivity for longwave band []
    albe      !surface albedo []
 real,dimension(mxlyr) :: &  !adjusted one layer below from GCE to radiation routine
    t    ,& !layer temperature [K]
    p    ,& !layer pressure [Pa]
    exner,& !exner function pi []
    dz   ,& !layer depth [m]
    sh   ,& !specific humidity 
    qc1  ,& !cloud water (small mode) mixing ratio [g/g] or [kg/kg]
    qc2  ,& !cloud water (large mode)  mixing ratio [g/g] or [kg/kg]
    qi1  ,& !cloud ice (small mode) mixing ratio [g/g] or [kg/kg]
    qi2  ,& !cloud ice (large mode) mixing ratio [g/g] or [kg/kg]
    qr   ,& !rain mixing ratio [g/g] or [kg/kg]
    qs   ,& !snow mixing ratio [g/g] or [kg/kg]
    qg   ,& !graupel mixing ratio [g/g] or [kg/kg]
    qh   ,& !hail mixing ratio [g/g] or [kg/kg]
    reff_qc1,&  ! cloud water (small mode) re [micron]
    reff_qc2,&  ! cloud water (large mode) re [micron]
    reff_qi1,&  ! clouud ice (small mode) re [micron]
    reff_qi2,&  ! clouud ice (large mode) re [micron]
    reff_qr ,&  ! rain  snow re [micron]
    reff_qs ,&  ! snow re [micron] 
    reff_qg ,&  ! graupel re [micron]
    reff_qh ,&  ! hail re [micron]
    cldfra   ,& !cloud fraction (0 or 1) []
    taucldi,&   !Visible optical dpeth of cloud ice
    taucldc,&   !Visible optical depth of cloud water
    sw_thrate,&!theta tendency due to shortwave radiative heating [K/s]
    lw_thrate  !theta tendency due to longwave radiative heating [K/s]
 real,dimension(mxlyr+1) :: &  !adjusted one layer below from GCE to radiation routine
    p_lev    !level pressure [Pa]
 real,dimension(alev_strat) :: &
     p00_strat,&  !level pressure [mb]
     p0_strat, &  !layer pressure [mb]
     t_strat,  &  !layer temperature [K]
     sh_strat, &  !specific humidity []
     o3_strat     !ozone conc [g/g]
 real, dimension(1:8) :: ERBE_out  !earth radiation budget output (1 ~ 8 is index)
                                   ! 1-TOA LW down, 2-TOA LW up, 3-surface LW down, 4-surface LW up
                                    ! 5-TOA SW down, 6-TOA SW up, 7-surface SW down, 8-surface SW up
! D. Posselt put in common block with PSD parameters for 3ice bulk scheme
 real    tnw,tns,tng,roqs,roqg,roqr, cpi
 common/size/ tnw,tns,tng,roqs,roqg,roqr

 save
!
! ---------------------------------- Program Start --------------------------------------
!
!
! quick pi computation
  cpi =4.*atan(1.)

! getting julidan date
!
  call get_julian(2000.,rmonth,rday,julian)
!
! Initilization 
!
  rflux = 0. !energy budget output
!
! Mclatchy sounding climatologies are subdivided into different latitudinal region and seasons.
! So, the actual values are interpolated for a given julidan date and latitude
!
   call sounding_interp(rlat, INT(julian))
   if(alev_strat > 0) then
!
! get pressure for additional upper level fractional reduction 
!
     do k_strat = 1, alev_strat
        if(k_strat==1) then
           p00_strat(k_strat) = 0.5*p00(nz)*1e-3  !pressure level [mb] (ave)
           p0_strat(k_strat) = 0.5*(p00_strat(k_strat) + p00(nz)*1e-3) !pressure layer [mb] (ave)
        else
           p00_strat(k_strat) = frac_p*p00_strat(k_strat-1) !pressure level[mb]
           p0_strat(k_strat)  = 0.5*( p00_strat(k_strat-1) + p00_strat(k_strat) ) !pressure layer[mb] (ave)
        endif
     enddo 
!
!get temp, specific humidity, and ozone profiles profiles from sounding (actually ozone is not used here)
!
     call sounding_strat( alev_strat, p0_strat,t_strat,sh_strat, o3_strat) 
   endif
!
! Start horizontal loop here (1 dimensionalize input parameters)
!
 j_loop: do j=1,1
    i_loop: do i=1,1
!
! surface parameter
!
       albe=alb_const    ! surface albedo from constant parameter
       emiss = emiss_const  ! surface emissivity for IR band
       tskin = tland(i,j)   ! skin temp [K]
! -------------------------Vertical Profile of GCE and Radiation --------------------------
!                         GCE                 |                  Radiation
!  (gce k)     (Layer)            (Level)     |  (radiation k)   (layer)     (level)
!  k = nz                        top level    |
!  k = kles    top layer           .          |   k_rad = nz-1               surface             
!  k = .        .                  .          |   k_rad = nz-2   bottom  
!  k = 3        .                  .          |
!  k = 2       bottom layer      bottom level |
!  k = 1                                      |   k_rad = 1      top layer   top level
!                                             |
!--------------------------------------------------------------------------------------------
!
! GCE profile
!
       do k=2,kles !gce layer profile  -> radiation level profile (reverse)
          k_rad = kles-k+1 + alev_strat
          qv=max(qv_min,qb(k)+dqv(i,j,k)) !water vapor [g/g]
          sh(k_rad)= qv/(1.+qv)  ! specific humidity []
          t(k_rad) =  max(temp_min, (tb(k)+dpt(i,j,k))*pi(k) )     ! layer temperature [K]
          p(k_rad)   = p0(k) *1e-3     ! layer pressure [mb or hPa] <- [g*cm/s2/cm2]
          dz(k_rad)  = dz0(k)*1e-2     ! layer depth [m] <- [cm]
          exner(k_rad) = pi(k)  !(p0(k)*1e-1/1e+5)**(287./1004.)       ! exner function [-]
       enddo
       do k=2,kles+1 !gce level profile -> radiation level profile (reverse)
         k_rad = (kles+1)-k+1+alev_strat
         p_lev(k_rad) = p00(k)*1e-3     ! level pressure [mb or hPa]  <-  [g*cm/s2/cm2]
       enddo

       if ( RAMS_MICRO .eq. 0 ) then

       do k=2,kles !gce layer profile  -> radiation level profile (reverse)
          k_rad = kles-k+1 + alev_strat
          qc1(k_rad) = qcl(i,j,k)    ! cloud water (small mode) [g/g]
          qc2(k_rad) = 0.            ! cloud water (large mode) [g/g]
          qi1(k_rad) = qci(i,j,k)    ! cloud ice (small mode) [g/g]
          qi2(k_rad) = 0.            ! cloud ice (large mode) [g/g]
          qr (k_rad) = qrn(i,j,k)    ! rain [g/g]
          qs (k_rad) = qcs(i,j,k)    ! snow aggregate [g/g]
          qg (k_rad) = qcg(i,j,k)    ! graupel [g/g]
          qh (k_rad) = 0.            ! hail [g/g]
          reff_qc1(k_rad)=10.                    ! cloud water (small mode) re [micron]
          reff_qc2(k_rad)=0.                     ! cloud water (large mode) re [micron]
          reff_qi1(k_rad)=20.                    ! clouud ice (small mode) re [micron]
          reff_qi2(k_rad)=0.                     ! clouud ice (large mode) re [micron]
! Set these as a first guess--these are defaults if mixing ratios are less than qv_min
          reff_qr (k_rad)=1000.                  ! rain re [micron]
          reff_qs (k_rad)=1000.                  ! snow aggrefate re [micron] 
          reff_qg (k_rad)=1000.                  ! graupel re [micron]
! D.Posselt Try computing effective radius from assumed MP PSD... (as in opt4.F)
          if (qr(k_rad) .gt. 0.1*q_thresh) reff_qr (k_rad)=3./((cpi*tnw*roqr/(rho(k)*qr(k_rad)))**.25) * 1.e4 ! rain re [micron]
          if (qs(k_rad) .gt. 0.1*q_thresh) reff_qs (k_rad)=3./((cpi*tns*roqs/(rho(k)*qs(k_rad)))**.25) * 1.e4 ! snow aggregate re [micron] 
          if (qg(k_rad) .gt. 0.1*q_thresh) reff_qg (k_rad)=3./((cpi*tng*roqg/(rho(k)*qg(k_rad)))**.25) * 1.e4 ! graupel re [micron]
          reff_qh (k_rad)=0.                     ! hail re [micron]
!           print '(a,i5,3f20.10)','Re for rain, snow, graupel ',k_rad,reff_qr(k_rad),reff_qs(k_rad),reff_qg(k_rad)
          if( (qc1(k_rad)+qc2(k_rad)+qr(k_rad)+qi1(k_rad)+qi2(k_rad)+qs(k_rad)+qg(k_rad)+qh(k_rad))  > q_thresh ) then
               cldfra(k_rad) = 1.    !cloud fraction 1 
          else 
               cldfra(k_rad) = 0.    !no cloud
          endif
       enddo

       else if ( RAMS_MICRO .eq. 1 ) then

!(Toshi- This needs modification lator)
       do k=2,kles !gce layer profile  -> radiation level profile (reverse)
          k_rad = kles-k+1 + alev_strat
          qc1(k_rad) = rcp(i,j,k)     ! cloud water (small mode) [g/g]
          qc2(k_rad) = rc2p(i,j,k)    ! cloud water (large mode) [g/g]
          qi1(k_rad) = rpp(i,j,k)     ! cloud ice (small mode) [g/g]
          qi2(k_rad) = rsp(i,j,k)     ! cloud ice (large mode) [g/g]
          qr (k_rad) = rrp(i,j,k)     ! rain [g/g]
          qs (k_rad) = rap(i,j,k)     ! snow aggregate [g/g]
          qg (k_rad) = rgp(i,j,k)     ! graupel [g/g]
          qh (k_rad) = rhp(i,j,k)     ! hail [g/g]
          reff_qc1(k_rad)=rams_reff(i,j,k,1)     ! cloud water (small mode) re [micron]
          reff_qc2(k_rad)=rams_reff(i,j,k,8)     ! cloud water (large mode) re [micron]
          reff_qi1(k_rad)=rams_reff(i,j,k,3)     ! clouud ice (small mode) re [micron]
          reff_qi2(k_rad)=rams_reff(i,j,k,4)     ! clouud ice (large mode) re [micron]
          reff_qr (k_rad)=rams_reff(i,j,k,2)     ! rain re [micron]
          reff_qs (k_rad)=rams_reff(i,j,k,5)     ! snow aggrefate re [micron] 
          reff_qg (k_rad)=rams_reff(i,j,k,6)     ! graupel re [micron]
          reff_qh (k_rad)=rams_reff(i,j,k,7)     ! hail re [micron]
          if( (qc1(k_rad)+qc2(k_rad)+qr(k_rad)+qi1(k_rad)+qi2(k_rad)+qs(k_rad)+qg(k_rad)+qh(k_rad))  > q_thresh ) then
               cldfra(k_rad) = 1.    !cloud fraction 1 
          else
               cldfra(k_rad) = 0.    !no cloud
          endif
       enddo

       endif
!
! Stratospheric profile
!
       do k_strat = 1, alev_strat
          k_rad = alev_strat - k_strat + 1
          qc1(k_rad) = 0.    ! cloud water (small mode) [g/g]
          qc2(k_rad) = 0.    ! cloud water (large mode) [g/g]
          qi1(k_rad) = 0.    ! cloud ice (small mode) [g/g]
          qi2(k_rad) = 0.    ! cloud ice (large mode) [g/g]
          qr (k_rad) = 0.    ! rain [g/g]
          qs (k_rad) = 0.    ! snow aggregate [g/g]
          qg (k_rad) = 0.    ! graupel [g/g]
          qh (k_rad) = 0.    ! hail [g/g]
          reff_qc1(k_rad)=0.     ! cloud water (small mode) re [micron]
          reff_qc2(k_rad)=0.     ! cloud water (large mode) re [micron]
          reff_qi1(k_rad)=0.     ! clouud ice (small mode) re [micron]
          reff_qi2(k_rad)=0.     ! clouud ice (large mode) re [micron]
          reff_qr (k_rad)=0.     ! rain re [micron]
          reff_qs (k_rad)=0.     ! snow aggrefate re [micron] 
          reff_qg (k_rad)=0.     ! graupel re [micron]
          reff_qh (k_rad)=0.     ! hail re [micron]
          cldfra(k_rad) = 0. !cloud fraction zero.
          sh(k_rad)=sh_strat(k_strat)  !specific humidity []
          t (k_rad)=t_strat (k_strat)  ! layer temperature [K]
          p (k_rad)=p0_strat(k_strat)  ! layer pressure [mb or hPa]
          exner(k_rad) = (p0_strat(k_strat)*1e+2/1e+5)**(287.04/1004.)       ! exner function [-]
          p_lev(k_rad) = p00_strat(k_strat)  !level pressure [mb or hPa]
       enddo
       do k_strat = 1, alev_strat
          k_rad =  alev_strat - k_strat + 1
         dz(k_rad)=t(k_rad)*287.04/(-9.8)*log(p(k_rad)/p(k_rad+1)) !may modify lator
       enddo
       tsurf = 0.5*(tland(i,j) + t(mxlyr))  !near-surface air temperature needed for IR raddriv (carefull)
!
! One dimensional driver of goddard radiation (no common block below in this subroutine)
! Use key-ward for subroutine 
!
        call goddardrad(                                         &
               mxlyr=mxlyr,tskin=tskin,tsurf=tsurf               &
              ,t=t,p=p,p_lev=p_lev,pi=exner,dz=dz,sh=sh          & 
              ,emiss=emiss,alb=albe,cosz=cosz,fcld=cldfra        &
              ,xlat=rlat,solcon=1378.                            &
              ,qc1=qc1,qc2=qc2,qi1=qi1,qi2=qi2                 & 
              ,qr=qr,qs=qs,qg=qg,qh=qh                           &
              ,reff_qc1=reff_qc1 &  
              ,reff_qc2=reff_qc2 &  
              ,reff_qi1=reff_qi1 &  
              ,reff_qi2=reff_qi2 &  
              ,reff_qr=reff_qr   &  
              ,reff_qs=reff_qs   &  
              ,reff_qg=reff_qg   &  
              ,reff_qh=reff_qh   &  
              ,ERBE_out=ERBE_out                                & !output
              ,taucldi=taucldi,taucldc=taucldc                  & !output
              ,sw_thrate=sw_thrate,lw_thrate=lw_thrate          & !output       
                                                   )
!
! output heating rate for gcempi
!
        do k=2,kles
          k_rad = kles-k+1 + alev_strat
           rsw(i,j,k)=sw_thrate(k_rad)  ! theta tendency due to shortwave radiative heating [K/s] for gcempi
           rlw(i,j,k)=-lw_thrate(k_rad)  ! theta tendency due to longwave radiative heating [K/s] for gcempi
        enddo
!
! output energy buget [W/m2] (analysis purpose)
!
       rflux(i,j,1)=ERBE_out(7) ! downward SW surface (+:down)
       rflux(i,j,2)=ERBE_out(3) ! downward LW surface
       rflux(i,j,3)=ERBE_out(5) ! downward SW TOA     (+:down)
       rflux(i,j,4)=ERBE_out(6) ! upward SW TOA (-:up) : (net SW flux at TOA) - downward SW flux TOA
       rflux(i,j,5)=ERBE_out(2) ! upward LW TOA
       rflux(i,j,6)=ERBE_out(8) ! upward SW surface (-:up) : (net SW flux at surface) - downward SW flux surface
       rflux(i,j,7)=ERBE_out(4) ! upward LW surface
       rflux(i,j,8)=ERBE_out(1) ! downward LW TOA
!
!
! for PLACE
!
       swdown(i,j) = ERBE_out(7)   !downwelling shortwave radiation at surface [W/m2] for PLACE
       rldown(i,j) = ERBE_out(3)   !downwelling longwave radiation at surface [W/m2] for PLACE
    enddo i_loop
 enddo j_loop 
!  print*,rflux(10,10,:)
  return
 end subroutine update_radiation
 subroutine get_julian(yyyy,mm,dd,julian)
 implicit none
!--------------------------------------------------------------------------------------------
! Comments:  
! Compute julian date from Month, Day, and Year.  
! 
! History:
! 12/2007  Toshi Matsui@NASA GSFC ; Initial   
!           
! References: 
!-----------------------------------------------------------------------------------------------------
 real,intent(in) :: yyyy,mm,dd !year,month,day
 real,intent(out) :: julian    !julian day
 real,parameter :: day(12) = (/31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./) !for normal year
 real,parameter :: dayl(12)= (/31.,29.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./) !for leap year  
  if( mod(yyyy,4.) == 0. ) then !leap year
     julian=sum( dayl(1:(int(mm)-1) ) ) + dd
  else                           ! normal year
     julian=sum( day(1:(int(mm)-1) ) ) + dd
  endif
 end subroutine get_julian
