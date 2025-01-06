! This F90 subroutine serves as the driver for the CRM microphysics
! and radiation routines. It assumes that all of the forcing fields
! have been read in and interpolated to the CRM timestep.
!
! Inputs:
!
! zlevs       vertical levels (in cm)
! dpt_int     perturbation potential temperature
! ta_int      domain mean potential temperature
! dqv_int     perturbation vapor mixing ratio
! qa_int      domain mean vapor mixing ratio
! rho_int     density
! rrho_int    1./rho
! pi_int      exner function
! p0_int      pressure
! p00_int     exp(p0)
! tland_int   Skin temperature
!
! Outputs:
!
!
! External routines called
!
! consatrh    sets up parameters for microphysics
! saticerh    3ICE microphysics routine
! zang        computes cosine of the zenith angle for radiation
! pradrat     radiation driver
! update_radiation  new (Goddard) radiation driver
!
! Procedure:
!
! 1.  Allocate column arrays for driver. Note: saticerh and all radiation
!     routines must be modified so that loops are over i=1,1, j=1,1.
! 2.  Loop over all timesteps
! 3.  In each loop iteration:
!     --  compute the solar zenith angle
!     --  call microphysics
!     --  call radiation
!     --  update the time (hour, day, month) for radiation
! 4.  Depending on the desired output frequency, save model output to file
!     (or into an array for saving later?)
!
! Derek Posselt
! University of Michigan
! 25 July 2008
!
!-----------------------------------------------------------------------
! 22 December 2009
! D. Posselt
! 
! Modified this to call Toshi's updated Goddard radiation (radiation_2009)
! and started to put placeholders in to interface RAMS microphysics. Yet to do
! on the RAMS scheme is to modify the observation simulator routine and output files, 
! and put in the advection
! 
!-----------------------------------------------------------------------
! 12 November 2010
! D. Posselt
! 
! Modified this so that the environment (temperature, water vapor, pressure) varies
! with time. The convective case still has an invariant base state environment,
! but the frontal case varies with time.
! 
! I have also modified the code to read in a variable (still constant) forcing
! time interval "tend_dt" (previously this was hard-coded to an hour).
! 
!-----------------------------------------------------------------------
! 02 December 2010
! D. Posselt
! 
! Modified the driver routine to accept a second set of moisture forcings
! This is with the intent of adding drying at lower tropospheric levels 
! in the frontal case, but could in theory be used to add a second 
! positive source to a case. Note that there is no associated second 
! vertical motion profile at this time.
! 
!-----------------------------------------------------------------------
! 18 March 2013
! D. Posselt
!
! Added the option to use time-averaged observations (average_obs = .true.)
! with start time and end time set by obs_start_time and obs_end_time, 
! in units of seconds.
!
!-----------------------------------------------------------------------
subroutine drive_micro_rad(nz_in,dz_in,nt,nt_out,dt,rad_dt,out_dt,ntend, tend_dt, tsfc, &
                           t_lev, t_lay, p_lev, p_lay, q_lev, q_lay, &
                           rho_lev, rho_lay, w, z_lev_in, z_lay_in, &
                           w_max_f, q_max_f, z_w_f, z_q_f, half_w_f, half_q_f, &
                           q_max_f2, z_q_f2, half_q_f2, neg_2nd_f, &
                           start_month, start_day, start_hour, latitude, &
                           nparams, params_in, max_obs, n_obs_times, obs_times, &
                           nobs, obs_array, write_output, print_more, rams_micro, new_rad, &
                           ri_full, ar_full, lwp_full, iwp_full, lwrad_full, swrad_full, &
                           case_type, improve, &
                           average_obs, obs_start_time, obs_end_time )
implicit none

include 'include/dimensions.h'

! Input variables
integer, intent(in) :: nz_in, nt, nt_out, ntend, nparams, max_obs, n_obs_times, nobs, &
                       rams_micro, new_rad, case_type, improve
real,    intent(in) :: dt, rad_dt, out_dt, tend_dt, tsfc
real,    intent(in) :: start_month, start_day, start_hour, latitude
real,    intent(in) :: dz_in
real, dimension(nz_in,nt), intent(in) :: z_lev_in,   z_lay_in   & ! Heights       at layer boundaries and midpoints
                                       , T_lev,   T_lay   & ! Temperature   at layer boundaries and midpoints
                                       , q_lev,   q_lay   & ! Vapor mix rat at layer boundaries and midpoints
                                       , p_lev,   p_lay   & ! Pressure      at layer boundaries and midpoints
                                       , rho_lev, rho_lay   ! Density       at layer boundaries and midpoints

real, dimension(nparams), intent(in) :: params_in   ! Microphysics parameters
real, dimension(max_obs), intent(in) :: obs_times   ! Observation times (minutes)

logical, intent(in)                  :: average_obs    ! whether to average observations over "nobs_times" time intervals or not
real,    dimension(max_obs), intent(in) :: obs_start_time & ! Start times for accumulating observations 
                                         , obs_end_time     ! End times for accumulating observations 

! real, dimension(nz_in), intent(in) :: w                & ! Vertical velocity (layer mid-points)
!                                     , qv_sorc            ! Vapor mixing ratio source term (layer mid-points)
real, dimension(nz_in) :: w                & ! Vertical velocity (layer mid-points)
                        , qv_sorc          & ! Vapor mixing ratio source term (layer mid-points)
                        , qv_sorc2           ! Vapor mixing ratio source term (layer mid-points)

real, dimension(ntend), intent(in) :: w_max_f   & ! Max vertical velocity (m/s)
                                    , q_max_f   & ! Max qv forcing (g/kg/hour)
                                    , z_w_f     & ! Height of max vertical velocity (km)
                                    , z_q_f     & ! Height of max qv forcing (km)
                                    , half_w_f  & ! Half-width of vertical velocity profile (km)
                                    , half_q_f  & ! Half-width of qv forcing profile (km)
! Second set of moisture forcings...
                                    , q_max_f2  & ! Max of 2nd qv forcing (g/kg/hour)
                                    , z_q_f2    & ! Height of 2nd max qv forcing (km)
                                    , half_q_f2   ! Half-width of 2nd qv forcing profile (km)

logical, intent(in) :: neg_2nd_f      ! Whether 2nd set of water vapor forcing fields should be negative

logical, intent(in) :: print_more, write_output

! Output from the model--observation array
real    :: obs_array(nobs, n_obs_times)

! Input/Output variables
real, dimension(1,1,nz_in) :: qcl,  qrn,  qci,  qcs,  qcg,  &
                           qcl1, qrn1, qci1, qcs1, qcg1, &
                           rsw, rlw, dpt, dqv, dpt_save, dqv_save, &
                           rcp, rc2p, rrp, rpp, rsp, rap, rgp, rhp, &
                           rcp1,rc2p1,rrp1,rpp1,rsp1,rap1,rgp1,rhp1, &
                           tair_out, tairc_out, theta_out, &
                           delta2, delta3, delta4
real, dimension(1,1,8) :: rflux
! real, dimension(1,1)    :: ri, ar
real    :: ri, ar

! Output variables
real    :: pwv, lwp, iwp, ctop
real, dimension(nz_in,nt_out) :: qcl_full, qrn_full, qci_full &
                               , qcs_full, qcg_full, rsw_full, rlw_full &
                               , dpt_full_micro, dpt_full, dqv_full &
                               , dpt_full_rad, w_full, qv_sorc_full, qv_sorc2_full &
                               , rcp_full, rc2p_full, rrp_full, rpp_full &
                               , rsp_full, rap_full, rgp_full, rhp_full &
                               , ta_full, qa_full, p0_full, p00_full, pi_full &
                               , rho_full, tair_full, tairc_full, theta_full &
                               , delta2_full, delta3_full, delta4_full
real, dimension(nt_out)       :: ri_full, ar_full &
                               , pwv_full, lwp_full, iwp_full &
                               , lwrad_full, swrad_full, ctop_full

! Local variables
real   :: dz
real, dimension(nz_in) :: z_lev, z_lay  ! Heights at layer boundaries and midpoints

! Forcing increments and current values
real :: w_max_inc, q_max_inc, z_w_inc, z_q_inc, half_w_inc, half_q_inc
real :: w_max, q_max, z_w, z_q, half_w, half_q
real :: q_max_inc2, z_q_inc2, half_q_inc2
real :: q_max2, z_q2, half_q2
! Microphysics variables
real, dimension(nz_in) :: f0, dz0, y1, y2, pi, p0, p00, ta, qa, &
                       rho, rrho, rho1, rrho1, am, am1, fv, fv1
real, dimension(1,1,nz_in) :: df0, vtp, w1, ak, ak1, ak2, ak_f, dfc, xxw, scr3d
! Radiation variables
real, dimension(1,1)    :: tland
real    :: hrl, riday, rmonth, cosz, rlat
real    :: al, cp, half_pi
! For computation of ak (diffusion) coefficient
real    :: y1k, qak, tak, a1, a2, xy, xym, xyp, ta1k, qa1k
! Indices, flags, etc
integer :: i, j, k, km, kp, t, tt, f1, f2, ob, iopcloud, iradave, iflag, rad_freq
integer :: irsg, ismg
integer :: ibeg, iend
real    :: nt_obs

! Threshold for detecting cloud in simobs
real, parameter :: twcz = 1.e-6 ! in g/g

! Check vertical dimension immediately
if ( nz_in /= nz ) then
  print '(a)','Computed nz not equal to input nz!'
  print '(a,2i10)','nz input, nz computed: ',nz, nz_in
  stop 'Check file include/dimensions.h'
endif

! Set integer values
iopcloud = 0
iradave  = 0

! Set values of constants
  cp = 1.004e7
  al = 2.5e10
  half_pi = 1.0 * acos(0.0)

! Compute radiation call frequency
rad_freq = int(rad_dt / dt)
if ( rad_freq .lt. 1 ) rad_freq = 1

! Set initial w_max, q_max, z_w, z_q, half_w, half_q
w_max   = w_max_f(1)
q_max   = q_max_f(1)
z_w     = z_w_f(1) * 100.e3 ! convert from km to cm
z_q     = z_q_f(1) * 100.e3 ! convert from km to cm
half_w  = half_w_f(1) * 100.e3 ! convert from km to cm
half_q  = half_q_f(1) * 100.e3 ! convert from km to cm
q_max2  = q_max_f2(1)
z_q2    = z_q_f2(1) * 100.e3 ! convert from km to cm
half_q2 = half_q_f2(1) * 100.e3 ! convert from km to cm

! Compute GCE profile variables for the first timestep
do k = 1, nz
  z_lev(k)  = z_lev_in(k,1) * 100.    ! Convert heights to cm (GCE works in cgs)
  z_lay(k)  = z_lay_in(k,1) * 100.    ! Convert heights to cm (GCE works in cgs)
  rho(k)    = rho_lay(k,1) / 1000. ! rho is density in cgs at layer midpoints
  rrho(k)   = 1.0 / rho(k)       ! rrho is 1.0/density in cgs at layer midpoints
  rho1(k)   = rho_lev(k,1) / 1000. ! rho1 is density in cgs at layer boundaries
  rrho1(k)  = 1.0 / rho1(k)      ! rrho1 is 1.0/density in cgs at layer boundaries
  ta(k)     = t_lay(k,1) * ( 1.e5 / p_lay(k,1) ) ** ( 0.286 ) ! Base state potential temperature
  qa(k)     = q_lay(k,1)           ! Base state vapor mixing ratio
  p0(k)     = p_lay(k,1) * 10.0    ! Base state pressure (cgs) at layer midpoints
  p00(k)    = p_lev(k,1) * 10.0    ! Base state pressure (cgs) at layer boundaries
  pi(k)     = ( p0(k) / 1000.e3 )**.286 ! Exner function at layer midpoints
  fv(k)     = sqrt(rho(2) * rrho(k))   ! Fall velocity constant at layer midpoints
  fv1(k)    = sqrt(rho1(2) * rrho1(k)) ! Fall velocity constant at layer boundaries
enddo

! tland(1,1) = tsfc ! Surface temperature
tland(1,1) = t_lev(1,1) ! Surface temperature

! Set indices on forcing
f1 = 1
f2 = 2

! Set initial increments on the above 4 variables
w_max_inc   = ( w_max_f(f2)  - w_max_f(f1) )  / (tend_dt / dt + 1.)
q_max_inc   = ( q_max_f(f2)  - q_max_f(f1) )  / (tend_dt / dt + 1.)
z_w_inc     = ( z_w_f(f2)    - z_w_f(f1) )    / (tend_dt / dt + 1.) * 100.e3 ! convert from km to cm
z_q_inc     = ( z_q_f(f2)    - z_q_f(f1) )    / (tend_dt / dt + 1.) * 100.e3 ! convert from km to cm
half_w_inc  = ( half_w_f(f2) - half_w_f(f1) ) / (tend_dt / dt + 1.) * 100.e3 ! convert from km to cm
half_q_inc  = ( half_q_f(f2) - half_q_f(f1) ) / (tend_dt / dt + 1.) * 100.e3 ! convert from km to cm

! Second set of water vapor forcing
q_max_inc2  = ( q_max_f2(f2) - q_max_f2(f1) )   / (tend_dt / dt + 1.)
z_q_inc2    = ( z_q_f2(f2) - z_q_f2(f1) )       / (tend_dt / dt + 1.) * 100.e3 ! convert from km to cm
half_q_inc2 = ( half_q_f2(f2) - half_q_f2(f1) ) / (tend_dt / dt + 1.) * 100.e3 ! convert from km to cm

if ( print_more ) print '(a,9f20.10)','increments: ', &
          w_max_inc, q_max_inc, z_w_inc, z_q_inc, half_w_inc, half_q_inc, q_max_inc2, z_q_inc2, half_q_inc2

! Compute initial vertical velocity and vapor mixing ratio tendency
w(1)        = 0.0
qv_sorc(1)  = 0.0
qv_sorc2(1) = 0.0
do k = 2, nz-1
  if ( abs(z_lay(k) - z_w) <= half_w ) then
    if ( case_type .eq. 1 ) then ! Convective case uses cos^4
      w(k) = w_max * ( cos ( half_pi * ( z_lay(k) - z_w ) / half_w ) )**4
    else if ( case_type .eq. 2 ) then ! Frontal case uses cos^6
      w(k) = w_max * ( cos ( half_pi * ( z_lay(k) - z_w ) / half_w ) )**6
    endif
  else
    w(k)      = 0.0
  endif
  if ( abs(z_lay(k) - z_q) <= half_q ) then
    if ( case_type .eq. 1 ) then ! Convective case uses cos^2
      qv_sorc(k)  = q_max  * ( cos ( half_pi * ( z_lay(k) - z_q  ) / half_q  ) )**2   &
                / ( 1.e3 * tend_dt ) ! Convert from g/kg/tend_dt to kg/kg/s
    else if ( case_type .eq. 2 ) then ! Frontal case uses cos^6
      qv_sorc(k)  = q_max  * ( cos ( half_pi * ( z_lay(k) - z_q  ) / half_q  ) )**6   &
                / ( 1.e3 * tend_dt ) ! Convert from g/kg/tend_dt to kg/kg/s
    endif
  else
    qv_sorc(k)  = 0.0
  endif
  if ( abs(z_lay(k) - z_q2) <= half_q2 ) then
    if ( case_type .eq. 1 ) then ! Convective case uses cos^2
      qv_sorc2(k) = q_max2 * ( cos ( half_pi * ( z_lay(k) - z_q2 ) / half_q2 ) )**2   &
                / ( 1.e3 * tend_dt ) ! Convert from g/kg/tend_dt to kg/kg/s
      if ( neg_2nd_f ) qv_sorc2(k) = qv_sorc2(k) * (-1.0)
    else if ( case_type .eq. 2 ) then ! Frontal case uses cos^6
      qv_sorc2(k) = q_max2 * ( cos ( half_pi * ( z_lay(k) - z_q2 ) / half_q2 ) )**6   &
                / ( 1.e3 * tend_dt ) ! Convert from g/kg/tend_dt to kg/kg/s
      if ( neg_2nd_f ) qv_sorc2(k) = qv_sorc2(k) * (-1.0)
    endif
  else
    qv_sorc2(k) = 0.0
  endif
enddo

! Set GCE vertical velocity
do k = 1, nz
  w1(1,1,k) = w(k) * 100.        ! Convert from m/s to cm/s
enddo

! Set dz0 (distance between full levels)
dz     = dz_in * 100. ! Convert heights to cm (GCE works in cgs)
dz0(1) = dz
do k = 2, nz-1
  dz0(k) = z_lev(k+1) - z_lev(k)
  am(k)=dz/(z_lev(k+1) - z_lev(k))
  am1(k)=dz/(z_lay(k)-z_lay(k-1))
enddo
am(1) = am(2)
am1(2) = am1(3)
am1(1) = am1(2)
am(nz) = am(nz-1)
am1(nz) = am1(nz-1)
dz0(nz) = dz0(nz-1)

if (print_more) then
  print '(a)','   k        rho      rho1      rrho     rrho1        ta        qa        p0       p00      pi        fv       fv1        am       am1       dz0'
  do k = 1, nz
    print '(i5,2F10.6,2F10.2,F10.2,F10.6,2F10.2,F10.3,2F10.5,3F10.3)', &
            k, rho(k), rho1(k), rrho(k), rrho1(k), ta(k), qa(k), &
            p0(k), p00(k), pi(k), fv(k), fv1(k), am(k), am1(k), dz0(k)
  enddo
endif

! stop

do k=2,nz-1
  f0(k) = al / (pi(k) * cp)
enddo
f0(1) = f0(2)
f0(nz) = 2. * f0(nz-1) - f0(nz-2)

y1(1) = 0.0
y2(1) = 0.0
do k=2,nz
  y1(k)=f0(k)+f0(k-1)
  y2(k)=.5*y1(k)
enddo

! Set initial values of variables
dpt = 0.0
dqv = 0.0
qcl = 0.0
qrn = 0.0
qci = 0.0
qcs = 0.0
qcg = 0.0
qcl1 = 0.0
qrn1 = 0.0
qci1 = 0.0
qcs1 = 0.0
qcg1 = 0.0
rcp = 0.0
rcp1 = 0.0
rc2p = 0.0
rc2p1 = 0.0
rrp = 0.0
rrp1 = 0.0
rpp = 0.0
rpp1 = 0.0
rsp = 0.0
rsp1 = 0.0
rap = 0.0
rap1 = 0.0
rgp = 0.0
rgp1 = 0.0
rhp = 0.0
rhp1 = 0.0
vtp  = 0.0
ar   = 0.0
ri   = 0.0

! Call consat routine
if (print_more) print '(a)','Calling consatrh'
call consatrh(improve,params_in)

! If we will be calling radiation at all...
if ( rad_freq .le. nt .or. n_obs_times .ge. 1 ) then 

  ! Set initial times for solar zenith angle computation
  hrl    = start_hour
  riday  = start_day
  rmonth = start_month
  rlat   = latitude
  
  ! Call zang
  ! print '(a)','Calling zang'
  ! print *,'rmonth,riday,hrl,rlat: ',rmonth,riday,hrl,rlat
  call zang (rmonth,riday,cosz,rlat,hrl)
  if ( print_more ) print*,'cosz: ',cosz

  ! Call pradrat with flag of zero if old radiation is called for
  if ( new_rad .ne. 1 ) then
    if ( print_more ) print '(a)','Calling pradrat'

    iflag = 0
    call pradrat (iflag,iopcloud,iradave,cosz,ta,qa,rho,p0,p00,pi,dz0, &
                  tland,qcl,qrn,qci,qcs,qcg,dpt,dqv,rsw,rlw,rflux,improve)
    if ( print_more ) print*,'After pradrat, line 342, rflux = ',rflux
  endif

endif

! Fill first position of output arrays
tt = 1
do k = 1, nz
  ta_full(k,tt)        = ta(k)
  qa_full(k,tt)        = qa(k)
  pi_full(k,tt)        = pi(k)
  p0_full(k,tt)        = p0(k)
  p00_full(k,tt)       = p00(k)
  rho_full(k,tt)       = rho(k)
  dpt_full_micro(k,tt) = dpt(1,1,k)
  dpt_full_rad(k,tt)   = dpt(1,1,k)
  dpt_full(k,tt)       = dpt(1,1,k)
  dqv_full(k,tt)       = dqv(1,1,k)
  qv_sorc_full(k,tt)   = qv_sorc(k)
  qv_sorc2_full(k,tt)  = qv_sorc2(k)
  w_full(k,tt)         = w1(1,1,k)
  qcl_full(k,tt)       = qcl(1,1,k)
  qrn_full(k,tt)       = qrn(1,1,k)
  qci_full(k,tt)       = qci(1,1,k)
  qcs_full(k,tt)       = qcs(1,1,k)
  qcg_full(k,tt)       = qcg(1,1,k)
  rlw_full(k,tt)       = rlw(1,1,k)
  rsw_full(k,tt)       = rsw(1,1,k)
  rcp_full(k,tt)       = rcp(1,1,k)
  rc2p_full(k,tt)      = rc2p(1,1,k)
  rrp_full(k,tt)       = rrp(1,1,k)
  rpp_full(k,tt)       = rpp(1,1,k)
  rsp_full(k,tt)       = rsp(1,1,k)
  rap_full(k,tt)       = rap(1,1,k)
  rgp_full(k,tt)       = rgp(1,1,k)
  rhp_full(k,tt)       = rhp(1,1,k)
enddo
ri_full(tt) = ri
ar_full(tt) = ar

! Call observation simulation routine for the first time
call simobs(nz, twcz, pi, dpt, ta, dqv, qa, &
            qcl, qrn, qci, qcs, qcg, pwv, lwp, iwp, ctop)
pwv_full(tt)   = pwv
lwp_full(tt)   = lwp
iwp_full(tt)   = iwp
! TOA LW and SW
! Set long and shortwave radiation positive
lwrad_full(tt) = abs(rflux(1,1,5))
swrad_full(tt) = abs(rflux(1,1,4))

! Set index in observation array
ob = 1

! Start loop over timesteps
do t = 1, nt

  if ( real(t)/1000. - t/1000 .eq. 0. .and. print_more ) print '(a,i10)','Timestep: ',t

  ! If we are past the first timestep, compute GCE profile variables
  if ( t .gt. 1 ) then
    do k = 1, nz
      z_lev(k)  = z_lev_in(k,t) * 100.    ! Convert heights to cm (GCE works in cgs)
      z_lay(k)  = z_lay_in(k,t) * 100.    ! Convert heights to cm (GCE works in cgs)
      rho(k)    = rho_lay(k,t) / 1000. ! rho is density in cgs at layer midpoints
      rrho(k)   = 1.0 / rho(k)       ! rrho is 1.0/density in cgs at layer midpoints
      rho1(k)   = rho_lev(k,t) / 1000. ! rho1 is density in cgs at layer boundaries
      rrho1(k)  = 1.0 / rho1(k)      ! rrho1 is 1.0/density in cgs at layer boundaries
      ta(k)     = t_lay(k,t) * ( 1.e5 / p_lay(k,t) ) ** ( 0.286 ) ! Base state potential temperature
      qa(k)     = q_lay(k,t)           ! Base state vapor mixing ratio
      p0(k)     = p_lay(k,t) * 10.0    ! Base state pressure (cgs) at layer midpoints
      p00(k)    = p_lev(k,t) * 10.0    ! Base state pressure (cgs) at layer boundaries
      pi(k)     = ( p0(k) / 1000.e3 )**.286 ! Exner function at layer midpoints
      fv(k)     = sqrt(rho(2) * rrho(k))   ! Fall velocity constant at layer midpoints
      fv1(k)    = sqrt(rho1(2) * rrho1(k)) ! Fall velocity constant at layer boundaries
      tland(1,1) = t_lev(1,t) ! Surface temperature
    enddo

  endif

  ! Add vapor tendencies in -- include second set of forcing...
  do k = 1, nz
    dqv(1,1,k) = dqv(1,1,k) + ( ( qv_sorc(k) + qv_sorc2(k) ) * dt )
    ! Ensure non-negative qv...
    if ( qa(k)+dqv(1,1,k) .lt. 0. ) dqv(1,1,k) = dqv(1,1,k) - (qa(k)+dqv(1,1,k)) ! The value of dqv that makes the total vapor = zero...
  enddo

  ! Set up diffusion coefficient
!  print '(a)','Setting up diffusion coefficient'
  xym=(dpt(1,1,1)+ta(1))*(1.+.61*(dqv(1,1,1)+qa(1))          &
        -qcl1(1,1,1)-qrn1(1,1,1)-qci1(1,1,1)-qcs1(1,1,1)-qcg1(1,1,1))

  xy=(dpt(1,1,2)+ta(2))*(1.+.61*(dqv(1,1,2)+qa(2))           &
        -qcl1(1,1,2)-qrn1(1,1,2)-qci1(1,1,2)-qcs1(1,1,2)-qcg1(1,1,2))

  do k=2,nz-1
    kp=k+1

    xyp=(dpt(1,1,kp)+ta(kp))*(1.+.61*(dqv(1,1,kp)+qa(kp))      &
    -qcl1(1,1,kp)-qrn1(1,1,kp)-qci1(1,1,kp)-qcs1(1,1,kp)-qcg1(1,1,kp))
    if(xyp-xym .ge. 0.0) then
      xxw(1,1,k)=0.0
    else
      xxw(1,1,k)=0.5
    endif
    if(k.ne.nz-1) then
      xym=xy
      xy=xyp
    endif

  enddo

  df0(1,1,2)=0.0
  dfc(1,1,2)=0.0
  xxw(1,1,1)=xxw(1,1,2)
  xxw(1,1,nz)=xxw(1,1,nz-1)

  do k=3,nz
    km=k-1
    y1k=2.17375*y1(k)
    qak=qa(k)+qa(k-1)
    tak=ta(k)+ta(k-1)
    ta1k=ta(k)-ta(k-1)
    qa1k=qa(k)-qa(k-1)

    df0(1,1,k)=y1k*(dqv(1,1,k)+dqv(1,1,km)+qak)                    &
    /(dpt(1,1,k)+dpt(1,1,km)+tak)**2

    dfc(1,1,k)=-am1(k)*(ak2(1,1,k)*(1.+xxw(1,1,k))+ak2(1,1,km)    &
      *(1.+xxw(1,1,km)))*(dpt(1,1,k)-dpt(1,1,km)+ta1k             &
      +y2(k)*(dqv(1,1,k)-dqv(1,1,km)+qa1k))                       &
      /((1.+df0(1,1,k)*y2(k))*dz)

    dfc(1,1,k)=dfc(1,1,k)*df0(1,1,k)

    dfc(1,1,k)=-dfc(1,1,k)-(am1(k)/dz)*(ak2(1,1,k)+ak2(1,1,km))    &
              *(dqv(1,1,k)-dqv(1,1,km)+qa(k)-qa(km)                &
              +qcl1(1,1,k)-qcl1(1,1,km)+qci1(1,1,k)-qci1(1,1,km))

  enddo

  if (print_more) then
    print '(a,2f20.10)','Min/max dqv ',minval(dqv), maxval(dqv)
    print '(a,2f20.10)','Min/max dpt ',minval(dpt), maxval(dpt)
  endif

  ! Do vertical advection for each condensate species (neglect dpt,dqv)

  ! Cloud liquid water
  ismg=3
  ! Do vertical advection
  if (print_more) print '(a)','Calling advection for condensate'
  scr3d = qcl
  call advect(qcl,qcl1,w1,xxw,ak,dfc,rho,rho1,rrho,scr3d,am, am1, vtp, ismg, dz, dt, 1.0)

  if (print_more) then
    print '(a,2f20.10)','Min/max qcl ',minval(qcl), maxval(qcl)
  endif

  ! Pristine ice
  irsg=4
  ismg=4
  ! Compute mass-weighted terminal velocity
!   print '(a)','Calling tervrh for qci'
  call tervrh (dpt,qcl,qrn,qci,qcs,qcg,ta,pi,vtp,irsg,rho,fv,improve)
  ! Do vertical advection
!   print '(a)','Calling advection for qci'
  scr3d = qci
  call advect(qci,qci1,w1,xxw,ak,dfc,rho,rho1,rrho,scr3d,am, am1, vtp, ismg, dz, dt, 1.0)

  if (print_more) then
    print '(a,2f20.10)','ri, ar      ',ri, ar
    print '(a,2f20.10)','qrn, vtp    ',qrn(1,1,2), vtp(1,1,2)
  endif

  if (print_more) then
    print '(a,2f20.10)','Min/max qci ',minval(qci), maxval(qci)
  endif

  ! Rain
  irsg=0
  ! Compute mass-weighted terminal velocity
!   print '(a)','Calling tervrh for qrn'
  call tervrh (dpt,qcl,qrn,qci,qcs,qcg,ta,pi,vtp,irsg,rho,fv,improve)

  ! Compute rainfall constants
!   a1 = 36.e3 * rho_int(2,t)
  a1 = 36.e3 * rho1(2) !  = 3600 * 10mm/cm * surface layer density
  a2 = 2.77778e-4 * dt !  = dt/3600 = timestep length in hours

  ! Compute rainfall
!   do j=1,1
!     do i=1,1
!       ri(i,j)=a1 * qrn(i,j,2) * vtp(i,j,2)
!       ar(i,j)=ar(i,j)+a2*ri(i,j)
!     enddo
!   enddo
  ri = a1 * qrn(1,1,2) * vtp(1,1,2) ! ri is in units of mm/hour due to conversion via multiplication with a1
  ar = ar + a2 * ri ! ar is in units of mm

  if (print_more) then
    print '(a,2f20.10)','ri, ar      ',ri, ar
    print '(a,2f20.10)','qrn, vtp    ',qrn(1,1,2), vtp(1,1,2)
  endif

  ismg=5
  ! Do vertical advection
!   print '(a)','Calling advection for qrn'
  scr3d(1,1,:) = 0.0
  call advect(qrn,qrn1,w1,xxw,ak,dfc,rho,rho1,rrho,scr3d,am, am1, vtp, ismg, dz, dt, 0.0)

  if (print_more) then
    print '(a,2f20.10)','Min/max qrn ',minval(qrn), maxval(qrn)
  endif

  ! Snow
  irsg=1
  ! Compute mass-weighted terminal velocity
!   print '(a)','Calling tervrh for qcs'
  call tervrh (dpt,qcl,qrn,qci,qcs,qcg,ta,pi,vtp,irsg,rho,fv,improve)
  ! Do vertical advection
!   print '(a)','Calling advection for qcs'
  scr3d(1,1,:) = 0.0
  call advect(qcs,qcs1,w1,xxw,ak,dfc,rho,rho1,rrho,scr3d,am, am1, vtp, ismg, dz, dt, 0.0)

  if (print_more) then
    print '(a,2f20.10)','Min/max qcs ',minval(qcs), maxval(qcs)
  endif

  ! Graupel
  irsg=2
  ! Compute mass-weighted terminal velocity
!   print '(a)','Calling tervrh for qcg'
  call tervrh (dpt,qcl,qrn,qci,qcs,qcg,ta,pi,vtp,irsg,rho,fv,improve)
  ! Do vertical advection
!   print '(a)','Calling advection for qcg'
  scr3d(1,1,:) = 0.0
  call advect(qcg,qcg1,w1,xxw,ak,dfc,rho,rho1,rrho,scr3d,am, am1, vtp, ismg, dz, dt, 0.0)

  if (print_more) then
    print '(a,2f20.10)','Min/max qcg ',minval(qcg), maxval(qcg)
  endif

!   print '(a,2(2x,e10.5))','a1, a2 ',a1, a2
!   do k=1,nz
!     print '(a,i5,9(2x,e10.5))', &
!           'k, dpt, ta, dqv, qa, y1, pi, f0, dz0, df0 ', &
!            k,dpt(1,1,k),ta(k),dqv(1,1,k),qa(k),y1(k), &
!            pi(k),f0(k),dz0(k),df0(1,1,k)
!   enddo

  ! Save theta, qv before call to saticerh
  dpt_save = dpt
  dqv_save = dqv

!   print '(a)','Calling saticerh'
  call saticerh (dt,dpt,dqv,qcl,qrn,qci,qcs,qcg,w1, &
                 rho,rrho,ta,qa,pi,p0,fv,improve, &
                 tair_out, tairc_out, theta_out, &
                 delta2, delta3, delta4)

  if (print_more) then
    print '(a)','After call to saticerh '
    print '(a,2f20.10)','Min/max dpt ',minval(dpt), maxval(dpt)
    print '(a,2f20.10)','Min/max dqv ',minval(dqv), maxval(dqv)
    print '(a,2f20.10)','Min/max qcl ',minval(qcl), maxval(qcl)
    print '(a,2f20.10)','Min/max qci ',minval(qci), maxval(qci)
    print '(a,2f20.10)','Min/max qrn ',minval(qrn), maxval(qrn)
    print '(a,2f20.10)','Min/max qcs ',minval(qcs), maxval(qcs)
    print '(a,2f20.10)','Min/max qcg ',minval(qcg), maxval(qcg)
    print '(a,2f20.10)','Min/max ri  ',ri,ar
    print '(a)',' '
  endif

!   print '(a,2(2x,e10.5))','ri, ar ',ri(1,1), ar(1,1)
!   print '(a)', 'k, dpt, ta, dqv, qa, qcl, qrn, qci, qcs, qcg '
!   do k=1,nz
!     print '(i5,9(2x,e10.5))', &
!           k,dpt(1,1,k),ta(k),dqv(1,1,k),qa(k), &
!            qcl(1,1,k),qrn(1,1,k),qci(1,1,k),qcs(1,1,k),qcg(1,1,k)
!   enddo

! Save output if at the proper time
  if ( (real(t)*dt)/out_dt - (t*int(dt))/int(out_dt) .eq. 0. ) then
    if (print_more) print '(a,i10,F20.10)','Saving output at timestep, seconds: ',t,t*dt
    ! Make sure we do not exceed dimensions of output arrays...
    if ( tt .lt. nt_out ) tt = tt + 1 ! Increment index in output arrays
    do k = 1, nz
      ta_full(k,tt)        = ta(k)
      qa_full(k,tt)        = qa(k)
      pi_full(k,tt)        = pi(k)
      p0_full(k,tt)        = p0(k)
      p00_full(k,tt)       = p00(k)
      rho_full(k,tt)       = rho(k)
      dpt_full_micro(k,tt) = dpt(1,1,k) - dpt_save(1,1,k) ! potential temp perturbation in this timestep due only to microphysics
      dpt_full(k,tt)       = dpt(1,1,k)
      dqv_full(k,tt)       = dqv(1,1,k)
      qv_sorc_full(k,tt)   = qv_sorc(k)
      qv_sorc2_full(k,tt)  = qv_sorc2(k)
      w_full(k,tt)         = w1(1,1,k)
      qcl_full(k,tt)       = qcl(1,1,k)
      qrn_full(k,tt)       = qrn(1,1,k)
      qci_full(k,tt)       = qci(1,1,k)
      qcs_full(k,tt)       = qcs(1,1,k)
      qcg_full(k,tt)       = qcg(1,1,k)
      rlw_full(k,tt)       = rlw(1,1,k)
      rsw_full(k,tt)       = rsw(1,1,k)
      tair_full(k,tt)      = tair_out(1,1,k)
      tairc_full(k,tt)     = tairc_out(1,1,k)
      theta_full(k,tt)     = theta_out(1,1,k)
      delta2_full(k,tt)    = delta2(1,1,k)
      delta3_full(k,tt)    = delta3(1,1,k)
      delta4_full(k,tt)    = delta4(1,1,k)
    enddo
    ri_full(tt) = ri
    ar_full(tt) = ar

    ! Call observation simulation routine
    call simobs(nz, twcz, pi, dpt, ta, dqv, qa, &
                qcl, qrn, qci, qcs, qcg, pwv, lwp, iwp, ctop)
    pwv_full(tt)   = pwv
    lwp_full(tt)   = lwp
    iwp_full(tt)   = iwp

  endif

  ! Restore potential temperature perturbation
  dpt = dpt_save

  ! If we have reached a radiation timestep, call radiation (do this every time for averaged obs)
  if ( (real(t)*dt)/out_dt - (t*int(dt))/int(out_dt) .eq. 0. .or. &
       real(t)/real(rad_freq) - t/rad_freq .eq. 0. .or. &
       nint ( ( real(t) * dt ) / 60. )  .eq. nint(obs_times(ob)) ) then 

    if ( print_more ) print '(a,i10)','Calling radiation at timestep: ',t

    ! Compute zenith angle
    call zang (rmonth,riday,cosz,rlat,hrl)
    if ( print_more ) print '(a,f10.5)','Zenith angle: ',acos(cosz)*180./3.1415

    ! Save theta before call to radiation
    dpt_save = dpt

    ! Choose a radiation scheme
    if ( new_rad .eq. 1 ) then ! Call new (Goddard) radiation

      call update_radiation (rmonth,riday,rlat,cosz,ta,qa,rho,p0,p00,    &
                            pi,dz0,tland,qcl,qrn,qci,qcs,qcg,dpt,dqv,   &
                            rcp,rc2p,rrp,rpp,rsp,rap,rgp,rhp,           &
                            rsw,rlw,rflux,rams_micro)

    else ! call old radiation

      ! Call pradrat
!       print '(a)','Calling pradrat'
      iflag = 1
      call pradrat(iflag,iopcloud,iradave,cosz,ta,qa,rho,p0,p00,pi,dz0,&
                    tland,qcl,qrn,qci,qcs,qcg,dpt,dqv,rsw,rlw,rflux,improve)

      if ( print_more ) print*,'After pradrat, line 673, rflux = ',rflux

    endif

    ! Save output
    if ( (real(t)*dt)/out_dt - (t*int(dt))/int(out_dt) .eq. 0. ) then
      do k = 1, nz
        rsw_full(k,tt) = rsw(1,1,k)
        rlw_full(k,tt) = rlw(1,1,k)
        dpt_full_rad(k,tt) = dpt(1,1,k) - dpt_save(1,1,k)
      enddo
      ! TOA LW and SW
      lwrad_full(tt) = abs(rflux(1,1,5))
      swrad_full(tt) = abs(rflux(1,1,4))
    endif

    if (print_more) then
      print '(a,2f20.10)','lwrad/swrad    ',abs(rflux(1,1,5)), abs(rflux(1,1,4))
      print '(a,8f20.10)','rflux          ',rflux(1,1,1:8)
    endif

    ! Restore potential temperature perturbation
    dpt = dpt_save

  endif ! End if-statement for whether to call radiation

  ! Fill observation array if not storing averaged obs
  if ( .not. average_obs ) then
    if ( nint ( ( real(t) * dt ) / 60. )  .eq. nint(obs_times(ob)) ) then
      if ( print_more ) then
        print '(a,i5,2(1x,i10))','Timestep, time, obs time: ',t, nint ( ( real(t) * dt ) / 60. ), nint(obs_times(ob))
      endif
      call simobs(nz, twcz, pi, dpt, ta, dqv, qa, &
                  qcl, qrn, qci, qcs, qcg, pwv, lwp, iwp, ctop)
      obs_array(1,ob) = ri           ! (1) Precip rate
      obs_array(2,ob) = ar           ! (2) Accumulated precip
      obs_array(3,ob) = lwp          ! (3) LWP
      obs_array(4,ob) = iwp          ! (4) IWP
      obs_array(5,ob) = abs(rflux(1,1,5)) ! (5) LWrad
      obs_array(6,ob) = abs(rflux(1,1,4)) ! (6) SWrad
      ob = ob + 1 ! Increment observation time index
    endif
 
  endif

  ! Update local time for solar zenith angle computation
!   hrl = hrl + dt / 3600.

  ! Update microphysics variables
  qcl1 = qcl
  qrn1 = qrn
  qci1 = qci
  qcs1 = qcs
  qcg1 = qcg

!  stop 'Check model'

  ! Compute diffusion coefficient at the very end
  call akcoef (ak,ak1,ak_f,dpt,ta,dqv,qa,pi,qcl1,qci1, &
               rho1,rrho,am,am1,w1,dz,dt)
  call satdt_ak (ak,ak1,ak_f,dt)
  call two_ak (ak,ak2)

  ! Update vertical profiles of w and qv'
  ! Update max values, height of max value and half-width
  w_max   = w_max  + w_max_inc
  q_max   = q_max  + q_max_inc
  z_w     = z_w    + z_w_inc
  z_q     = z_q    + z_q_inc
  half_w  = half_w + half_w_inc
  half_q  = half_q + half_q_inc

  q_max2  = q_max2  + q_max_inc2
  z_q2    = z_q2    + z_q_inc2
  half_q2 = half_q2 + half_q_inc2

  ! Recompute forcing profiles
  w(1)        = 0.0
  qv_sorc(1)  = 0.0
  qv_sorc2(1) = 0.0
  do k = 2, nz-1
    if ( abs(z_lay(k) - z_w) <= half_w ) then
      if ( case_type .eq. 1 ) then ! Convective case uses cos^4
        w(k) = w_max * ( cos ( half_pi * ( z_lay(k) - z_w ) / half_w ) )**4
      else if ( case_type .eq. 2 ) then ! Frontal case uses cos^6
        w(k) = w_max * ( cos ( half_pi * ( z_lay(k) - z_w ) / half_w ) )**6
      endif
    else
      w(k)      = 0.0
    endif
    if ( abs(z_lay(k) - z_q) <= half_q ) then
      if ( case_type .eq. 1 ) then ! Convective case uses cos^2
        qv_sorc(k)  = q_max  * ( cos ( half_pi * ( z_lay(k) - z_q  ) / half_q  ) )**2   &
                  / ( 1.e3 * tend_dt ) ! Convert to kg/kg/s
      else if ( case_type .eq. 2 ) then ! Frontal case uses cos^6
        qv_sorc(k)  = q_max  * ( cos ( half_pi * ( z_lay(k) - z_q  ) / half_q  ) )**6   &
                  / ( 1.e3 * tend_dt ) ! Convert to kg/kg/s
      endif
    else
      qv_sorc(k)  = 0.0
    endif
    if ( abs(z_lay(k) - z_q2) <= half_q2 ) then
      if ( case_type .eq. 1 ) then ! Convective case uses cos^2
        qv_sorc2(k) = q_max2 * ( cos ( half_pi * ( z_lay(k) - z_q2 ) / half_q2 ) )**2   &
                  / ( 1.e3 * tend_dt ) ! Convert to kg/kg/s
        if ( neg_2nd_f ) qv_sorc2(k) = qv_sorc2(k) * (-1.0)
      else if ( case_type .eq. 2 ) then ! Frontal case uses cos^6
        qv_sorc2(k) = q_max2 * ( cos ( half_pi * ( z_lay(k) - z_q2 ) / half_q2 ) )**6   &
                  / ( 1.e3 * tend_dt ) ! Convert to kg/kg/s
        if ( neg_2nd_f ) qv_sorc2(k) = qv_sorc2(k) * (-1.0)
      endif
    else
      qv_sorc2(k) = 0.0
    endif
  enddo

  ! Update GCE vertical velocity
  do k = 1, nz
    w1(1,1,k) = w(k) * 100.        ! Convert from m/s to cm/s
  enddo

  ! If we have hit an even forcing time, update indices and increments
  if ( (real(t)*dt) / tend_dt .eq. real( (t*int(dt)) / int(tend_dt) ) &
      .and. t .gt. 2 .and. t .lt. nt-2 ) then

    ! Update indices on forcing
    f1 = f1 + 1
    f2 = f2 + 1

    if (print_more) print '(a,i10)', 'Updating forcing at timestep ',t
    if (print_more) print '(a,2i5)', 'Forcing indices: ',f1,f2

    ! Set increments on the above 4 variables
    w_max_inc = ( w_max_f(f2) - w_max_f(f1) ) / (tend_dt / dt + 1.)
    q_max_inc = ( q_max_f(f2) - q_max_f(f1) ) / (tend_dt / dt + 1.)
    z_w_inc = ( z_w_f(f2) - z_w_f(f1) ) / (tend_dt / dt + 1.) * 100.e3 ! convert from km to cm
    z_q_inc = ( z_q_f(f2) - z_q_f(f1) ) / (tend_dt / dt + 1.) * 100.e3 ! convert from km to cm
    half_w_inc = ( half_w_f(f2) - half_w_f(f1) ) / (tend_dt / dt + 1.) * 100.e3 ! convert from km to cm
    half_q_inc = ( half_q_f(f2) - half_q_f(f1) ) / (tend_dt / dt + 1.) * 100.e3 ! convert from km to cm

    q_max_inc2  = ( q_max_f2(f2)  - q_max_f2(f1)  ) / (tend_dt / dt + 1.)
    z_q_inc2    = ( z_q_f2(f2)    - z_q_f2(f1)    ) / (tend_dt / dt + 1.) * 100.e3 ! convert from km to cm
    half_q_inc2 = ( half_q_f2(f2) - half_q_f2(f1) ) / (tend_dt / dt + 1.) * 100.e3 ! convert from km to cm

  endif

enddo ! End loop over timesteps


! If called for, fill observation array with averaged obs
if ( average_obs ) then
  ! Loop over each obs interval
  do ob = 1, n_obs_times

    ! Find start and indices corresponding to obs averaging interval
    ibeg = nint ( ( obs_start_time(ob) * 60.0 ) / out_dt )
    iend = nint ( ( obs_end_time(ob)   * 60.0 ) / out_dt )
    ! Make sure indices do not exceed bounds
    if ( ibeg .le. 0  ) ibeg = 1
    if ( iend .gt. nt / (out_dt / dt) ) iend = nint(nt / (out_dt / dt))
    nt_obs = real( iend - ibeg + 1 )

    if ( print_more ) then
      print '(a,2(1x,F6.1),2(1x,i10),f6.1)','Beg time, end time, beg index, end index, nt: ', &
        obs_start_time(ob), obs_end_time(ob), ibeg, iend, nt_obs
    endif
    obs_array(1,ob) = sum(ri_full (ibeg:iend)) / nt_obs          ! (1) Precip rate
    obs_array(2,ob) = sum(ar_full (ibeg:iend)) / nt_obs          ! (2) Accumulated precip
    obs_array(3,ob) = sum(lwp_full(ibeg:iend)) / nt_obs          ! (3) LWP
    obs_array(4,ob) = sum(iwp_full(ibeg:iend)) / nt_obs          ! (4) IWP
    obs_array(5,ob) = sum(lwrad_full(ibeg:iend)) / nt_obs        ! (5) LWrad
    obs_array(6,ob) = sum(swrad_full(ibeg:iend)) / nt_obs        ! (6) SWrad
  enddo
endif


! Write output to file
if ( write_output ) then

  open ( 300, file='dpt_micro_out.dat', form='unformatted', action='write')
  open ( 301, file='dpt_rad_out.dat',   form='unformatted', action='write')
  open ( 302, file='dqv_out.dat', form='unformatted', action='write')
  open ( 303, file='qv_sorc_out.dat', form='unformatted', action='write')
  open ( 304, file='w_out.dat', form='unformatted', action='write')
  open ( 305, file='qcl_out.dat', form='unformatted', action='write')
  open ( 306, file='qrn_out.dat', form='unformatted', action='write')
  open ( 307, file='qci_out.dat', form='unformatted', action='write')
  open ( 308, file='qcs_out.dat', form='unformatted', action='write')
  open ( 309, file='qcg_out.dat', form='unformatted', action='write')
  open ( 310, file='rsw_out.dat', form='unformatted', action='write')
  open ( 311, file='rlw_out.dat', form='unformatted', action='write')
  open ( 312, file='ri_out.dat',  form='unformatted', action='write')
  open ( 313, file='ar_out.dat',  form='unformatted', action='write')
  open ( 314, file='pwv_out.dat', form='unformatted', action='write')
  open ( 315, file='lwp_out.dat', form='unformatted', action='write')
  open ( 316, file='iwp_out.dat', form='unformatted', action='write')
  open ( 317, file='lwrad_out.dat', form='unformatted', action='write')
  open ( 318, file='swrad_out.dat', form='unformatted', action='write')
  open ( 319, file='ta_out.dat', form='unformatted', action='write')
  open ( 320, file='qa_out.dat', form='unformatted', action='write')
  open ( 321, file='pi_out.dat', form='unformatted', action='write')
  open ( 322, file='p0_out.dat', form='unformatted', action='write')
  open ( 323, file='p00_out.dat', form='unformatted', action='write')
  open ( 324, file='rho_out.dat', form='unformatted', action='write')
  open ( 325, file='dpt_out.dat', form='unformatted', action='write')
  open ( 326, file='tair_out.dat', form='unformatted', action='write')
  open ( 327, file='tairc_out.dat', form='unformatted', action='write')
  open ( 328, file='theta_out.dat', form='unformatted', action='write')
  open ( 329, file='delta2_out.dat', form='unformatted', action='write')
  open ( 330, file='delta3_out.dat', form='unformatted', action='write')
  open ( 331, file='delta4_out.dat', form='unformatted', action='write')
  open ( 332, file='qv_sorc2_out.dat', form='unformatted', action='write')

  write ( 300 ) dpt_full_micro
  write ( 301 ) dpt_full_rad
  write ( 302 ) dqv_full
  write ( 303 ) qv_sorc_full
  write ( 304 ) w_full
  write ( 305 ) qcl_full
  write ( 306 ) qrn_full
  write ( 307 ) qci_full
  write ( 308 ) qcs_full
  write ( 309 ) qcg_full
  write ( 310 ) rsw_full
  write ( 311 ) rlw_full
  write ( 312 ) ri_full
  write ( 313 ) ar_full
  write ( 314 ) pwv_full
  write ( 315 ) lwp_full
  write ( 316 ) iwp_full
  write ( 317 ) lwrad_full
  write ( 318 ) swrad_full
  write ( 319 ) ta_full
  write ( 320 ) qa_full
  write ( 321 ) pi_full
  write ( 322 ) p0_full
  write ( 323 ) p00_full
  write ( 324 ) rho_full
  write ( 325 ) dpt_full
  write ( 326 ) tair_full
  write ( 327 ) tairc_full
  write ( 328 ) theta_full
  write ( 329 ) delta2_full
  write ( 330 ) delta3_full
  write ( 331 ) delta4_full
  write ( 332 ) qv_sorc2_full

  do i = 300, 332
    close (i)
  enddo

endif

return
end subroutine drive_micro_rad
