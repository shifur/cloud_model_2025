! This FORTRAN 90 subroutine runs an ensemble of model forecasts for a 
! single-column version of the NASA Goddard Cumulus Ensemble (GCE) model. 
! The microphysics, radiation and vertical advection have been stripped 
! from the GCE, and are driven by (time-varying) vertical motion and 
! water vapor increments. This column version of the GCE is based on the 
! unified 3D version of the model current as of 3 November 2008. It 
! includes the Chou and Suarez radiation code and the 3ice microphysics 
! with updates by Steve Lang (GSFC).
!
! The ensemble consists of a set of cloud microphysical parameters read 
! in from file.
!
! Setup is similar to Ben Shipway's idealized CRM tests, which are 
! documented more fully at:
!
!         http://www.convection.info/microphysics/
!
! Note that this routine assumes that the user has already prepared 
! all of the nessary inputs to the GCE model by calling "ideal_setup" in
! the main driver routine "crm_model_error".
!
! Procedure:
! 1.  Open the file containing parameters and read in the set
! 2.  Check parameters vs. the allowable range
! 3.  Run simulations for each set of parameters storing output in output files
!
! NOTE: The ensemble file is assumed to contain the parameters to be run, one experiment per line (row)
! Each row is assumed to contain the full list of parameters. 
! 
! INPUT VARIABLES
!
! ENSEMBLE VARIABLES
!
! nmembers        Number of ensemble members (integer)
! ens_file        Name of the file containing ensemble of parameters (character)
! print_less      Logical, whether to print a little diagnostic output
! print_more      Logical, whether to print a lot of diagnostic output
! exp_name        Character, prefix for output file names
!
! VARIABLES NEEDED TO RUN THE COLUMN MODEL
!
! OUTPUT VARIABLES
!
! NONE
!
! External:
!
!   subroutine drive_micro_rad    Driver for column GCE model, uses:
!
!     microphysics.a              library containing 3ice microphysics routines
!     radiation.a                 library containing radiation routines
!     dynamics.a                  library containing dynamics (advection) routines
!
! In vectors dimensioned (nparams), variables are ordered:
!   as, bs, ag, bg, tnw, tns, tng, roqs, roqg, bnd21, bnd1
!
! In vectors dimensioned (nobs), variables are ordered:
!   Prate  Ptot  LWP  IWP   LW   SW
!
! Derek Posselt
! University of Michigan
! 17 January 2011
!
!-----------------------------------------------------------------------
! 18 March 2013
! D. Posselt
!
! Added the option to use time-averaged observations (average_obs = .true.)
! with start time and end time set by obs_start_time and obs_end_time, 
! in units of seconds.
!
! ----------------------------------------------------------------------
! 20 January 2023
! D. Posselt
!
! The number of ensemble members is now determined from the length of the ens_file
!
! ----------------------------------------------------------------------
subroutine ensemble_fcst( nparams, ens_file, pmin, pmax, &
                          print_less, print_more, exp_name, &
                          max_obs, n_obs_times, obs_times, nobs, &
                          average_obs, obs_start_time, obs_end_time, &
! Variables needed to drive the column model ("drive_micro_rad" routine)
                          nz,dz,nt,nt_out,dt,rad_dt,out_dt,ntend,tend_dt,tsfc, &
                          t_lev, t_lay, p_lev, p_lay, q_lev, q_lay, &
                          rho_lev, rho_lay, w, z_lev, z_lay, &
                          w_max_f, q_max_f, z_w_f, z_q_f, half_w_f, half_q_f, &
                          q_max_f2, z_q_f2, half_q_f2, neg_2nd_f, &
                          start_month, start_day, start_hour, latitude, &
                          write_output, rams_micro, new_rad, case_type, improve, &
                          output_file )
implicit none

!-----------------------------------------------------------------------
!                        VARIABLE DECLARATIONS
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!                         ENSEMBLE VARIABLES
!-----------------------------------------------------------------------

! Input-only variables

character (len=200), intent(in) :: exp_name ! Name of the experiment (for output file names)

logical, intent(in) :: print_less      & ! Write basic output to screen (keeping track of inversion)
                     , print_more        ! Write extended output to screen (model run diagnostics, etc)

character (len=255), intent(in) :: ens_file        ! File containing ensemble of parameters

real,    dimension(nparams), intent(in) :: pmin  & ! Minimum allowable parameter values
                                         , pmax    ! Maximum allowable parameter values

! Parameter values
real,    dimension(nparams)       :: params        ! Parameter set

character (len=255), intent(in) :: output_file    ! Name of output file

!-----------------------------------------------------------------------
!               VARIABLES USED TO DRIVE THE COLUMN MODEL
!-----------------------------------------------------------------------
! Input-only
integer, intent(in) :: nz, nt, nt_out, ntend, nparams, max_obs, n_obs_times, nobs, &
                       rams_micro, new_rad, case_type, improve
real,    intent(in) :: dt, rad_dt, out_dt, tend_dt, tsfc
real,    intent(in) :: start_month, start_day, start_hour, latitude
real                :: dz
real, dimension(nz,nt_out), intent(in) :: &
                       z_lev,   z_lay   & ! Heights       at layer boundaries and midpoints
                     , T_lev,   T_lay   & ! Temperature   at layer boundaries and midpoints
                     , q_lev,   q_lay   & ! Vapor mix rat at layer boundaries and midpoints
                     , p_lev,   p_lay   & ! Pressure      at layer boundaries and midpoints
                     , rho_lev, rho_lay   ! Density       at layer boundaries and midpoints

real, dimension(max_obs), intent(in) :: obs_times   ! Observation times (minutes)

logical, intent(in)                     :: average_obs      ! whether to average observations over "nobs_times" time intervals or not
real,    dimension(nobs), intent(in)    :: obs_start_time & ! Start times for accumulating observations 
                                         , obs_end_time     ! End times for accumulating observations 

real, dimension(nz), intent(in)   :: w              ! Vertical velocity (layer mid-points)

real, dimension(ntend), intent(in) :: w_max_f   & ! Max vertical velocity (m/s)
                                    , q_max_f   & ! Max qv forcing (g/kg/hour)
                                    , z_w_f     & ! Height of max vertical velocity (km)
                                    , z_q_f     & ! Height of max qv forcing (km)
                                    , half_w_f  & ! Half-width of vertical velocity profile (km)
                                    , half_q_f  & ! Half-width of qv forcing profile (km)
                                    , q_max_f2  & ! Max of 2nd qv forcing (g/kg/hour)
                                    , z_q_f2    & ! Height of 2nd max qv forcing (km)
                                    , half_q_f2   ! Half-width of 2nd qv forcing profile (km)

logical, intent(in) :: neg_2nd_f      ! Whether 2nd set of water vapor forcing fields should be negative

! Output from column model at full time resolution (nt_out)
real, dimension(nt_out)       :: ri_full, ar_full &
                               , lwp_full, iwp_full &
                               , lwrad_full, swrad_full

! Whether to write model output to file--should be false for MCMC...
logical, intent(in) :: write_output

! Output from the column model
real, dimension(nobs, n_obs_times) :: fwd_obs_array

!-----------------------------------------------------------------------
!                          LOCAL VARIABLES
!-----------------------------------------------------------------------

! File unit numbers
integer, parameter :: unit_ens      = 300
integer, parameter :: unit_fwd_obs  = 301
integer, parameter :: unit_params   = 302

! Loop indices
integer :: iens, p, t

! Ensemble variables
integer :: nmembers   ! Number of ensemble members

! Flag for negative parameter values
logical :: found_negative

! Capture end of file in the ensemble parameter read
integer :: read_status

!-----------------------------------------------------------------------
!                           INITIAL SETUP
!-----------------------------------------------------------------------

! Check for file containing ensemble of parameters
call check_file ( trim(adjustl(ens_file)),      .true. )

! Check for existence of output files--if any exist, stop the program
! call check_file ( trim(adjustl(output_file)),  .false. )
! call check_file ( trim(adjustl(exp_name))//'_ens_full_params.dat',   .false. )

! Files have passed the check--open each output file for write
! Unformatted Fortran binary
! open ( unit_fwd_obs , file = trim(adjustl(exp_name))//'_ens_fwd_obs.dat', form = 'unformatted', action = 'write', access = 'sequential' )
! open ( unit_params ,  file = trim(adjustl(exp_name))//'_ens_full_params.dat',  form = 'unformatted', action = 'write', access = 'sequential' )
! Text
open ( unit_fwd_obs , file = trim(adjustl(output_file)), form = 'formatted', action = 'write', access = 'sequential' )
! open ( unit_params ,  file = trim(adjustl(exp_name))//'_ens_full_params.dat',  form = 'formatted', action = 'write', access = 'sequential' )

! Get the number of ensemble members from the length of the ensemble input file
open ( unit_ens      , file = trim(adjustl(ens_file)), form = 'formatted', action = 'read', access = 'sequential' )
nmembers = 0 ! Initialize the number of ensemble members
! Loop until we find EOF
DO
  read_status = 0
  read( unit_ens, *, IOSTAT = read_status )
  if ( read_status .ne. 0 ) exit
  nmembers = nmembers + 1
ENDDO
if ( print_less ) print '(a,i10,a)','Found ',nmembers, ' ensemble members in the file'

! Rewind the file back to the beginning
REWIND ( unit_ens )

!-----------------------------------------------------------------------
!                    LOOP OVER ENSEMBLE MEMBERS
!-----------------------------------------------------------------------

do iens = 1, nmembers

  if ( print_less .and. float(iens/50) .eq. float(iens)/50.0 ) &
    print '(a,i10)','Ensemble number ',iens

  ! Read parameters from file
  read_status = 0
  read ( unit_ens, *, IOSTAT = read_status ) params

  ! If read has failed, fill parameters and model output with missing values and move on
  if ( read_status .ne. 0 ) then

    found_negative = .true.
    params(:) = -9999.0
  
  else

    ! Check parameter values vs. min and max
    do p = 1, nparams
      if ( params(p) .lt. pmin(p) .or. params(p) .gt. pmax(p) .and. print_less ) then
        print '(a)','Warning! Parameter outside allowable range'
        print '(a,i10,1x,i3,3(1x,f20.10))','Ensemble member, parameter number, min, value, max: ', &
                                     iens, p, pmin(p), params(p), pmax(p)
  !       stop 'Stopping -- check values in parameter input file'
      endif
    enddo

    ! Check for negative parameter values--if found, skip the model run 
    ! and fill the forward array with -9999.
    found_negative = .false.
    do p = 1, nparams
      if ( params(p) .le. 0.0 ) found_negative = .true.
    enddo

  endif ! End if statement checking for valid read of parameters

  if ( found_negative ) then

    fwd_obs_array = -9999.0

    if ( print_less ) then 
      print '(a)','Found negative parameter in set'
      do p = 1, nparams
        write ( 6,'(1x,f15.6)',advance='no') params(p)
      enddo
      write ( 6, * )
    endif

  else

    ! Run the model with the current set of parameter values
    if ( print_less ) then
      print '(a)','Running model with parameter set'
      do p = 1, nparams
        write ( 6,'(1x,f11.5)',advance='no') params(p)
      enddo
      write ( 6, * )
    endif

    call drive_micro_rad( nz,dz,nt,nt_out,dt,rad_dt,out_dt,ntend,tend_dt,tsfc, &
                          t_lev, t_lay, p_lev, p_lay, q_lev, q_lay, &
                          rho_lev, rho_lay, w, z_lev, z_lay, &
                          w_max_f, q_max_f, z_w_f, z_q_f, half_w_f, half_q_f, &
                          q_max_f2, z_q_f2, half_q_f2, neg_2nd_f, &
                          start_month, start_day, start_hour, latitude, &
                          nparams, params, max_obs, n_obs_times, obs_times, &
                          nobs, fwd_obs_array, write_output, print_more, rams_micro, new_rad, &
                          ri_full, ar_full, lwp_full, iwp_full, lwrad_full, swrad_full, &
                          case_type, improve, &
                          average_obs, obs_start_time, obs_end_time )


    ! Print diagnostics
    if ( print_more ) then
      print '(a,i10)','Output from column model, ensemble member ',iens
      do t = 1, n_obs_times ! (n_obs_times is the number of observation times)
        print '(a,i6,6f12.2)','Time, forward observations ', t, fwd_obs_array(:,t)
      enddo
    endif

  endif ! End check for negative values

  ! Write parameters and forward output to file
  ! write ( unit_params   ) params
  ! write ( unit_fwd_obs  ) fwd_obs_array
  ! Write to formatted text
  ! do p = 1, nparams
  !   write ( unit_params, '(1X,F20.10)', advance='no' ) params(p)
  ! enddo
  ! write ( unit_params, * ) ! End the line
  do t = 1, n_obs_times 
    do p = 1, nobs     
      write ( unit_fwd_obs, '(1X,F20.10)', advance='no'  ) fwd_obs_array(p,t)
    enddo
  enddo
  write ( unit_fwd_obs, * ) ! End the line

enddo ! End loop over ensemble members

! Close file containing input parameters
close ( unit_ens )

! Close output files
close (  unit_fwd_obs  )
! close (  unit_params   )

return
end subroutine ensemble_fcst




