! This F90 program serves as the main driver for the CRM column radiation
! and microphysics model error experiments.
!
! Function
!
! This program is tasked with reading in the forcing fields and interpolating
! them to the model timestep, then calling the microphysics/radiation driver
! for the requested experiment. Possible experiments include:
!
! 1. Test run of the column microphysics--single integration
! 2. Ensemble forecast
!
! Procedure
!
! 1. Read in the namelist
! 2. Read in forcing fields and interpolate in time
! 3. Decide which experiment to perform
! 4. Call the appropriate subroutine
!
! Derek Posselt
! University of Michigan
! 3 October 2008
!
!-----------------------------------------------------------------------
! 12 November 2010
! D. Posselt
! 
! Modified this so that the environment (temperature, water vapor, pressure) varies
! with time. The convective case still has an invariant base state environment,
! but the frontal case varies with time. Forcing fields are contained in the 
! variable "tend_file" and read in and interpolated to model timestep in
! subroutine "read_forcing"
! 
! I have also modified the code to read in a variable (still constant) forcing
! time interval "tend_dt" (previously this was hard-coded to an hour).
! 
!-----------------------------------------------------------------------
! 17 January 2011
! D. Posselt
! 
! Added another case to the mix (exp_type = 4). This case reads in a set
! of parameter values from a file and runs an ensemble of forecasts.
! The ensemble of parameters and output forward observations are saved 
! in the same format as in the MCMC case in files with the same prefix as
! in the namelist variable exp_name
!
! Input required consists of the name of the file containing parameters
! (one set per line in the file), the parameters varied (pflag), and the number
! of ensemble runs.
! 
! ----------------------------------------------------------------------
! 18 March 2013
! D. Posselt
!
! Added the option to use time-averaged observations (average_obs = .true.)
! with start time and end time set by obs_start_time and obs_end_time, 
! in units of seconds.
!
!-----------------------------------------------------------------------
! 20 January 2023
! D. Posselt
!
! The input parameter file name (ens_file), output file name, and namelist file name
! are now all optional command line parameters. They default to:
!
! ens_file    = 'run_one.txt'
! output_file = 'crm_column_model_output.txt'
! nl_file     = 'namelist.f90'
!
! Also, the number of ensemble members is now determined from the length of the ensemble input file
!
!-----------------------------------------------------------------------

program crm_model_error
implicit none

! Parameters
integer, parameter :: data_unit = 100, unit_nml = 200, unit_param = 300
integer, parameter :: ntend = 20    ! Number of tendency times (max number of forcing times)
integer, parameter :: nparams = 11  ! Number of perturbable parameters
integer, parameter :: nobs = 6      ! Number of observations (precip rate, accumulated precip, iwp, lwp, lwrad, swrad)
integer, parameter :: max_obs = 720 ! Maximum number of observation times

! Command line inputs
INTEGER :: n_arguments
character (len=255) :: ens_file    ! Name of the file containing ensemble parameters
character (len=255) :: output_file ! Name of the file containing forward obs
character (len=255) :: nl_file     ! Name of the namelist file

! Namelist variables
! Record 0: Model definition
integer :: exp_type       ! Type of experiment to run (1 = single model run, 2 = pdf map, 3 = mcmc)
integer :: case_type      ! Type of case to run (1 = convection, 2 = frontal)
integer :: improve        ! Version of the Goddard 3-ice microphysics to be used. This is the "improve" flag in the 3-ice microphysics. Should be set to -1, 0, or 3
integer :: rams_micro     ! Whether or not to use RAMS microphysics (1 = yes)
integer :: new_rad        ! Whether or not to use new (Goddard) radiation scheme (1 = yes)
real    :: dt           & ! Model timestep (seconds)
         , rad_dt       & ! Radiation timestep (seconds)
         , out_dt       & ! File output timestep (seconds)
         , dz           & ! Vertical grid spacing (first layer, in meters)
         , ztop         & ! Model top (kilometers)
         , run_time       ! Model output frequency (seconds)
real    :: Tsfc         & ! Fixed model surface temperature (K)
         , qsfc         & ! Surface mass mixing ratio (g/kg)
         , psfc         & ! Surface pressure (hPa)
         , tend_dt        ! Forcing time tendency interval (s)
logical :: neg_2nd_f      ! Whether 2nd set of water vapor forcing fields should be negative
character ( len = 255 ) :: tend_file  ! File containing forcing for T, P, qv, rho--only used in case_type 2
real, dimension(ntend) :: w_max_f   & ! Max vertical velocity (m/s)
                        , q_max_f   & ! Max qv forcing (g/kg/hour)
                        , z_w_f     & ! Height of max vertical velocity (km)
                        , z_q_f     & ! Height of max qv forcing (km)
                        , half_w_f  & ! Half-width of vertical velocity profile (km)
                        , half_q_f  & ! Half-width of qv forcing profile (km)
! Add a second set of water vapor forcing terms--no associated vertical velocity for now
                        , q_max_f2  & ! Max 2nd qv forcing (g/kg/hour)
                        , z_q_f2    & ! Height of max 2nd qv forcing (km)
                        , half_q_f2   ! Half-width of 2nd qv forcing profile (km)
real    :: start_month  & ! Month of model start time (imonth)
         , start_day    & ! Day of the month of start time (iday)
         , start_hour   & ! Local time (hour) of model start (hrl)
         , latitude       ! Latitude of model domain
logical :: write_output   ! Whether to write model output to binary files
logical :: print_less     ! Whether to produce more diagnostic output
logical :: print_more     ! Whether to produce a ridiculous amount of diagnostic output
character (len=200) :: exp_name ! Name of the experiment (for output file names)

! Record 1: Parameter bounds settings
real,    dimension(nparams) :: pmin          & ! Minimum parameter values
                             , pmax            ! Maximum parameter values

! Record 2: Observations
integer                     :: n_obs_times    ! Number of observation times
real, dimension(max_obs)    :: obs_times      ! Observation times (in minutes after model start)
logical                     :: average_obs    ! whether to average observations over "nobs_times" time intervals or not
real,    dimension(max_obs) :: obs_start_time & ! Start times for accumulating observations 
                             , obs_end_time     ! End times for accumulating observations 

! Record 3: Parameters and MCMC settings
! integer, dimension(nparams) :: pflag      ! Flag for whether to perturb parameters (1=pert, 0=set to "true")
real, dimension(nparams)    :: ptrue      ! "True" values of parameters

! Data used to drive column model
! Height levels and number of timesteps
integer :: nz, nt, nt_out
real, allocatable, dimension(:,:)   :: z_lev,   z_lay   & ! Heights       at layer boundaries and midpoints
                                     , T_lev,   T_lay   & ! Temperature   at layer boundaries and midpoints
                                     , q_lev,   q_lay   & ! Vapor mix rat at layer boundaries and midpoints
                                     , p_lev,   p_lay   & ! Pressure      at layer boundaries and midpoints
                                     , rho_lev, rho_lay   ! Density       at layer boundaries and midpoints

real, allocatable, dimension(:,:)   :: w                  ! Vertical velocity (layer mid-points)

real, dimension(nparams) :: params_in ! Microphysics parameters

real, allocatable, dimension(:,:) :: obs_array        ! Array (nobs, n_obs_times) containing observations

! Output from column model at full time resolution (nt_out)
real, allocatable, dimension(:)   :: ri_full, ar_full &
                                   , lwp_full, iwp_full &
                                   , lwrad_full, swrad_full

! Miscellaneous variables
logical :: file_exists = .false.
integer :: i, j, k, t, tt, ist, ien, ob, p

! Declare namelist records
namelist /record0/ exp_type, case_type, improve, rams_micro, new_rad, &
                   dt, rad_dt, out_dt, dz, ztop, run_time, &
                   Tsfc, qsfc, psfc, tend_dt, tend_file, neg_2nd_f, &
                   w_max_f, q_max_f, z_w_f, z_q_f, half_w_f, half_q_f, &
                   q_max_f2, z_q_f2, half_q_f2, &
                   start_month, start_day, start_hour, latitude, &
                   write_output, print_less, print_more, exp_name
namelist /record1/ pmin, pmax
namelist /record2/ n_obs_times, obs_times, average_obs, obs_start_time, obs_end_time
namelist /record3/ ptrue

! Read command line arguments, if provided
         
! Check for command line arguments
n_arguments = iargc()     

! Set default values    
ens_file    = 'run_one.txt'
output_file = 'crm_column_model_output.txt'
nl_file     = 'namelist.f90'

! Loop over the number of command line arguments
DO i = 1, n_arguments

  ! (1) If the user has provided the name of the param input file, use it
  IF ( i .eq. 1 ) CALL getarg(i,ens_file)

  ! (2) If the user has provided the name of the model output file, use it
  IF ( i .eq. 2 ) CALL getarg(i,output_file)

  ! (3) If the user has provided the name of the namelist file, use it
  IF ( i .eq. 3 ) CALL getarg(i,nl_file)

ENDDO

! Read in namelist
file_exists = .false.
inquire ( file = trim(adjustl(nl_file)), exist = file_exists )
if ( file_exists ) then
  open ( unit_nml, file = trim(adjustl(nl_file)), form = 'formatted', &
         action = 'read', access = 'sequential', status = 'old' )
  read  ( unit_nml, nml = record0 )
  if (print_less) write ( 6, nml = record0 )
  read  ( unit_nml, nml = record1 )
  if (print_less) write ( 6, nml = record1 )
  read  ( unit_nml, nml = record2 )
  if (print_less) write ( 6, nml = record2 )
  read  ( unit_nml, nml = record3 )
  if (print_less) write ( 6, nml = record3 )
else
  print '(a)','Could not find namelist file'
  stop 'Check for existence of namelist'
endif

! Check compatibility between microphysics and radiation--RAMS micro MUST be
! called with new radiation scheme...
if ( rams_micro .eq. 1 .and. new_rad .ne. 1 ) then
  print '(a)', 'Found incompatible microphysics and radiation: RAMS micro with old radiation'
  print '(a)', 'RAMS micro must be called with new_rad = 1'
  stop 'Stopping--check microphysics and radiation flags'
endif

! Zero obs times after n_obs_times
obs_times(n_obs_times+1:max_obs) = 0.0
obs_start_time(n_obs_times+1:max_obs) = 0.0
obs_end_time(n_obs_times+1:max_obs) = 0.0

! Convert input to SI units
qsfc = qsfc * 1.0e-3
psfc = psfc * 100.0

! Compute number of vertical levels
nz = ((ztop * 1000.) / dz ) + 2
if (print_less) print '(a,i10)','Number of vertical levels: ',nz

! Based on the desired timestep, compute number of output times
nt     = ( run_time * 3600. ) / dt + 1     ! Number of model timesteps
nt_out = ( run_time * 3600. ) / out_dt + 1 ! Number of model output times

! Allocate full obs output vectors
allocate ( ri_full    ( nt_out ) )
allocate ( ar_full    ( nt_out ) )
allocate ( lwp_full   ( nt_out ) )
allocate ( iwp_full   ( nt_out ) )
allocate ( lwrad_full ( nt_out ) )
allocate ( swrad_full ( nt_out ) )

if (print_less) print '(a,i10)','Number of timesteps:        ',nt
if (print_less) print '(a,i10)','Number of output timesteps: ',nt_out

! Allocate column arrays
allocate ( z_lev   (nz, nt) )
allocate ( z_lay   (nz, nt) )
allocate ( T_lev   (nz, nt) )
allocate ( T_lay   (nz, nt) )
allocate ( q_lev   (nz, nt) )
allocate ( q_lay   (nz, nt) )
allocate ( p_lev   (nz, nt) )
allocate ( p_lay   (nz, nt) )
allocate ( rho_lev (nz, nt) )
allocate ( rho_lay (nz, nt) )
allocate ( w       (nz, nt) )

! Allocate observation array
allocate ( obs_array (nobs, n_obs_times) )

! Call the routine that computes or reads in the vertical profiles
if ( case_type .eq. 1 ) then ! Convective case
  ! Convection case
  call ideal_setup (Tsfc,qsfc,psfc,z_lev, z_lay, T_lev, T_lay, &
                    q_lev, q_lay, p_lev, p_lay, rho_lev, rho_lay,    &
                    print_more, dz, nz, nt)

else if ( case_type .eq. 2 ) then 
  ! Frontal case
  call read_forcing ( tend_file, nz, nt, dt, tend_dt, run_time, &
                      z_lev, z_lay, T_lev, T_lay, &
                      q_lev, q_lay, p_lev, p_lay, rho_lev, rho_lay,    &
                      print_more )

  ! If we are running the front case, convert vertical velocity from cm/s to m/s
  w_max_f = w_max_f / 100.

else
! Unrecognized case type

  print '(a)', 'Unrecognized case_type--please check namelist'
  stop 'Stopping'

endif ! End if statement for whether to compute profiles or to read them in

! Write vector of true parameters out to file
! open ( unit_param , file = trim(adjustl(exp_name))//'_param_true.dat', form = 'unformatted', action = 'write', access = 'sequential' )
! write ( unit_param ) ptrue
! close ( unit_param )

!-----------------------------------------------------------------------
!                             EXPERIMENT TYPE 1
!                     SINGLE RUN OF THE COLUMN MODEL
!-----------------------------------------------------------------------

if ( exp_type .eq. 1 ) then

  ! For now, set parameters = true
  params_in = ptrue

  print '(a)','Calling microphysics and radiation driver'
  call drive_micro_rad(nz,dz,nt,nt_out,dt,rad_dt,out_dt,ntend,tend_dt,tsfc, &
                      t_lev, t_lay, p_lev, p_lay, q_lev, q_lay, &
                      rho_lev, rho_lay, w, z_lev, z_lay, &
                      w_max_f, q_max_f, z_w_f, z_q_f, half_w_f, half_q_f, &
                      q_max_f2, z_q_f2, half_q_f2, neg_2nd_f, &
                      start_month, start_day, start_hour, latitude, &
                      nparams, params_in, max_obs, n_obs_times, obs_times, &
                      nobs, obs_array, write_output, print_more, rams_micro, new_rad, &
                      ri_full, ar_full, lwp_full, iwp_full, lwrad_full, swrad_full, &
                      case_type, improve, &
                      average_obs, obs_start_time, obs_end_time, & 
                      output_file )

  print '(a)','Simulated observations '
  do t = 1, n_obs_times
    write (6, '(a,f6.0,a)', advance='no') 'Time ',obs_times(t),' obs '
    do ob = 1, nobs
      write (6, '(f20.10)', advance='no')  obs_array(ob,t)
    enddo
    write ( 6, * ) 
  enddo

!-----------------------------------------------------------------------
!                             EXPERIMENT TYPE 2
!-----------------------------------------------------------------------

elseif ( exp_type .eq. 2 ) then

  if (print_less) print '(a)','Exp type 2 not supported'

!-----------------------------------------------------------------------
!                             EXPERIMENT TYPE 3
!-----------------------------------------------------------------------

elseif ( exp_type .eq. 3 ) then

  if (print_less) print '(a)','Exp type 3 not supported'

!-----------------------------------------------------------------------
!                             EXPERIMENT TYPE 4
!                 RUN ENSEMBLE OF COLUMN MODEL SIMULATIONS
!                  TO GENERATE ENSEMBLE OF FORWARD OBS
!-----------------------------------------------------------------------

elseif ( exp_type .eq. 4 ) then

  ! Run an ensemble of model forecasts
  if (print_less) print '(a)','Calling ensemble forecast routine'
call ensemble_fcst  ( nparams, ens_file, pmin, pmax, &
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

endif

! Check to make sure output file does not yet exist
!   call check_file ( trim(adjustl(exp_name))//'_param_true.dat', .false. )

! stop 'End program crm_model_error'
end program crm_model_error
