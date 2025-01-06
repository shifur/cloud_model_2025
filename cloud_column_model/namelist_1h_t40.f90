! Namelist for the CRM model error experiments

! Model run parameters
&RECORD0
  exp_type   =       4          ! Type of experiment to run (1 = single model run, 4 = ensemble)
  case_type  =       1          ! Type of case to run (1 = convection, 2 = frontal)
  improve    =       3          ! Version of the Goddard 3-ice microphysics to be used. This is the "improve" flag in the 3-ice microphysics. Should be set to -1, 0, or 3
  rams_micro =       0          ! Whether or not to use RAMS microphysics (1 = yes, 0 = use 3ice)
  new_rad    =       0          ! Whether or not to use new (Goddard) radiation (1 = yes, 0 = use old)
  dt         =       5.         ! Model timestep in seconds
  rad_dt     =  100000.         ! Radiation timestep in seconds (set > run_time * 3600. for no radiation)
  out_dt     =      60.         ! Output timestep in seconds -- also the time interval used in storing obs for avaraging if "average_obs" is turned on!
  dz         =     250.         ! Vertical grid spacing (meters)
  ztop       =      15.         ! Model top (kilometers)
!   run_time   =       3.         ! Model run time (hours) - MCMC, ETKF
  run_time   =       1.         ! Model run time (hours) - GIG
  Tsfc       =     300.         ! Fixed model surface temperature (K) (used only in case_type = 1)
  qsfc       =      18.         ! Surface mass mixing ratio (g/kg)    (used only in case_type = 1)
  psfc       =    1000.         ! Surface pressure (hPa)              (used only in case_type = 1)
  tend_dt    =    3600.         ! Forcing time tendency interval (s)
  tend_file  = 'Frontal_Forcing_dz250m_Lat38_5-45_5_AdjQv.txt' ! Name of file containing background fields (T,P,qv,rho,heights) for the frontal case
  neg_2nd_f  =   .true.         ! Whether 2nd set of water vapor forcing fields should be negative -- note that regardless the values should be set positive below...
!  Periodic ("tend_dt") values for the specified w and qv tendency fields (note--user needs to fill the first "run_time" + 1 values)
! "Squall-line" thunderstorm case
!  Hour:       0,  1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15
  w_max_f  =  2.5, 5.,  2.5, 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.  ! Max vertical velocity (m/s)
  q_max_f  =  5.,  5.,  5.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.  ! Max qv forcing (g/kg/dt) -- Note that this is in units of per forcing time (e.g., per hour)
  z_w_f    =  2.5, 5.,  7.5, 7.5, 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.  ! Height of max vertical velocity (km)
  z_q_f    =  2.,  7.,  5.,  5.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.  ! Height of max qv forcing (km)
  half_w_f =  1.,  5.,  2.5, 2.5, 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.  ! Half-width of vertical velocity profile (km)
  half_q_f =  1.,  7.,  5.,  5.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.  ! Half-width of qv forcing profile (km)
! Forcing time: 0,    1,    2,    3,    4,    5,    6,    7,    8,    9,   10,   11,   12,   13,   14,   15,   16,   17,   18,   19
  q_max_f2  =  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0  ! Max qv forcing (g/kg/dt)    Note that this is in units of per forcing time (e.g., per hour)
  z_q_f2    =  1.5,  1.5,  1.5,  1.5,  1.5,  1.5,  1.5, 2.25,  3.0,  3.5,  4.0,  4.0,  4.0,  4.0,  4.0,  4.0,  4.0,  4.0,  4.0,  4.0  ! Height of max qv forcing (km)
  half_q_f2 =  3.0,  3.0,  3.0,  3.0,  3.0,  3.0,  3.0,  3.5,  4.0,  5.0,  5.5,  6.0,  6.0,  6.0,  6.0,  6.0,  6.0,  6.0,  6.0,  6.0  ! Half-width of qv forcing profile (km)
  start_month  =    6.         ! Month of model start time (imonth) -- for radiation
  start_day    =   22.         ! Day of the month of start time (iday) -- for radiation
  start_hour   =   12.         ! Local time (hour) of model start (hrl) -- for radiation
  latitude     =    0.         ! Latitude for radiative transfer calcs. -- for radiation
! File and output information
  write_output = .false.       ! Whether to write output to binary files for each model run
  print_less   = .false.       ! Whether to produce more diagnostic output
  print_more   = .false.       ! Whether to produce a ridiculous amount of diagnostic output
  exp_name = 'cloud_column_model_experiment'
/

! Parameter bounds settings
&RECORD1
! Parameters:      as,      bs,      ag,      bg,     tnw,     tns,     tng,     roqs,    roqg,    bnd21,   bnd1
  pmin        =   63.0,    0.10,    75.0,    0.1,    0.0001,   0.001,   0.001,    0.05,    0.05,    1.e-6,  1.e-6  ! Reasonable min from literature
  pmax        = 1000.0,     1.0,   1200.0,   0.85,     5.0,     5.0,     5.0,     1.0,     1.0,    3.e-3,  2.e-3  ! Reasonable max from literature
/

! Observations
&RECORD2
!   n_obs_times   = 6  ! Number of observation times (max is hard coded in "crm_model_error.f90") - MCMC and ETKF
!   obs_times     =  30.,  60.,  90., 120., 150., 180.   ! Observation times (in minutes since start of model run) - MCMC and ETKF
  n_obs_times   = 1  ! Number of observation times (max is hard coded in "crm_model_error.f90") = GIG precip
  obs_times     =  40.,   ! Observation times (in minutes since start of model run) - GIG precip
  average_obs     = .false. ! whether to average observations over "nobs_times" time intervals or not
  obs_start_time  =  0.0,  50.0, 100.0 ! Start times for accumulating observations 
  obs_end_time    = 50.0, 100.0, 150.0 ! End times for accumulating observations 
/

! Microphysics parameter settings
&RECORD3
! Parameters:      as,    bs,    ag,    bg,   tnw,   tns,   tng,   roqs,  roqg,  bnd21,  bnd1
! use (1) or skip (0) parameter
  ! pflag       =     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1
  ptrue       =  200.0,  0.3,  400.0,  0.4,   0.5,    0.5,   0.5,   0.2,   0.4,  1.e-3, 6.e-4 ! Values used in PV10, PB12, and PHB14
  ! ptrue       = 100.0,   0.11, 125.0,  0.15,  2.714,  0.05,  0.05,   0.1,   0.15, 2.186e-3,  1.e-5 ! Produces small precip early
!   ptrue       = 100.0,   0.11, 125.0,  0.15,  2.201,  0.05,  0.05,   0.1,   0.15, 4.859e-4,  1.e-5 ! Produces moderate (5 mm/h) precip early
!   ptrue       = 100.0,   0.11, 125.0,  0.15,  0.13095,  0.05,  0.05,   0.1,   0.15, 1.6415e-3,  1.e-5 ! Produces large (10 mm/h) precip early
!   ptrue       = 100.0,   0.11, 125.0,  0.15,  0.05,  0.05,  0.05,   0.1,   0.15, 5.e-4,  1.e-5 ! Low set
!   ptrue       = 400.0,   0.4,  500.0,  0.5,    2.0,  2.0,   2.0,    0.4,   0.5,  1.5e-3, 1.e-4 ! Middle set
!   ptrue       = 700.0,   0.7,  900.0,  0.75  , 4.0,  4.0,   4.0,    0.7,   0.9,  2.5e-3, 1.e-3 ! High set
/
