! This fortran 90 subroutine does the work of reading in profiles for the 
! idealized CRM error experiments. It assumes the user has already prepared 
! a file containing the following variables:
! 
! header (text)
! nforcing, nlevs
! header (text)
! layer boundary heights (meters)
! header (text)
! layer midpoint heights (meters)
! header (text)
! Layer boundary pressures (Pa), time 1
! Layer boundary pressures (Pa), time 2
! ...
! Layer boundary pressures (Pa), time nforcing
! header (text)
! Layer midpoint pressures (Pa), time 1
! Layer midpoint pressures (Pa), time 2
! ...
! Layer midpoint pressures, time nforcing
! header (text)
! Layer boundary temperatures (K), time 1
! Layer boundary temperatures (K), time 2
! ...
! Layer boundary temperatures (K), time nforcing
! header (text)
! Layer midpoint temperatures (K), time 1
! Layer midpoint temperatures (K), time 2
! ...
! Layer midpoint temperatures (K), time nforcing
! header (text)
! Layer boundary vapor mixing ratios (kg/kg), time 1
! Layer boundary vapor mixing ratios (kg/kg), time 2
! ...
! Layer boundary vapor mixing ratios (kg/kg), time nforcing
! header (text)
! Layer midpoint vapor mixing ratios (kg/kg), time 1
! Layer midpoint vapor mixing ratios (kg/kg), time 2
! ...
! Layer midpoint vapor mixing ratios (kg/kg), time nforcing
! header (text)
! Layer boundary density (kg/m3), time 1
! Layer boundary density (kg/m3), time 2
! ...
! Layer boundary density (kg/m3), time nforcing
! header (text)
! Layer midpoint density (kg/m3), time 1
! Layer midpoint density (kg/m3), time 2
! ...
! Layer midpoint density (kg/m3), time nforcing
! 
! Procedure
! 1. Check and open file
! 2. Read data
! 3. Interpolate to model timestep
! 4. Close file
!
! Derek Posselt
! Naval Research Laboratory/University of Michigan
! 12 November 2010
!
!-----------------------------------------------------------------------
subroutine read_forcing ( tend_file, nz, nt, dt, tend_dt, run_time, &
                          z_lev, z_lay, T_lev, T_lay, &
                          q_lev, q_lay, p_lev, p_lay, rho_lev, rho_lay,    &
                          print_more )
implicit none

! Input variables
character ( len = 255 ), intent(in) :: tend_file
integer, intent(in) :: nz, nt
real,    intent(in) :: dt, tend_dt, run_time
logical, intent(in) :: print_more

! Output variables
real, dimension(nz,nt), intent(out) :: z_lev,   z_lay   & ! Heights       at layer boundaries and midpoints
                                     , T_lev,   T_lay   & ! Temperature   at layer boundaries and midpoints
                                     , q_lev,   q_lay   & ! Vapor mix rat at layer boundaries and midpoints
                                     , p_lev,   p_lay   & ! Pressure      at layer boundaries and midpoints
                                     , rho_lev, rho_lay   ! Density       at layer boundaries and midpoints

! Local variables
integer :: f1, f2, k, t, nforcing, nforcing_test, nlevs, nt_test
integer, parameter :: file_unit = 100     ! Unit number for text file containing profile data
real, allocatable, dimension(:)     :: z_lev_in,   z_lay_in     ! Input Heights       at layer boundaries and midpoints
real, allocatable, dimension(:,:)   :: T_lev_in,   T_lay_in   & ! Input Temperature   at layer boundaries and midpoints
                                     , q_lev_in,   q_lay_in   & ! Input Vapor mix rat at layer boundaries and midpoints
                                     , p_lev_in,   p_lay_in   & ! Input Pressure      at layer boundaries and midpoints
                                     , rho_lev_in, rho_lay_in   ! Input Density       at layer boundaries and midpoints

! Increments on variables for time interpolation
real, allocatable, dimension(:)    :: z_lev_inc, z_lay_inc    & ! Time increments on heights
                                    , p_lev_inc, p_lay_inc    & ! Time increments on pressure
                                    , t_lev_inc, t_lay_inc    & ! Time increments on temperature
                                    , q_lev_inc, q_lay_inc    & ! Time increments on mixing ratio
                                    , rho_lev_inc, rho_lay_inc  ! Time increments on density

! Check and open file
call check_file ( trim(adjustl(tend_file)),      .true. )
open ( unit = file_unit, file = trim(adjustl(tend_file)), form = 'formatted', &
       status = 'old', action = 'read' )

if ( print_more ) print '(a,a)', 'Reading forcing fields from file ',trim(adjustl(tend_file))

! Read in dimensions: number of forcing times and number of levels
read ( file_unit, * )
read ( file_unit, * ) nforcing, nlevs

! Check the number of levels and make sure we have enough forcing times
if ( nlevs .ne. nz ) then
  print '(a)','Number of levels in forcing /= number of levels defined in include/dimensions.h'
  print '(a,2i10)','Forcing, dimensions.h: ',nlevs, nz
  stop 'Check levels'
endif

nforcing_test = int(nt / tend_dt)
if ( nforcing_test .gt. nforcing ) then
  print '(a)','Model run time exceeds forcing time period'
  print '(a,2(2x,f10.2))','Run time, forcing times ',(nt*dt)/3600., (nforcing*tend_dt)/3600.
  stop 'Check number of forcing times'
endif

! Allocate input data and increments
allocate ( z_lev_in (nlevs) )
allocate ( z_lay_in (nlevs) )
allocate ( t_lev_in (nlevs,nforcing) )
allocate ( t_lay_in (nlevs,nforcing) )
allocate ( p_lev_in (nlevs,nforcing) )
allocate ( p_lay_in (nlevs,nforcing) )
allocate ( q_lev_in (nlevs,nforcing) )
allocate ( q_lay_in (nlevs,nforcing) )
allocate ( rho_lev_in (nlevs,nforcing) )
allocate ( rho_lay_in (nlevs,nforcing) )

allocate ( z_lev_inc   (nlevs) )
allocate ( z_lay_inc   (nlevs) )
allocate ( t_lev_inc   (nlevs) )
allocate ( t_lay_inc   (nlevs) )
allocate ( p_lev_inc   (nlevs) )
allocate ( p_lay_inc   (nlevs) )
allocate ( q_lev_inc   (nlevs) )
allocate ( q_lay_inc   (nlevs) )
allocate ( rho_lev_inc (nlevs) )
allocate ( rho_lay_inc (nlevs) )

! Read in height layer boundaries and midpoints
read ( file_unit, * )
read ( file_unit, * ) z_lev_in
read ( file_unit, * )
read ( file_unit, * ) z_lay_in

! Read in layer boundary and midpoint pressures
read ( file_unit, * )
do t = 1, nforcing
  read ( file_unit, * ) p_lev_in(:,t)
enddo
read ( file_unit, * )
do t = 1, nforcing
  read ( file_unit, * ) p_lay_in(:,t)
enddo

! Read in layer boundary and midpoint temperature
read ( file_unit, * )
do t = 1, nforcing
  read ( file_unit, * ) t_lev_in(:,t)
enddo
read ( file_unit, * )
do t = 1, nforcing
  read ( file_unit, * ) t_lay_in(:,t)
enddo

! Read in layer boundary and midpoint mixing ratio
read ( file_unit, * )
do t = 1, nforcing
  read ( file_unit, * ) q_lev_in(:,t)
enddo
read ( file_unit, * )
do t = 1, nforcing
  read ( file_unit, * ) q_lay_in(:,t)
enddo

! Read in layer boundary and midpoint density
read ( file_unit, * )
do t = 1, nforcing
  read ( file_unit, * ) rho_lev_in(:,t)
enddo
read ( file_unit, * )
do t = 1, nforcing
  read ( file_unit, * ) rho_lay_in(:,t)
enddo

! do k = 1,nz
!   z_lev_inc(k) = ( z_lev_in(k,2) - z_lev_in(k,1) ) / (tend_dt / dt + 1.)
! enddo
! Set initial indices on forcing
f1 = 1
f2 = 2

! Set initial increments on all variables
z_lev_inc(:) = ( z_lev_in(:) - z_lev_in(:) ) / (tend_dt / dt)
z_lay_inc(:) = ( z_lay_in(:) - z_lay_in(:) ) / (tend_dt / dt)
t_lev_inc(:) = ( t_lev_in(:,f2) - t_lev_in(:,f1) ) / (tend_dt / dt)
t_lay_inc(:) = ( t_lay_in(:,f2) - t_lay_in(:,f1) ) / (tend_dt / dt)
p_lev_inc(:) = ( p_lev_in(:,f2) - p_lev_in(:,f1) ) / (tend_dt / dt)
p_lay_inc(:) = ( p_lay_in(:,f2) - p_lay_in(:,f1) ) / (tend_dt / dt)
q_lev_inc(:) = ( q_lev_in(:,f2) - q_lev_in(:,f1) ) / (tend_dt / dt)
q_lay_inc(:) = ( q_lay_in(:,f2) - q_lay_in(:,f1) ) / (tend_dt / dt)
rho_lev_inc(:) = ( rho_lev_in(:,f2) - rho_lev_in(:,f1) ) / (tend_dt / dt)
rho_lay_inc(:) = ( rho_lay_in(:,f2) - rho_lay_in(:,f1) ) / (tend_dt / dt)

! Fill first time with initial profiles
z_lev(:,1) = z_lev_in(:)
z_lay(:,1) = z_lay_in(:)
t_lev(:,1) = t_lev_in(:,1)
t_lay(:,1) = t_lay_in(:,1)
p_lev(:,1) = p_lev_in(:,1)
p_lay(:,1) = p_lay_in(:,1)
q_lev(:,1) = q_lev_in(:,1)
q_lay(:,1) = q_lay_in(:,1)
rho_lev(:,1) = rho_lev_in(:,1)
rho_lay(:,1) = rho_lay_in(:,1)

! Interpolate in time to model timesteps
do t = 2,nt

  ! Fill interpolated variables
  z_lev(:,t) = z_lev(:,t-1) + z_lev_inc(:)
  z_lay(:,t) = z_lay(:,t-1) + z_lay_inc(:)
  t_lev(:,t) = t_lev(:,t-1) + t_lev_inc(:)
  t_lay(:,t) = t_lay(:,t-1) + t_lay_inc(:)
  p_lev(:,t) = p_lev(:,t-1) + p_lev_inc(:)
  p_lay(:,t) = p_lay(:,t-1) + p_lay_inc(:)
  q_lev(:,t) = q_lev(:,t-1) + q_lev_inc(:)
  q_lay(:,t) = q_lay(:,t-1) + q_lay_inc(:)
  rho_lev(:,t) = rho_lev(:,t-1) + rho_lev_inc(:)
  rho_lay(:,t) = rho_lay(:,t-1) + rho_lay_inc(:)

  ! If necessary, update forcing increments
  if ( (real(t)*dt) / tend_dt .eq. real( (t*int(dt)) / int(tend_dt) ) &
      .and. t .gt. 2 .and. t .lt. nt-2 ) then

    if ( print_more ) &
      print '(a,3i10)','Updating increment at timestep ',t, f1, f2

    ! Update indices on forcing
    f1 = f1 + 1
    f2 = f2 + 1

    ! do k = 1,nz
    !   z_lev_inc(k) = ( z_lev_in(k,2) - z_lev_in(k,1) ) / (tend_dt / dt + 1.)
    ! enddo
    ! Compute new forcing increments
    z_lev_inc(:) = ( z_lev_in(:) - z_lev_in(:) ) / (tend_dt / dt)
    z_lay_inc(:) = ( z_lay_in(:) - z_lay_in(:) ) / (tend_dt / dt)
    t_lev_inc(:) = ( t_lev_in(:,f2) - t_lev_in(:,f1) ) / (tend_dt / dt)
    t_lay_inc(:) = ( t_lay_in(:,f2) - t_lay_in(:,f1) ) / (tend_dt / dt)
    p_lev_inc(:) = ( p_lev_in(:,f2) - p_lev_in(:,f1) ) / (tend_dt / dt)
    p_lay_inc(:) = ( p_lay_in(:,f2) - p_lay_in(:,f1) ) / (tend_dt / dt)
    q_lev_inc(:) = ( q_lev_in(:,f2) - q_lev_in(:,f1) ) / (tend_dt / dt)
    q_lay_inc(:) = ( q_lay_in(:,f2) - q_lay_in(:,f1) ) / (tend_dt / dt)
    rho_lev_inc(:) = ( rho_lev_in(:,f2) - rho_lev_in(:,f1) ) / (tend_dt / dt)
    rho_lay_inc(:) = ( rho_lay_in(:,f2) - rho_lay_in(:,f1) ) / (tend_dt / dt)

  endif


enddo

! If called for, write out forcing fields
if ( print_more ) then

  call check_file ( 'z_lev_forcing.txt',   .false. )
  call check_file ( 'z_lay_forcing.txt',   .false. )
  call check_file ( 't_lev_forcing.txt',   .false. )
  call check_file ( 't_lay_forcing.txt',   .false. )
  call check_file ( 'p_lev_forcing.txt',   .false. )
  call check_file ( 'p_lay_forcing.txt',   .false. )
  call check_file ( 'q_lev_forcing.txt',   .false. )
  call check_file ( 'q_lay_forcing.txt',   .false. )
  call check_file ( 'rho_lev_forcing.txt', .false. )
  call check_file ( 'rho_lay_forcing.txt', .false. )
  open ( unit = file_unit+1,  file = 'z_lev_forcing.txt', form = 'formatted', status = 'new', action = 'write' )
  open ( unit = file_unit+2,  file = 'z_lay_forcing.txt', form = 'formatted', status = 'new', action = 'write' )
  open ( unit = file_unit+3,  file = 't_lev_forcing.txt', form = 'formatted', status = 'new', action = 'write' )
  open ( unit = file_unit+4,  file = 't_lay_forcing.txt', form = 'formatted', status = 'new', action = 'write' )
  open ( unit = file_unit+5,  file = 'p_lev_forcing.txt', form = 'formatted', status = 'new', action = 'write' )
  open ( unit = file_unit+6,  file = 'p_lay_forcing.txt', form = 'formatted', status = 'new', action = 'write' )
  open ( unit = file_unit+7,  file = 'q_lev_forcing.txt', form = 'formatted', status = 'new', action = 'write' )
  open ( unit = file_unit+8,  file = 'q_lay_forcing.txt', form = 'formatted', status = 'new', action = 'write' )
  open ( unit = file_unit+9,  file = 'rho_lev_forcing.txt', form = 'formatted', status = 'new', action = 'write' )
  open ( unit = file_unit+10, file = 'rho_lay_forcing.txt', form = 'formatted', status = 'new', action = 'write' )
  do t = 1,nt
    do k = 1,nz
      write (file_unit+1, '(2X,F10.2)',advance='no') z_lev(k,t)
      write (file_unit+2, '(2X,F10.2)',advance='no') z_lay(k,t)
      write (file_unit+3, '(2X,F10.2)',advance='no') t_lev(k,t)
      write (file_unit+4, '(2X,F10.2)',advance='no') t_lay(k,t)
      write (file_unit+5, '(2X,F10.2)',advance='no') p_lev(k,t)
      write (file_unit+6, '(2X,F10.2)',advance='no') p_lay(k,t)
      write (file_unit+7, '(2X,F12.8)',advance='no') q_lev(k,t)
      write (file_unit+8, '(2X,F12.8)',advance='no') q_lay(k,t)
      write (file_unit+9, '(2X,F10.5)',advance='no') rho_lev(k,t)
      write (file_unit+10,'(2X,F10.5)',advance='no') rho_lay(k,t)
    enddo
    write (file_unit+1, *)
    write (file_unit+2, *)
    write (file_unit+3, *)
    write (file_unit+4, *)
    write (file_unit+5, *)
    write (file_unit+6, *)
    write (file_unit+7, *)
    write (file_unit+8, *)
    write (file_unit+9, *)
    write (file_unit+10,*)
  enddo
  close ( file_unit+1  )
  close ( file_unit+2  )
  close ( file_unit+3  )
  close ( file_unit+4  )
  close ( file_unit+5  )
  close ( file_unit+6  )
  close ( file_unit+7  )
  close ( file_unit+8  )
  close ( file_unit+9  )
  close ( file_unit+10 )

endif

return
end subroutine read_forcing
