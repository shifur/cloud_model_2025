! This fortran 90 subroutine does the work of setting up profiles for the 
! idealized CRM error experiments. It is loosely based on the deep convection
! set up described on www.convection.info/microphysics with modifications
! to accommodate time-varying forcing.
!
! Procedure
!
! 1.  Set up vertical height vectors
! 2.  Integrate upward in height, computing first guess T, qv, and pressure
!     at layer boundaries. The procedure is as follows:
!     a.  Guess new T and P from a dry (moist) adiabatic lapse rate if the
!         the previous layer was unsaturated (saturated).
!     b.  Compute saturation vapor mixing ratio. If it is less than the 
!         surface value, set qv to qsat.
!     c.  If the level is saturated, re-compute the temperature based on the
!         moist adiabatic lapse rate.
!     d.  If temperature has been re-computed, re-compute the saturation
!         mixing ratio.
!     e.  Recompute the pressure from the layer arithmetic mean temperature
! 3.  Compute the dry air density at the layer boundaries from T and P
! 4.  Compute the temperature, pressure, and dry air density at
!     the layer mid-points via linear interpolation (linear in log-pressure).
! 5.  Integrate upward through the layer mid-points, computing the saturation
!     vapor mixing ratio from the temperature and pressure.
! 6.  Based on user-specifications, compute the vertical velocity profile
!     and the profile of vapor mixing ratio tendencies.
!
! Derek Posselt
! University of Michigan
! 3 October 2008
!
!-----------------------------------------------------------------------
subroutine ideal_setup (Tsfc, qsfc, psfc, z_lev, z_lay, T_lev, T_lay, &
                        q_lev, q_lay, p_lev, p_lay, rho_lev, rho_lay,    &
                        print_more, dz, nz, nt)
implicit none

! Input variables
integer, intent(in) :: nz, nt
real,    intent(in) :: dz
real,    intent(in) :: Tsfc, qsfc, psfc
logical, intent(in) :: print_more

! Output variables
real, dimension(nz,nt), intent(out) :: z_lev,   z_lay   & ! Heights       at layer boundaries and midpoints
                                     , T_lev,   T_lay   & ! Temperature   at layer boundaries and midpoints
                                     , q_lev,   q_lay   & ! Vapor mix rat at layer boundaries and midpoints
                                     , p_lev,   p_lay   & ! Pressure      at layer boundaries and midpoints
                                     , rho_lev, rho_lay   ! Density       at layer boundaries and midpoints

! Local variables
integer :: k, t
integer, parameter :: file_unit = 100     ! Unit number for text file containing profile data
real, parameter :: g       =    9.81    & ! Gravitational acceleration (m/s2)
                 , Rd      =  287.06    & ! Gas constant of dry air (J/K/kg)
                 , Rv      =  461.5     & ! Gas constant of water vapor (J/K/kg)
                 , cp      = 1004.5     & ! Specific heat of dry air at constant pressure (J/K/kg)
                 , Lv      = 2.501e6    & ! Latent heat of vaporization (J/kg)
                 , d_lapse =    9.77e-3 & ! Dry adiabatic lapse rate (K/m)
                 , eps     = 0.622        ! Constant in computation of saturation vapor mixing ratio

real, parameter :: z0w =  5.0e3  & ! Vertical velocity center height (meters)
                 , z0q =  7.0e3  & ! Vapor mixing ratio tendency center height (meters)
                 , Dw  =  5.0e3  & ! Vertical velocity half-width
                 , Dq  =  7.0e3    ! Vapor mixing ratio tendency half-width

real    :: e_sat, q_sat, m_lapse, t_mid, pi, half_pi

pi = 2.0 * acos(0.0)
half_pi = 0.5 * pi

! if ( print_more ) print '(a,2(1X,i5))', 'Input dimensions (nz, nt): ',nz,nt
! if ( print_more ) print*, size(z_lay)

! Compute vertical height vectors
! print*,'Z_lev(1,1): ',z_lev(1,1)
! if ( print_more ) print '(a)', 'Checkpoint 00'
do t = 1, nt
z_lev(1,t) = 0.0
! if ( print_more ) print '(a,i5)', 'Checkpoint t = ',t
z_lay(1,t) = 0.0
z_lev(2,t) = 0.0
z_lay(2,t) = 0.5 * dz
!if ( print_more ) print '(a)', 'Checkpoint 0'
do k = 3, nz
  z_lev(k,t) = z_lev(k-1,t) + dz
  z_lay(k,t) = z_lay(k-1,t) + dz
!   if (print_more) print*,'K = ',k
enddo
enddo ! End loop over times

! if ( print_more ) print '(a)', 'Checkpoint 1'

! Loop over vertical height levels, computing T, P, and qv
t_lev(1,:) = Tsfc
t_lev(2,:) = Tsfc
q_lev(1,:) = qsfc
q_lev(2,:) = qsfc
p_lev(1,:) = psfc
p_lev(2,:) = psfc
! if ( print_more ) print '(a)', 'Checkpoint 2'

do k = 3, nz

  ! Guess temperature using the dry adiabatic lapse rate
  t_lev(k,:) = t_lev(k-1,:) - (dz * d_lapse)
  t_mid      = 0.5 * (t_lev(k-1,1) + t_lev(k,1))

!   print*, 't_lev',t_lev(k)

  ! Compute pressure based on layer mean temperature and lower level pressure
  p_lev(k,:) = p_lev(k-1,:) * exp ( ( -g * dz ) / ( Rd  * t_mid ) )
!   print*, 'p_lev',p_lev(k), p_lev(k-1), -g, dz, Rd, t_mid, exp ( ( -g * dz ) / ( Rd  * t_mid ) )

  ! Compute saturation vapor pressure from T and P
  e_sat = 6.112 * exp ( ( 17.67 * (t_lev(k,1) - 273.15) ) / &
                        ( 243.5 +  t_lev(k,1) - 273.15  ) ) &
                * 100. ! Convert to Pa from hPa...
!   print*, 'e_sat',e_sat
  ! Adjust saturation vapor pressure if T < freezing
  if ( t_lev(k,1) < 273.15 ) &
    e_sat = e_sat * ( t_lev(k,1) / 273.15 )**2.66
!   print*, 'e_sat',e_sat

  ! Compute saturation vapor mixing ratio from saturation vapor pressure
  q_sat = eps * e_sat / (p_lev(k,1) - e_sat * ( 1.0 - eps ) )
!   print*, 'q_sat',q_sat

  ! If saturation vapor mixing ratio < surface value, then set qv = q_sat and recompute T and P
  if ( q_sat < qsfc ) then

    q_lev(k,:) = q_sat

    ! Compute moist adiabatic lapse rate (based on temp obtained from dry adiabatic lapse--gives a small error)
    m_lapse = d_lapse * ( 1.0 + (Lv    * q_sat) / (Rd * t_mid)         ) / &
                        ( 1.0 + (Lv**2 * q_sat) / (Rv * t_mid**2 * cp) )
    ! Re-compute temperature, layer-mid temperature, pressure, and saturation vapor mixing ratio
    t_lev(k,:) = t_lev(k-1,:) - (dz * m_lapse)
    t_mid    = 0.5 * (t_lev(k-1,1) + t_lev(k,1))
!     print*, 't_lev',t_lev(k,1)
    p_lev(k,:) = p_lev(k-1,:) * exp ( ( -g * dz ) / ( Rd  * t_mid ) )
!     print*, 'p_lev',p_lev(k,1), p_lev(k-1,1), -g, dz, Rd, t_mid, exp ( ( -g * dz ) / ( Rd  * t_mid ) )
    ! Compute saturation vapor pressure from T and P
    e_sat = 6.112 * exp ( ( 17.67 * (t_lev(k,1) - 273.15) ) / &
                          ( 243.5 +  t_lev(k,1) - 273.15  ) ) &
                  * 100. ! Convert to Pa from hPa...
!     print*, 'e_sat',e_sat
    ! Adjust saturation vapor pressure if T < freezing
    if ( t_lev(k,1) < 273.15 ) &
      e_sat = e_sat * ( t_lev(k,1) / 273.15 )**2.66
!     print*, 'e_sat',e_sat

    ! Compute saturation vapor mixing ratio from saturation vapor pressure
    q_lev(k,:) = eps * e_sat / (p_lev(k,:) - e_sat * ( 1.0 - eps ) )
!     print*, 'q_sat',q_sat

  else

    q_lev(k,:) = qsfc

  endif

enddo

! if ( print_more ) print '(a)', 'Checkpoint 3'


! Loop over vertical height levels, computing rho (use virtual temperature...)
rho_lev(1,:) = p_lev(1,:) / ( Rd * t_lev(1,:) * ( 1.0 + 0.61 * q_lev(1,:) ) )
rho_lev(2,:) = p_lev(2,:) / ( Rd * t_lev(2,:) * ( 1.0 + 0.61 * q_lev(2,:) ) )
do k = 1, nz
  rho_lev(k,:) = p_lev(k,:) / ( Rd * t_lev(k,:) * ( 1.0 + 0.61 * q_lev(k,:) ) )
enddo

! if ( print_more ) print '(a)', 'Checkpoint 4'

! Loop over layers, computing T and P at layer mid-points
p_lay(1,:) = psfc
t_lay(1,:) = tsfc
do k = 2, nz-1
  t_lay(k,:) = 0.5 * ( t_lev(k+1,:) + t_lev(k,:) )
  p_lay(k,:) = exp ( 0.5 * ( log(p_lev(k+1,:)) + log(p_lev(k,:)) ) )
enddo

! if ( print_more ) print '(a)', 'Checkpoint 5'

! Do T and P at the top via linear extrapolation
t_lay(nz,:) = 2.0 * t_lay(nz-1,:) - t_lay(nz-2,:)
p_lay(nz,:) = exp ( 2.0 * log(p_lay(nz-1,:)) - log(p_lev(nz-2,:)) )

! if ( print_more ) print '(a)', 'Checkpoint 6'

! Loop over layers, computing qv at layer mid-points
q_lay(1,:) = qsfc
do k = 2, nz-1

  ! Compute saturation vapor pressure from T and P
  e_sat = 6.112 * exp ( ( 17.67 * (t_lay(k,1) - 273.15) ) / &
                        ( 243.5 +  t_lay(k,1) - 273.15  ) ) &
                * 100. ! Convert to Pa from hPa...
!   print*, 'e_sat',e_sat
  ! Adjust saturation vapor pressure to saturation wrt ice if T < freezing
  if ( t_lay(k,1) < 273.15 ) &
    e_sat = e_sat * ( t_lay(k,1) / 273.15 )**2.66
!   print*, 'e_sat',e_sat

  ! Compute saturation vapor mixing ratio from saturation vapor pressure
  q_sat = eps * e_sat / (p_lay(k,1) - e_sat * ( 1.0 - eps ) )
!   print*, 'q_sat',q_sat

  ! If saturation vapor mixing ratio < surface value, then set qv = q_sat and recompute T and P
  if ( q_sat < qsfc ) then

    q_lay(k,:) = q_sat

    ! Compute mid-level temperature
    t_mid    = 0.5 * (t_lay(k-1,1) + t_lay(k,1))
    ! Compute saturation vapor pressure from T and P
    e_sat = 6.112 * exp ( ( 17.67 * (t_lay(k,1) - 273.15) ) / &
                          ( 243.5 +  t_lay(k,1) - 273.15  ) ) &
                  * 100. ! Convert to Pa from hPa...
!     print*, 'e_sat',e_sat
    ! Adjust saturation vapor pressure if T < freezing
    if ( t_lay(k,1) < 273.15 ) &
      e_sat = e_sat * ( t_lay(k,1) / 273.15 )**2.66
!     print*, 'e_sat',e_sat

    ! Compute saturation vapor mixing ratio from saturation vapor pressure
    q_lay(k,:) = eps * e_sat / (p_lay(k,:) - e_sat * ( 1.0 - eps ) )
!     print*, 'q_sat',q_sat

  else

    q_lay(k,:) = qsfc

  endif

enddo

! Do qv at the top via linear extrapolation
q_lay(nz,:) = 2.0 * q_lay(nz-1,:) - q_lay(nz-2,:)

! if ( print_more ) print '(a)', 'Checkpoint 7'

! Loop over layers, computing rho at layer mid-points
rho_lay(1,:) = p_lay(1,:) / ( Rd * t_lay(1,:) * ( 1.0 + 0.61 * q_lay(1,:) ) )
do k = 2, nz-1
  rho_lay(k,:) = p_lay(k,:) / ( Rd * t_lay(k,:) * ( 1.0 + 0.61 * q_lay(k,:) ) )
enddo

! Do T, P, and rho at the top via linear extrapolation
rho_lay(nz,:) = 2.0 * rho_lay(nz-1,:) - rho_lay(nz-2,:)

! if ( print_more ) print '(a)', 'Checkpoint 8'

if (print_more) then
  print '(a)','   k      z_lev     z_lay     t_lev     t_lay     p_lev     p_lay     q_lev     q_lay   rho_lev   rho_lay'
  do k = 1, nz
    print '(i5,2F10.2,2F10.3,2F10.2,2F10.7,2F10.6)', &
            k, z_lev(k,1), z_lay(k,1), t_lev(k,1), t_lay(k,1), p_lev(k,1), p_lay(k,1), &
            q_lev(k,1), q_lay(k,1), rho_lev(k,1), rho_lay(k,1)
  enddo
endif

open ( file_unit, file='Ideal_Profile_Data.txt', form='formatted', &
       action='write', status='unknown', access='sequential' )
write ( file_unit, '(a)') '   k      z_lev     z_lay     t_lev     t_lay     p_lev     p_lay     q_lev     q_lay   rho_lev   rho_lay'
do k = 1, nz
  write ( file_unit, '(i5,2F10.2,2F10.3,2F10.2,2F10.7,2F10.6)')  &
          k, z_lev(k,1), z_lay(k,1), t_lev(k,1), t_lay(k,1), p_lev(k,1), p_lay(k,1), &
          q_lev(k,1), q_lay(k,1), rho_lev(k,1), rho_lay(k,1)
enddo


return
end subroutine ideal_setup
