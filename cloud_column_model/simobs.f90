!              SUBROUTINE TO ACCUMULATE OBS TOTALS FOR MCMC
!
! This routine is written in F90 to conform to other portions of the
! GCE model. It computes the "observations" used in the MCMC procedure.
! The list of observations is placed at the end of the driver program
! and is kept in array "obs_gce".
!
subroutine simobs(nz, twcz, pi, dpt, ta, dqv, qa, &
                  qcl, qrn, qci, qcs, qcg, pwv, lwp, iwp, ctop)
implicit none

! Input variables
integer, intent(in) :: nz
real,    intent(in) :: pi(nz), dpt(1,1,nz), ta(nz), dqv(1,1,nz), qa(nz), &
                       qcl(1,1,nz), qrn(1,1,nz), qci(1,1,nz), qcs(1,1,nz), qcg(1,1,nz)
real,    intent(in) :: twcz   ! cloud threshold (g/g) for radiation calcs...

! Output variables (the "observations")
real, intent(out)   :: pwv, lwp, iwp, ctop

! Local variables
integer :: k
real    :: ptop, pbot
real    :: cpor
logical found_top

! Constants
real, parameter :: grav=9.81
real, parameter :: rgas=287.
real, parameter :: cp=1004.

! Missing value setting
real, parameter :: missing = -999.

! Get constant for pressure exponentiation
  cpor = cp / rgas

! Reset counters, logical flags, and observations
  pwv = 0.
  lwp = 0.
  iwp = 0.
  ctop= missing

found_top = .false.

! Set end point in GCE grid

! Extrapolate top pressure
ptop = 0.
ptop = 100000. * ( 0.5 * ( 3.*pi(nz-1) - pi(nz-2) ) )**cpor   !Compute pressure


! Search downward through model layers computing pwv, lwp, iwp
do k=nz,2,-1

  ! Compute layer bottom pressure
  ! interp/extrap pi first, then convert to pressure
  pbot = 0.
  if (k.gt.2) then !Interp pi, then compute pressure
    pbot = 100000. * ( 0.5 * (   pi(k) + pi(k-1) ) )**cpor
  else 
    pbot = 100000. * ( 0.5 * ( 3*pi(2) - pi(3)   ) )**cpor !extrapolate pi to bottom layer
  endif

  ! Compute PWV, LWP, IWP
  pwv = pwv+((pbot-ptop)/grav)*(qa(k)+dqv(1,1,k))
  lwp = lwp+((pbot-ptop)/grav)*(qcl(1,1,k)+qrn(1,1,k))
  iwp = iwp+((pbot-ptop)/grav)*(qci(1,1,k)+qcs(1,1,k)+qcg(1,1,k))

  ! Look for top of ANY cloud for cloud top temp calculation
  if ( ((qcl(1,1,k)+qrn(1,1,k)+qci(1,1,k)+qcs(1,1,k)+qcg(1,1,k)) .ge. twcz) .and. &
        (.not. found_top) ) then
    found_top = .true.
    ctop = pi(k)*(dpt(1,1,k)+ta(k)) ! T=pi*theta
  endif

  ! Set top pressure equal to bottom pressure
  ptop = pbot

! End k-loop
enddo

! Return to calling program
return
end subroutine simobs

