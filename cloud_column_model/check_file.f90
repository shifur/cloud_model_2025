! This f90 subroutine checks for the existence of a file for input or output
! If input  is called for, and the file does not exist, it generates an error
! If output is called for, and the file exists,         it generates an error
!
! D. Posselt
! Colorado State University
! 9 May 2007
!
!-----------------------------------------------------------------------
subroutine check_file ( fname, read_file )
implicit none

! Input variables
character (*), intent(in) :: fname
logical, intent(in) :: read_file

! Local variables
logical :: file_exists

! Set flag and inquire about the file
file_exists = .false.
inquire ( file = trim(adjustl(fname)), exist = file_exists )

! For a file to be read from
if ( read_file ) then

  ! If file does not exist, generate an error
  if ( .not. file_exists ) then
    print '(a,a)','Could not find requested input file ',trim(adjustl(fname))
    stop 'Stopping'
  endif

! For a file to be written to
else

  ! If file already exists, generate an error
  if ( file_exists ) then
    print '(a,a)','Found existing output file ',trim(adjustl(fname))
    stop 'Stopping--please rename or delete this file'
  endif

endif

return
end subroutine check_file
