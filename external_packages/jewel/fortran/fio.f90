! JEWEL 2.6.0 C++ port -- minimal Fortran I/O shim.
! PYTHIA 6 writes its messages to Fortran unit MSTU(11) and PYRGET/PYRSET
! read/write the random-number state on unformatted Fortran units. These
! routines let the C++ code open, write to and close such units so that the
! JEWEL log file receives both JEWEL and PYTHIA output in the original order.
subroutine jewel_fopen(unit, path, pathlen, unformatted, status_old, ierr)
  implicit none
  integer, intent(in) :: unit, pathlen, unformatted, status_old
  character(len=pathlen), intent(in) :: path
  integer, intent(out) :: ierr
  character(len=7) :: st
  if (status_old /= 0) then
    st = 'old'
  else
    st = 'unknown'
  end if
  if (unformatted /= 0) then
    open(unit=unit, file=path, access='sequential', form='unformatted', status=trim(st), iostat=ierr)
  else
    open(unit=unit, file=path, status=trim(st), iostat=ierr)
  end if
end subroutine jewel_fopen

subroutine jewel_fclose(unit)
  implicit none
  integer, intent(in) :: unit
  close(unit, status='keep')
end subroutine jewel_fclose

subroutine jewel_fwrite(unit, text, len)
  implicit none
  integer, intent(in) :: unit, len
  character(len=len), intent(in) :: text
  write(unit, '(A)') text
end subroutine jewel_fwrite

subroutine jewel_fflush(unit)
  implicit none
  integer, intent(in) :: unit
  flush(unit)
end subroutine jewel_fflush
