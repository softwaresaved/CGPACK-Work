!$Id: cgca_m2out.f90 228 2016-02-19 14:23:28Z mexas $

!*robodoc*m* CGPACK/cgca_m2out
!  NAME
!    cgca_m2out
!  SYNOPSIS

module cgca_m2out

!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  DESCRIPTION
!    Module dealing with output 
!  CONTAINS
!    Subroutines: cgca_swci, ==> cgca_csvi <== obsolete remove later!
!    Submodules: m2out_sm1, m2out_sm2_mpi.
!  USES
!    cgca_m1co
!  USED BY
!    cgca
!  SOURCE

use cgca_m1co
implicit none

private
public :: cgca_swci, cgca_pc, cgca_pswci

interface
  module subroutine cgca_pc( coarray, stype, fname )
    ! coarray - what array to dump
    ! stype - what cell state type to dump
    ! fname - what file name to use
    integer( kind=iarr ), allocatable, intent( in ) ::                 &
      coarray(:,:,:,:)[:,:,:]
    integer( kind=idef ),intent( in ) :: stype
    character( len=* ), intent( in ) :: fname
  end subroutine cgca_pc

  module subroutine cgca_pswci( coarray, stype, fname )
    ! Parallel Stream Write Coarray of Integers:
    ! - coarray - what array to dump
    ! - stype - what cell state type to dump
    ! - fname - what file name to use
    integer( kind=iarr ), allocatable, intent( in ) ::                 &
      coarray(:,:,:,:)[:,:,:]
    integer( kind=idef ),intent( in ) :: stype
    character( len=* ), intent( in ) :: fname
  end subroutine cgca_pswci

end interface

contains

!*roboend*


!*robodoc*s* cgca_m2out/cgca_swci
!  NAME
!    cgca_swci
!  SYNOPSIS

subroutine cgca_swci( coarray, stype, iounit, fname )

!  INPUTS
 
integer( kind=iarr ),allocatable,intent( in ) :: coarray(:,:,:,:)[:,:,:]
integer( kind=idef ),intent( in ) :: stype, iounit
character( len=* ),intent( in ) :: fname

!  SIDE EFFECTS
!    A single binary file is created on image 1 with contents of coarray.
!  DESCRIPTION
!    Stream Write Coarray of Integers:
!    - coarray - what array to dump
!    - stype - what cell state type to dump
!    - iounit - which I/O unit to use
!    - fname - what file name to use
!  NOTES
!    All images call this routine!
!    However only image 1 does all the work.
!    The other images are waiting. 
!  USES
!    none
!  USED BY
!    none, end user.
!  SOURCE

integer :: errstat, coi1, coi2, coi3, i2, i3, &
  lb(4),   & ! lower bounds   of the coarray
  ub(4),   & ! upper bounds   of the coarray
  lcob(3), & ! lower cobounds of the coarray
  ucob(3)    ! upper cobounds of the coarray

! Only image1 does this. All other images do nothing.
! So sync all probably should be used after a call to
! this routine in the program.

main: if ( this_image() .eq. 1 ) then
 errstat = 0

 ! Assume the coarray has halos. Don't write those.
 lb   = lbound( coarray ) + 1
 ub   = ubound( coarray ) - 1
 lcob = lcobound( coarray )
 ucob = ucobound( coarray )

!write (*,*) "DEBUG: cgca_swci: lb: " , lb
!write (*,*) "DEBUG: cgca_swci: ub: " , ub
!write (*,*) "DEBUG: cgca_swci: lcob: " , lcob
!write (*,*) "DEBUG: cgca_swci: ucob: " , ucob

 open( unit=iounit, file=fname, form="unformatted", access="stream", &
       status="replace", iostat=errstat )
 if ( errstat .ne. 0 ) then
  write (*,'(a)') "ERROR: cgca_swci: cannot open file for writing"
  write (*,'(a,i0)') "ERROR: cgca_swci: error code: ", errstat
  error stop
 end if

!write (*,*) "DEBUG: cgca_swci: starting data output"

 ! nested loops for writing in correct order from all images
 do coi3 = lcob(3), ucob(3)
   do i3 = lb(3), ub(3)
     do coi2 = lcob(2), ucob(2)
       do i2 = lb(2), ub(2)
         do coi1 = lcob(1), ucob(1)

  write(unit=iounit, iostat=errstat) &
    coarray(lb(1):ub(1),i2,i3,stype)[coi1,coi2,coi3]
  if (errstat .ne. 0) then
    write (*,'(a)') "ERROR: cgca_swci: cannot write to file"
    write (*,'(a,i0)') "ERROR: cgca_swci: error code: ", errstat
    error stop
  end if

!write (*,*) "DEBUG: cgca_swci: wrote cells with: i2, i3, coi1, coi2, coi3", i2, i3, coi1, coi2, coi3

         end do
       end do
     end do
   end do
 end do

!write (*,*) "DEBUG: cgca_swci: finished data output"

 close( unit=iounit, iostat=errstat )
 if ( errstat .ne. 0 ) then
  write (*,'(a)') "ERROR: cgca_swci: cannot close file"
  write (*,'(a,i0)') "ERROR: cgca_swci: error code: ", errstat
  error stop
 end if

end if main

end subroutine cgca_swci

!*roboend*


!*robodoc*s* cgca_m2out/cgca_csvi
!  NAME
!    cgca_csvi
!  SYNOPSIS
 
subroutine cgca_csvi(coarray,iounit,fname)

!  INPUTS

integer(kind=iarr),allocatable,intent(in) :: coarray(:,:,:)[:,:,:]
integer(kind=idef),intent(in) :: iounit
character(len=*),intent(in) :: fname

!  SIDE EFFECTS
!    creates a text file from image 1 and writes coarray in it
!  DESCRIPTION
!    Comma Separated Values of Integers.
!    All images call this routine.
!    However only image 1 does all the work.
!    The other images are waiting. 
!  WARNINGS
!    DO NOT USE!!!!
!    This subroutine is not very useful, and therefore
!    the work on it stopped.
!    In fact it might not work at all.
!    So for now it is just a dummy.
!    If you want to use it, do so on your own risk.
!    In fact it is not even accessible from outside of the enclosing module.
!    It might be removed in future version with no notice!
!  USES
!    none
!  USED BY
!    none, DO NOT USE!!!
!  SOURCE

integer :: errstat, coi1, coi2, coi3, i1, i2, i3, &
  lbv(3)      , & ! lower bounds of the "virtual" coarray
  ubv(3)      , & ! upper bounds of the "virtual" coarray
  lbr(3)      , & ! lower bounds of the "real" coarray, lbv+1
  ubr(3)      , & ! upper bounds of the "real" coarray, ubv-1
  lcob(3)     , & ! lower cobounds of the coarray
  ucob(3)     , & ! upper cobounds of the coarray
  dimr(3)         ! shape of the "real" coarray

! leave immediately, 
write (*,*) "PROBLEM: cgca_csvi is not ready for use, skipping"
return

! Only image1 does this. All other images do nothing.
! So sync all probably should be used after a call to
! this routine in the program.

main: if (this_image() .eq. 1) then
  errstat = 0

  ! Assume the coarray has halos. Don't write those.
  lbv=lbound(coarray)
  ubv=ubound(coarray)
  lbr=lbv+1
  ubr=ubv-1
  dimr=ubr-lbr+1
  lcob=lcobound(coarray)
  ucob=ucobound(coarray)

  open(unit=iounit,file=fname, form="formatted", status="replace", &
    iostat=errstat)
  if (errstat.ne.0) then
    write (*,*) &
      "ERROR: cgca_csvi: cannot open file for writing"
    error stop
  end if

  ! write in any order
  do concurrent ( i1 =  lbr(1): ubr(1), &
                  i2 =  lbr(2): ubr(2), &
                  i3 =  lbr(3): ubr(3), &
                coi1 = lcob(1):ucob(1), &
                coi2 = lcob(2):ucob(2), &
                coi3 = lcob(3):ucob(3) )
    ! need to calculate the global cell coordinates
    write (unit=iounit, iostat=errstat, fmt='(3(i0,a),i0)' ) &
      i1 + dimr(1)*(coi1-1), " , ", & 
      i2 + dimr(2)*(coi2-1), " , ", & 
      i3 + dimr(3)*(coi3-1), " , ", & 
      coarray(i1,i2,i3)[coi1,coi2,coi3]
    if (errstat.ne.0) then
      write (*,'(a)') "ERROR: cgca_csvi: cannot write to file"
      write (*,'(a,i0)') "ERROR: cgca_csvi: error code: ", errstat
      ! ideally should issue "error stop" here, but it is an image
      ! control statement, which are not allowed inside do concurrent 
    end if
  end do

  close(unit=iounit,iostat=errstat)
  if (errstat.ne.0) error stop "ERROR: cgca_csvi: cannot close file"

end if main

end subroutine cgca_csvi

!*roboend*

end module cgca_m2out
