!$Id: cgca_m2glm.f90 99 2015-06-10 13:47:07Z mexas $

!*robodoc*m* CGPACK/cgca_m2glm
!  NAME
!    cgca_m2glm
!  SYNOPSIS

module cgca_m2glm

!  DESCRIPTION
!    Module dealing with Global to Local Mapping (glm) and vice versa
!  AUTHOR 
!    Anton Shterenlikht
!  COPYRIGHT
!    See CGPACK_Copyright
!  CONTAINS
!    cgca_gl, cgca_lg
!  USES
!    cgca_m1co
!  USED BY
!    cgca_m3nucl
!  SOURCE

use cgca_m1co
implicit none

private
public :: cgca_gl, cgca_lg

contains

!*roboend*


!*robodoc*s* cgca_m2glm/cgca_gl
!  NAME
!    cgca_gl
!  SYNOPSIS

subroutine cgca_gl(super,coarray,imgpos,local)

!  INPUTS

integer(kind=idef),intent(in) :: super(3)
integer(kind=iarr),allocatable,intent(inout) :: coarray(:,:,:,:)[:,:,:]

! OUTPUT

integer(kind=idef),intent(out) :: imgpos(3),local(3)

! DESCRIPTION
!   This routine converts a cell coordinate from a global, super, array
!   to the image coordinates in the coarray grid and the local cell
!   coordinates in this image :
!   - super(3) are cell coordinates in a super array
!   - coarray is the model
!   - imgpos(3) is the image position in the grid
!   - local(3) are cell coordinates in that image's array
! NOTES
!   The global coordinates must start from 1!
!
!   Any image can call this routine
! USES
! USED BY
! SOURCE

integer :: & 
 lbr(4)   ,& ! lower bounds of the "real" coarray, lbv+1
 ubr(4)   ,& ! upper bounds of the "real" coarray, ubv-1
 szr(3)   ,& ! size or the "real" coarray, ubr-lbr+1
 lcob(3)  ,& ! lower cobounds of the coarray
 ucob(3)  ,& ! upper cobounds of the coarray
 usup(3)  ,& ! upper bound of the super array, szr*(ucob-lcob+1)
 thisimage

thisimage = this_image()

! check for allocated

if (.not. allocated(coarray)) then
 write (*,'(a,i0)') "ERROR: cgca_gl: image", thisimage
 write (*,'(a)') "ERROR: cgca_gl: coarray is not allocated"
 error stop
end if

lbr=lbound(coarray)+1
ubr=ubound(coarray)-1

! the 4th dimension is to do with the number of cell state
! types. This is not relevant here.
szr=ubr(1:3)-lbr(1:3)+1

lcob=lcobound(coarray)
ucob=ucobound(coarray)
usup=szr*(ucob-lcob+1)

! check for bounds

if (any(super .gt. usup) .or. any(super .lt. 1)) then
 write (*,'(a,i0)') "ERROR: cgca_gl: image", thisimage
 write (*,'(a)') "ERROR: cgca_gl: one or more super array&
                  & coordinate(s) are ouside the bounds"
 write (*,'(a,3(i0,tr1))') "ERROR: cgca_gl: super array coord: ",super
 write (*,'(a)') "ERROR: cgca_gl: lower bound must be 1"
 write (*,'(a,3(i0,tr1))') "ERROR: cgca_gl: upper bounds: ", usup
 error stop
end if

! actual calculation

imgpos = lcob + (super-1)/szr
local = lbr(1:3) + super-szr*(imgpos-lcob) - 1

! checks after

if (any(imgpos .gt. ucob) .or. any(imgpos .lt. lcob)) then
 write (*,'(a,i0)') "ERROR: cgca_gl: image", thisimage
 write (*,'(a)') "ERROR: cgca_lg: one or more image positions&
                  & are ouside the bounds"
 write (*,'(a,3(i0,tr1))') &
  "ERROR: cgca_gl: image positions: ",imgpos
 write (*,'(a,3(i0,tr1))') &
  "ERROR: cgca_gl: lower image grid bounds: ", lcob
 write (*,'(a,3(i0,tr1))') &
  "ERROR: cgca_gl: upper image grid bounds: ", ucob
 error stop
end if

if (any(local .gt. ubr(1:3)) .or. any(local .lt. lbr(1:3))) then
 write (*,'(a,i0)') "ERROR: cgca_gl: image", thisimage
 write (*,'(a)') "ERROR: cgca_lg: one or more local coordinates &
                  & are ouside the bounds"
 write (*,'(a,3(i0,tr1))') "ERROR: cgca_gl: local coordinates: ",local
 write (*,'(a,3(i0,tr1))') "ERROR: cgca_gl: lower bounds: ", lbr
 write (*,'(a,3(i0,tr1))') "ERROR: cgca_gl: upper bounds: ", ubr
 error stop
end if

end subroutine cgca_gl

!*roboend*


!*robodoc*s* cgca_m2glm/cgca_lg
!  NAME
!    cgca_lg
!  SYNOPSIS

subroutine cgca_lg(imgpos,local,coarray,super)

!  INPUTS

integer(kind=idef),intent(in) :: imgpos(3),local(3)
integer(kind=iarr),allocatable,intent(inout) :: coarray(:,:,:,:)[:,:,:]

! OUTPUT

integer(kind=idef),intent(out) :: super(3)

! DESCRIPTION
!  This routine converts the image coordinates in the grid and the local
!  cell coordinates in this image into the global cell coordinates in
!  the super array:
!   - imgpos(3) is the image position in the grid
!   - local(3) are cell coordinates in that image's array
!   - coarray is the model
!   - super(3) are cell coordinates in a super array
! NOTES
!   The global, super, coordinates must start from 1!
!
!   Any image can call this routine
! USES
! USED BY
!   cgca_gbf1f
! SOURCE

integer :: & 
 lbr(4)   ,& ! lower bounds of the "real" coarray, lbv+1
 ubr(4)   ,& ! upper bounds of the "real" coarray, ubv-1
 szr(3)   ,& ! size or the "real" coarray, ubr-lbr+1
 lcob(3)  ,& ! lower cobounds of the coarray
 ucob(3)  ,& ! upper cobounds of the coarray
 usup(3)  ,& ! upper bound of the super array, szr*(ucob-lcob+1)
 thisimage

thisimage = this_image()

! check for allocated

if (.not. allocated(coarray)) then
 write (*,'(a,i0)') "ERROR: cgca_lg: image", thisimage
 write (*,'(a)') "ERROR: cgca_lg: coarray is not allocated"
 error stop
end if

lbr=lbound(coarray)+1
ubr=ubound(coarray)-1

! the 4th dimension is to do with the number of cell state
! types. This is not relevant here.
szr=ubr(1:3)-lbr(1:3)+1

lcob=lcobound(coarray)
ucob=ucobound(coarray)
usup=szr*(ucob-lcob+1)

! check for bounds

if (any(imgpos .gt. ucob) .or. any(imgpos .lt. lcob)) then
 write (*,'(a,i0)') "ERROR: cgca_lg: image", thisimage
 write (*,'(a)') "ERROR: cgca_lg: one or more image positions&
                  & are ouside the bounds"
 write (*,'(a,3(i0,tr1))') "ERROR: cgca_lg: image positions: ",imgpos
 write (*,'(a,3(i0,tr1))') &
  "ERROR: cgca_lg: lower image grid bounds: ", lcob
 write (*,'(a,3(i0,tr1))') &
  "ERROR: cgca_lg: upper image grid bounds: ", ucob
 error stop
end if

if (any(local .gt. ubr(1:3)) .or. any(local .lt. lbr(1:3))) then
 write (*,'(a,i0)') "ERROR: cgca_lg: image", thisimage
 write (*,'(a)') "ERROR: cgca_lg: one or more local coordinates&
                  & are ouside the bounds"
 write (*,'(a,3(i0,tr1))') &
  "ERROR: cgca_lg: local coordinates: ", local
 write (*,'(a,3(i0,tr1))') &
  "ERROR: cgca_lg: lower bounds: ", lbr
 write (*,'(a,3(i0,tr1))') &
  "ERROR: cgca_lg: upper bounds: ", ubr
 error stop
end if

! actual calculation

super = szr*(imgpos-lcob) + local-lbr(1:3)+1

! check for bounds

if (any(super .gt. usup) .or. any(super .lt. 1)) then
 write (*,'(a,i0)') "ERROR: cgca_lg: image", thisimage
 write (*,'(a)') "ERROR: cgca_lg: one or more super array &
                  & coordinates are ouside the bounds"
 write (*,'(a,3(i0,tr1))') &
  "ERROR: cgca_lg: super array coord: ",super
 write (*,'(a)') "ERROR: cgca_lg: lower bound must be 1"
 write (*,'(a,3(i0,tr1))') "ERROR: cgca_lg: upper bounds: ", usup
 error stop
end if

end subroutine cgca_lg

!*roboend*

end module cgca_m2glm
