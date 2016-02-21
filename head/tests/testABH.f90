!$Id: testABH.f90 8 2014-12-01 09:13:54Z mexas $

!*robodoc*u* tests/testABH
!  NAME
!    testABH
!  SYNOPSIS

program testABH

!  PURPOSE
!    Checking: cgca_redand, part of cgca_m2red
!  DESCRIPTION
!    Checking collective AND reduction over a logical coarray.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See CGPACK_Copyright
!  USES
!    cgca testaux
!  USED BY
!    Part of CGPACK test suite
!  SOURCE

use testaux

implicit none

real, parameter :: l2 = log(real(2))
logical, parameter :: nodebug = .false.
real :: num
integer(kind=idef) :: p, nimages, img, codim(3)[*]
logical(kind=ldef) :: z[*]

!**********************************************************************73
! first executable statement

nimages=num_images()
img = this_image()

! check than n is a power of 2
p = nint(log(real(nimages))/l2)
if ( 2**p .ne. nimages) error stop "number of images is not a power of 2"

! do a check on image 1
if (img .eq. 1) then
 call getcodim(nimages,codim)
 ! print a banner
 call banner("ABH")
 write (*,'(a,i0,a)') "running on ", nimages, " images in a 3D grid"
end if

! initialise random number seed
call cgca_irs(nodebug)

! assign z
call random_number(num)

if (num .gt. 0.5) then
 z = .true.
else
 z = .false.
end if

z = .true.
if (img .eq. nimages) z = .false.

write (*,*) "image", img, "z", z

! call collective AND
call cgca_redand(z,p)

write (*,*) "image", img, "answer", z

end program testABH

!*roboend*
