!$Id: testABC.f90 172 2015-11-04 11:59:43Z mexas $

!*robodoc*u* tests/testABC
!  NAME
!    testABC
!  SYNOPSIS

program testABC

!  PURPOSE
!    Checking: cgca_pdmp
!  DESCRIPTION
!    Dump the global CGPACK parameters to stdout.
!    Any or all images can call this routine.
!    Here only image 1 does it.
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

! only image 1 works

if ( this_image() .eq. 1) then
 call banner("ABC") ! print a banner
 call cgca_pdmp     ! print parameter values
end if

end program testABC

!*roboend*
