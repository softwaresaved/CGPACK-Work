!$Id: cgca_m2rnd.f90 14 2014-12-01 10:14:16Z mexas $

!*robodoc*m* CGPACK/cgca_m2rnd
!  NAME
!    cgca_m2rnd
!  SYNOPSIS

module cgca_m2rnd

!  DESCRIPTION
!    Module dealing with random number generation
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See CGPACK_Copyright
!  CONTAINS
!    cgca_irs
!  USES
!    cgca_m1co
!  USED BY
!    cgca
!  SOURCE

use cgca_m1co
implicit none

private
public :: cgca_irs

contains

!*roboend*

!*robodoc*s* cgca_m2rnd/cgca_irs
!  NAME
!    cgca_irs
!  SYNOPSIS

subroutine cgca_irs(debug)

!  INPUT

logical(kind=ldef),intent(in) :: debug

!  SIDE EFFECTS
!    initialise random seed on all images   
!  DESCRIPTION
!    Initialise random seed based on system_clock intrinsic,
!    adapted from:
!    http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html.
!    Note that the seed is based on THIS_IMAGE intrinsic, hence
!    each image uses a different seed.
!  USES
!    none
!  USED BY
!    none, end user
!  SOURCE

integer :: i,n,clock,errstat=0
integer,allocatable :: seed(:)
         
call random_seed(size = n)

allocate( seed(n), stat=errstat )
if ( errstat .ne. 0) stop "ERROR: cgca_irs: cannot allocate seed"

call system_clock(count=clock)
          
seed = int(real(clock)/real(this_image())) +  &
  999999937*(/ (i - 1, i = 1, n) /)

if (debug) write (*,*) "image:",this_image(), "; size:",n,"; seed",seed

call random_seed(put = seed)
          
deallocate( seed, stat=errstat )
if ( errstat .ne. 0 ) stop "ERROR: cgca_irs: cannot deallocate seed"

end subroutine cgca_irs

!*roboend*

end module cgca_m2rnd
