!$Id: testAAE.f90 8 2014-12-01 09:13:54Z mexas $

!*robodoc*u* tests/testAAE
!  NAME
!    testAAE
!  SYNOPSIS

program testAAE

!  PURPOSE
!    Checking: cgca_irs, cgca_nr
!  DESCRIPTION
!    Checking nucleation and seed initialisation.
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

 integer(kind=idef) :: l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3, &
   nimages,myimage,codim(3)[*]
 logical(kind=ldef),parameter :: yesdebug=.true., nodebug=.false.
 integer(kind=iarr),allocatable :: space1(:,:,:,:)[:,:,:]

!**********************************************************************73
! first executable statement

nimages=num_images()
myimage=this_image()

! do a check on image 1
if (myimage .eq. 1) then
 call getcodim(nimages,codim)
 ! print a banner
 call banner("AAE")
 ! print the parameter values
 call cgca_pdmp
 write (*,'(a,i0,a)') "running on ", nimages, " images in a 3D grid"
 write (*,*) "codim:", codim
end if

sync all

codim(:) = codim(:)[1]

if (myimage .eq. 2) call system("sleep 1")

l1=1
u1=10
l2=l1
u2=u1
l3=l1
u3=u1
col1=1
cou1=codim(1)
col2=1
cou2=codim(2)
col3=1

if (myimage .eq. 1) then
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "bounds: (",l1,u1,l2,u2,l3,u3
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "cobounds: (",col1,cou1,col2,cou2,col3, &
     nimages/((cou1-col1+1)*(cou2-col2+1))+col3-1
end if

! allocate coarray
call cgca_as(l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3,1,space1)

! initialise random number seed
call cgca_irs(yesdebug)

! initialise coarray to liquid state
space1 = cgca_liquid_state
sync all 

call cgca_swci(space1,cgca_state_type_grain,10,"z1.raw")
sync all

! nucleate, sync in the routine, no need to sync in the program
call cgca_nr(space1,10,yesdebug)

call cgca_swci(space1,cgca_state_type_grain,10,"z2.raw")
sync all

call cgca_ds(space1)

end program testAAE

!*roboend*
