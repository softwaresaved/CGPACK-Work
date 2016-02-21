!$Id: testAAH.f90 8 2014-12-01 09:13:54Z mexas $

!*robodoc*u* tests/testAAH
!  NAME
!    testAAH
!  SYNOPSIS

program testAAH

!  PURPOSE
!    Checking: cgca_sld, cgca_nr
!  DESCRIPTION
!    Checking solidification with 640 nuclei,
!    coarray (100,100,100)[4,4,4]. This gives the
!    desired resolution of 1e-5. Therefore run with 64
!    images... or try other values...
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

 logical(kind=ldef),parameter :: yesdebug=.true., nodebug=.false., &
   noperiodicbc=.false.
 integer(kind=idef) :: l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3, &
   nimages,myimage,codim(3)[*]
 integer(kind=iarr),allocatable :: space1(:,:,:,:)[:,:,:]
 logical(kind=ldef) :: solid

!**********************************************************************73
! first executable statement

nimages=num_images()
myimage=this_image()

! do a check on image 1
if (myimage .eq. 1) then
 call getcodim(nimages,codim)
 ! print a banner
 call banner("AAH")
 ! print the parameter values
 call cgca_pdmp
 write (*,'(a,i0,a)') "running on ", nimages, " images in a 3D grid"
 write (*,*) "codim:", codim
end if

sync all

codim(:) = codim(:)[1]

if (myimage .eq. 2) call system("sleep 1")

l1=1
u1=100
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

! initialise coarray to liquid phase
space1(:,:,:,cgca_state_type_grain) = cgca_liquid_state

! initialise random number seed
call cgca_irs(nodebug)

sync all

! nuclei, sync all inside
call cgca_nr(space1,640,yesdebug)

call cgca_swci(space1,cgca_state_type_grain,10,"z0.raw")
sync all

! solidify 1
call cgca_sld(space1,noperiodicbc,100,10,solid)

call cgca_swci(space1,cgca_state_type_grain,10,"z100.raw")
sync all

! if solid, issue a message from image 1 and stop
if (myimage .eq. 1 .and. solid) write (*,*) "all solid, stop"
if (solid) stop

! solidify 2
call cgca_sld(space1,noperiodicbc,0,10,solid)

call cgca_swci(space1,cgca_state_type_grain,10,"z9end.raw")

call cgca_ds(space1)

end program testAAH

!*roboend*
