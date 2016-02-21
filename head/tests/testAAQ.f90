!$Id: testAAQ.f90 8 2014-12-01 09:13:54Z mexas $

!*robodoc*u* tests/testAAQ
!  NAME
!    testAAQ
!  SYNOPSIS

program testAAQ

!  PURPOSE
!    Checking: cgca_rt, cgca_ckrt
!  DESCRIPTION
!    Checking rotation tensor calculation and orthogonality checking
!    routine.
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
   periodicbc=.true.
 real,parameter :: gigabyte=real(2**30), resolution=1.0e-5

 real(kind=rdef),allocatable :: grt(:,:,:)[:,:,:]

 integer(kind=idef) :: l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3, &
   cou3,   &
   nuc,    & ! number of nuclei in the model
   nimages,flag,codim(3)[*],i
 integer(kind=iarr),allocatable :: space(:,:,:,:)[:,:,:]
 integer(kind=ilrg),allocatable :: grainvol(:)[:,:,:]
 integer(kind=ilrg) :: icells,mcells

 logical(kind=ldef) :: solid,image1

 real :: image_storage

!**********************************************************************73
! first executable statement

nimages=num_images()
image1=.false.
if (this_image().eq.1) image1=.true.

! do a check on image 1
if (image1) then
 call getcodim(nimages,codim)
 ! print a banner
 call banner("AAQ")
 ! print the parameter values
 call cgca_pdmp
 write (*,'(a,i0,a)') "running on ", nimages, " images in a 3D grid"
 write (*,*) "codim:", codim
end if

sync all

codim(:) = codim(:)[1]

if (image1) call system("sleep 1")

l1=1
l2=l1
l3=l1

! The array size is only controlled by this value
! in this program.
u1=100
u2=u1
u3=u1

col1=1
cou1=codim(1)-col1+1
col2=1
cou2=codim(2)-col2+1
col3=1
cou3=codim(3)-col3+1

! total number of cells in a coarray
icells = int(u1-l1+1,kind=ilrg) * int(u2-l2+1,kind=ilrg) * &
  int(u3-l3+1,kind=ilrg)

! total number of cells in the model
mcells = icells * int(codim(1),kind=ilrg) * int(codim(2),kind=ilrg) * &
  int(codim(3),kind=ilrg)

! total number of nuclei
nuc = int( resolution * mcells)

if (image1) then
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "bounds: (",l1,u1,l2,u2,l3,u3
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "cobounds: (",col1,cou1,col2,cou2,col3,cou3

  ! An absolute minimum of storage, in GB, per image.
  image_storage = real( &
    ! 2nd space array is allocated in _sld
    2 * icells*storage_size(space) &
    ! grain volumes
    + nuc*storage_size(grainvol) &
    ! rotation tensors
    + nuc*9*storage_size(grt))/8/gigabyte 

  write (*,'(a,i0,a)') "Each image has ",icells, " cells"
  write (*,'(a,i0,a)') "The model has ", mcells, " cells"
  write (*,'(a,i0,a)') "The model has ", nuc, " nuclei"
  write (*,'(a,es9.2,a)') "Each image will use at least ", &
    image_storage, " GB memory"
end if

! initialise random number seed
call cgca_irs(nodebug)

! allocate coarray
call cgca_as(l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3,1,space)

! initialise coarray
space = this_image()
solid = .true.

! allocate grain volume
call cgca_av(1,nuc,col1,cou1,col2,cou2,col3,grainvol)

! calculate volumes, sync all inside
call cgca_gv(space,grainvol)

! allocate rotation tensors
call cgca_art(1,nuc,col1,cou1,col2,cou2,col3,grt)

! assign rotation tensors, sync all inside
call cgca_rt(grt)
sync all

! check all rotation tensors on image 1
if (image1) then
  write (*,*) "rotation tensors assigned"
  do i=1,nuc
    call cgca_ckrt(grt(i,:,:),yesdebug,flag)
    if (flag .ne. 0) then
      write (*,*) "problem with grain: ", i
      write (*,*) "failed test: ", flag
      write (*,*) "stopping!"
      error stop
    end if
  end do
end if
sync all

! deallocate all arrays
call cgca_ds(space)
call cgca_dv(grainvol)
call cgca_drt(grt)
if (image1) write (*,*) "successfully deallocated rotation tensor coarray"

end program testAAQ

!*roboend*
