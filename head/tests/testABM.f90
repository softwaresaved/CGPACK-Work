!$Id: testABM.f90 175 2015-12-15 12:31:30Z mexas $

!*robodoc*u* tests/testABM
!  NAME
!    testABM
!  SYNOPSIS

program testABM

!  PURPOSE
!    Testing MPI/IO, cgca_m2mpiio/cgca_pswci
!  DESCRIPTION
!    Timing output of MPI/IO (cgca_pswci) against
!    the serial version (cgca_swci).
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

real,parameter :: gigabyte=real(2**30), resolution=1.0e-5,             &
 loge2 = log(real(2))
logical(kind=ldef),parameter :: yesdebug = .true., nodebug = .false.

integer( kind=idef ) :: l1, u1, l2, u2, l3, u3, col1, cou1, col2,      &
 cou2,col3,cou3, &
 nuc,    & ! number of nuclei in the model
 nimages, codim(3)[*], p, img
integer( kind=iarr ), allocatable :: space(:,:,:,:)[:,:,:]
integer( kind=ilrg ) :: icells, mcells

real :: time1, time2, fsizeb, fsizeg, tdiff

!**********************************************************************73
! first executable statement

    img = this_image()
nimages = num_images()

! do a check on image 1
if ( img .eq. 1 ) then
 call getcodim(nimages,codim)
 ! print a banner
 call banner("ABM")
 ! print the parameter values
 call cgca_pdmp
 write (*,'(a,i0,a)') "running on ", nimages, " images in a 3D grid"
 ! dump the value of p
 write (*,"(a,i0)") "p=",p
end if

! all images read codim from image 1
sync all
codim(:) = codim(:)[1]

l1=1
l2=l1
l3=l1

! The array size is only controlled by this value
! in this program.
u1 = 2**6 ! 64
u2 = u1
u3 = u1

col1 = 1
col2 = col1
col3 = col1

cou1 = codim(1)-col1+1
cou2 = codim(2)-col2+1
cou3 = codim(3)-col3+1

! total number of cells in a coarray
icells = int(u1-l1+1,kind=ilrg) * int(u2-l2+1,kind=ilrg) *            &
 int(u3-l3+1,kind=ilrg)

! total number of cells in the model
mcells = icells * int(codim(1),kind=ilrg) * int(codim(2),kind=ilrg) * &
  int(codim(3),kind=ilrg)

! total number of nuclei
nuc = int( resolution * mcells )

if ( img .eq. 1 ) then
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "bounds: (",l1,u1,l2,u2,l3,u3
  write (*,'(a,2(i0,":",i0,","),i0,":",i0,")")') &
    "cobounds: (",col1,cou1,col2,cou2,col3,cou3

! Total output file size, in B and in GB.
fsizeb = real( mcells*storage_size(space,kind=ilrg)/8_ilrg )
fsizeg = fsizeb / gigabyte

  write (*,'(a,i0,a)') "Each image has ",icells, " cells"
  write (*,'(a,i0,a)') "The model has ", mcells, " cells"
  write (*,'(a,i0,a)') "The model has ", nuc, " nuclei"
  write (*,'(2(a,es9.2),a)') "The output file size is ", fsizeb, &
   " B, or ", fsizeg, "GB."
end if

! calculate file size in MB
! allocate coarray
call cgca_as(l1,u1,l2,u2,l3,u3,col1,cou1,col2,cou2,col3,1,space)

! initialise coarray to image number
space = int( img, kind=iarr )

! dump the model
call cpu_time(time1)
call cgca_pswci(space, cgca_state_type_grain, 'mpiio.raw')
call cpu_time(time2)
tdiff = time2-time1
if (img .eq. 1) write (*,*) "MPI/IO: ", tdiff, "s, rate: ", fsizeg/tdiff, "GB/s."

call cpu_time(time1)
call cgca_swci (space, cgca_state_type_grain, 10, 'serial.raw')
call cpu_time(time2)
tdiff = time2-time1
if (img .eq. 1) write (*,*) "Serial IO: ", tdiff, "s, rate: ", fsizeg/tdiff, "GB/s."

! deallocate all arrays
call cgca_ds(space)

end program testABM

!*roboend*
