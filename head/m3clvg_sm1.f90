!$Id: m3clvg_sm1.f90 216 2016-02-18 17:25:04Z mexas $

!*********************************************************************72

!*robodoc*f* cgca_m3clvg/m3clvg_sm1
!  NAME
!    m3clvg_sm1
!  SYNOPSIS

submodule ( cgca_m3clvg ) m3clvg_sm1

!  DESCRIPTION
!    Submodule of module cgca_m3clvg. It contains subroutines dealing
!    with updating the grain connectivity array.
!  AUTHOR
!    Anton Shterenlikht
!  COPYRIGHT
!    See LICENSE
!  CONTAINS
!    cgca_gcupd, cgca_gcupdn
!  USES
!    all variables and parameters of module cgca_m3clvg by host
!    association
!  USED BY
!    cgca_m3clvg
!  SOURCE

contains

!*roboend*


!*robodoc*s* m3clvg_sm1/cgca_gcupd
!  NAME
!    cgca_gcupd
!  SYNOPSIS  

module procedure cgca_gcupd

!  SIDE EFFECTS
!    State of GB array in module cgca_m2gb is updated
!  DESCRIPTION
!    This routine reads gcupd from *all* images and
!    adds the pairs to the local GB array on this image.
!    If you want to use just the nearest neighbouring images,
!    use cgca_gcupdn instead.
!    Synchronisation must be used before and after
!    calling this routine, to comply with the standard.
!  USES
!    cgca_gcf
!  SOURCE

integer( kind=idef ) :: i, j, img, nimgs, img_curr, rndint
integer( kind=kind(gcupd) ) :: gcupd_local( gculen, 3 )
real :: rnd

  img = this_image()
nimgs = num_images()

! choose the first image at random
call random_number( rnd )   ! [ 0 .. 1 )
rndint = int( rnd*nimgs )+1 ! [ 1 .. nimgs ]

! loop over all images, starting at a randomly chosen image
images: do j = rndint, rndint+nimgs-1

  ! Get the current image number.
  ! If it's > nimgs, subtract nimgs
  img_curr = j
  if ( img_curr .gt. nimgs ) img_curr = img_curr - nimgs

  ! Skip this image, because the GC array has already been updated
  ! in cgca_clvgsd.
  if ( img_curr .eq. img ) cycle images

  ! copy gcupd from image j into a local var
  gcupd_local( : , : ) = gcupd( : , : ) [img_curr]

  gcarray: do i = 1, gculen

    ! The gcupd array is filled with fractured pairs from the beginning
    ! so exit as soon as the GB state is intact.
    if ( gcupd_local( i , 3 ) .eq. cgca_gb_state_intact ) exit gcarray

!write (*,*) "DEBUG: cgca_gcupd: img:", img, &
!            "gcupd_local(i,:):", gcupd_local( i , : )

    ! add the pair to the GC array on this image
    call cgca_gcf( gcupd_local( i , 1 ), gcupd_local( i , 2 ) )

 end do gcarray
end do images

end procedure cgca_gcupd

!*roboend*


!*robodoc*s* m3clvg_sm1/cgca_gcupdn
!  NAME
!    cgca_gcupdn
!  SYNOPSIS  

module procedure cgca_gcupdn

!  INPUT
!    periodicbc - logical, .true. if the CA space has periodic BC,
!    and .false. otherwise.
!  SIDE EFFECTS
!    State of GB array in module cgca_m2gb is updated
!  DESCRIPTION
!    This routine reads gcupd_arr from the *nearest neighbouring*
!    images only, and adds the pairs to the local GB array on this
!    image. If you want to read from all images, use
!    cgca_gcupd.
!    Synchronisation must be used before and after
!    calling this routine, to comply with the standard.
!  NOTES
!    This routine must be used only after gcupd_arr has been
!    allocated. A runtime error will result if gcupd_arr has not
!    been allocated yet. 
!  USES
!    cgca_gcf
!  SOURCE

integer( kind=idef ) :: i, j, k, s, mycod( cgca_scodim ),              &
 neicod( cgca_scodim )
integer( kind=kind(gcupd_arr) ) ::                                     &
 gcupd_local( cgca_gcupd_size1, cgca_gcupd_size2 )

! Get my coindex set
if ( .not. allocated( gcupd_arr ) ) then
  write (*,'(a)')                                                      &
     "ERROR: cgca_m3clvg/cgca_gcupdn: gcupd_arr not allocated"
  error stop
end if
mycod = this_image( gcupd_arr )

! Loop over all nearest neighbours, taking special attention of
! the images at the edges of the model
do i = -1 , 1
do j = -1 , 1
inner: do k = -1 , 1

  ! Get the coindex set of the neighbour
  neicod = mycod + (/ i, j, k /)

  ! Skip this image
  if ( all( neicod .eq. mycod ) ) cycle inner

  ! Dealing with edges
  ! Loop over all codimensions
  edges: do s = 1 , cgca_scodim

    ! If the neighbour is below the lower edge
    if ( neicod( s ) .lt. cgca_slcob( s ) ) then
      if ( periodicbc ) then
        ! If periodic BC are in use, take the data from the opposite
        ! edge.
        neicod( s ) = cgca_sucob( s )       
      else
        ! Otherwise, do not pull data from this neighbour, move to
        ! the next one.
        cycle inner
      end if
    end if     

    ! If the neighbour is above the upper edge
    if ( neicod( s ) .gt. cgca_sucob( s ) ) then
      if ( periodicbc ) then
        ! If periodic BC are in use, take the data from the opposite
        ! edge.
        neicod( s ) = cgca_slcob( s )       
      else
        ! Otherwise, do not pull data from this neighbour, move to
        ! the next one.
        cycle inner
      end if
    end if     

  end do edges

  ! Now the coindex set of the neighbour has been obtained.
  ! Pull its data

  ! Copy gcupd_arr from the neighbouring image into a local var.
  ! Remote read.
  gcupd_local( : , : ) =                                               &
    gcupd_arr( : , : ) [ neicod(1), neicod(2), neicod(3) ]

  ! Scan all values in gcupd. Can reuse loop index s.
  gcarray: do s = 1, cgca_gcupd_size1

    ! The gcupd array is filled with fractured pairs from the beginning
    ! so exit as soon as the GB state is intact.
    if ( gcupd_local( s , 3 ) .eq. cgca_gb_state_intact ) exit gcarray

    !write (*,*) "DEBUG: cgca_gcupd: img:", img, &
    !            "gcupd_local(i,:):", gcupd_local( i , : )

    ! add the pair to the GC array on this image
    call cgca_gcf( gcupd_local( s , 1 ), gcupd_local( s , 2 ) )

  end do gcarray

end do inner
end do
end do

end procedure cgca_gcupdn

!*roboend*

end submodule m3clvg_sm1
