!$Id: ca_mpi.f90 8 2014-12-01 09:13:54Z mexas $

program ca_mpi

use mpi
implicit none

integer err, num_procs, rank
integer :: ca[*]

call MPI_INIT( err )
call MPI_COMM_RANK( MPI_COMM_WORLD, rank, err )
call MPI_COMM_SIZE( MPI_COMM_WORLD, num_procs, err )

write (*,"(4(a,i0))") "MPI rank ", rank , " : coar img ", this_image(), &
  " : MPI size ", num_procs, " : num_images ", num_images()

call MPI_FINALIZE ( err )

end program ca_mpi
