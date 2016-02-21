program z

use mpi
implicit none

integer ierr, num_procs, my_id

call MPI_INIT ( ierr )
call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

write (*,*) "Hello from ", my_id, " out of ", num_procs, " processes."
call execute_command_line( "hostname" )

call MPI_FINALIZE ( ierr )

end program z
