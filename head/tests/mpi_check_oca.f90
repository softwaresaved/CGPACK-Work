program z

use mpi
implicit none

integer :: ierr, num_procs, my_id, length
character*(MPI_MAX_PROCESSOR_NAME) :: name

!call MPI_INIT ( ierr )
call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

call execute_command_line( "hostname" )
call MPI_GET_PROCESSOR_NAME( name, length, ierr )
write(*,*) my_id, " out of ", num_procs, " Processor name: ", trim(name)

!call MPI_FINALIZE ( ierr )

end program z
