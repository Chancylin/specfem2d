subroutine calculate_accumulative_pnts(myrank, npnt_local, offset)

  use mpi

  use specfem_par, only: nproc 
  
  implicit none
  
  integer, intent(in) :: myrank, npnt_local
  integer, intent(inout):: offset
  integer :: offset_for_next
  integer :: tag_to_next,tag_from_prev,ier
  
  tag_to_next = 2*myrank + 1
  tag_from_prev = 2*myrank -1
  
  if( nproc == 1 ) then
     offset = 0
     return
  endif
  
  if( myrank == 0 )then
     !only send
     offset = 0 !I guess for this first partition, it is zero (or one)? 
     
     offset_for_next = offset + npnt_local
     call MPI_SEND (offset_for_next, 1, MPI_INTEGER, &
          myrank + 1, tag_to_next, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier )
  else if (myrank == nproc-1)then
     !only receive
     call MPI_RECV (offset, 1, MPI_INTEGER, &
          myrank -1, tag_from_prev , MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier )
  else
     !receive and send
     call MPI_RECV (offset, 1, MPI_INTEGER, &
          myrank -1, tag_from_prev , MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier )

     offset_for_next = offset + npnt_local
     call MPI_SEND (offset_for_next, 1, MPI_INTEGER, &
          myrank + 1, tag_to_next, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier )
  endif

  
end subroutine calculate_accumulative_pnts 
