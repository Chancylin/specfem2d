subroutine build_commu_bg_record(nspec_bd_pnt_elastic, nspec_bd_pnt_acoustic)

  use mpi

  use specfem_par, only: bg_record_elastic, bg_record_acoustic

  implicit none

  integer, intent(in) :: nspec_bd_pnt_elastic, nspec_bd_pnt_acoustic
  integer :: ierror
  integer :: color, key = 1

  color = 1
  if( nspec_bd_pnt_elastic == 0 ) color = 0
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,bg_record_elastic,ierror)

  color = 1
  if( nspec_bd_pnt_acoustic == 0 ) color = 0
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,bg_record_acoustic,ierror)
  
end subroutine build_commu_bg_record
