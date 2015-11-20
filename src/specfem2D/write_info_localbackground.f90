
! by lcx: this file contains subroutines to store displacements,velocity information 
! of GLL points at local/background boundaries.
!this file basically follows write_seismograms.F90
!=================================  
  
  subroutine write_info_localbackground_to_file()
  
  use specfem_par

  implicit none

  integer :: i, j, iglob, irecloc, irec, ispec
  logical :: node1_located, node2_located

  node1_located = .false.
  node2_located = .false.
  
  if ( any_local_background_edges ) then
    do inum = 1, num_local_background_edges
       ispec = localbackground_local_ispec(inum) 
       node1 = localbackground_nodes(2*inum-1)
       node2 = localbackground_nodes(2*inum)
       !check whether node locates on z axis
       do i = 1, NGLLX, NGLLX-1
          do j = 1, NGLLZ
             if( coord(1,ibool(i,j,ispec)) == coorg(1,node1) &
                 .and. coord(2,ibool(i,j,ispec)) == coorg(2,node1) ) then
                node1_located = .true.
                node1_i = i
                node1_j = j
             endif

             if( coord(1,ibool(i,j,ispec)) == coorg(1,node2) &
                 .and. coord(2,ibool(i,j,ispec)) == coorg(2,node2) ) then
                node2_located = .true.
                node2_i = i
                node2_j = j
             endif
           enddo
        enddo

       !check whether node locates on x axis
       do j = 1, NGLLZ, NGLLZ-1
          do i = 2, NGLLX-1
             if( coord(1,ibool(i,j,ispec)) == coorg(1,node1) &
                 .and. coord(2,ibool(i,j,ispec)) == coorg(2,node1) ) then
                node1_located = .true.
                node1_i = i
                node1_j = j
             endif

             if( coord(1,ibool(i,j,ispec)) == coorg(1,node2) &
                 .and. coord(2,ibool(i,j,ispec)) == coorg(2,node2) ) then
                node2_located = .true.
                node2_i = i
                node2_j = j
             endif
           enddo
        enddo
        
        if ( (.not. node1_located) .or. (.not. node2_located) ) then
           stop 'error, nodes on local/background can not be located' 
        endif
       
        
       
       coorg(1,node1)
       coorg(2,node1) 
       coorg(1,node2)
       coorg(2,node2) 
    enddo
   endif

  



  end subroutine write_info_localbackground_to_file
