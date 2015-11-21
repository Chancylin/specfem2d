
! by lcx: this file contains subroutines to store displacements,velocity information 
! of GLL points at local/background boundaries.
!this file basically follows write_seismograms.F90
!=================================  
 
  subroutine check_nodesToGLL() 

  use specfem_par
  
  implicit none

  integer :: i, j, ispec, node1, node2
  logical :: node1_located, node2_located
 
  allocate(localbackground_node1_gll(2,num_local_background_edges))
  allocate(localbackground_node2_gll(2,num_local_background_edges))
  allocate(localbackground_edges_type(num_local_background_edges))
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
                localbackground_node1_gll(1,inum) = i
                localbackground_node1_gll(2,inum) = j
             endif

             if( coord(1,ibool(i,j,ispec)) == coorg(1,node2) &
                 .and. coord(2,ibool(i,j,ispec)) == coorg(2,node2) ) then
                node2_located = .true.
                localbackground_node2_gll(1,inum) = i
                localbackground_node2_gll(2,inum) = j
             endif
           enddo
        enddo

       !check whether node locates on x axis
       do j = 1, NGLLZ, NGLLZ-1
          do i = 2, NGLLX-1
             if( coord(1,ibool(i,j,ispec)) == coorg(1,node1) &
                 .and. coord(2,ibool(i,j,ispec)) == coorg(2,node1) ) then
                node1_located = .true.
                localbackground_node1_gll(1,inum) = i
                localbackground_node1_gll(2,inum) = j
             endif

             if( coord(1,ibool(i,j,ispec)) == coorg(1,node2) &
                 .and. coord(2,ibool(i,j,ispec)) == coorg(2,node2) ) then
                node2_located = .true.
                localbackground_node2_gll(1,inum) = i
                localbackground_node2_gll(2,inum) = j
             endif
           enddo
        enddo
        
        if ( (.not. node1_located) .or. (.not. node2_located) ) then
           stop 'error, nodes on local/background element can not be located' 
        endif
        node1_i = localbackground_node1_gll(1,inum)
        node1_j = localbackground_node1_gll(2,inum)
        node2_i = localbackground_node2_gll(1,inum)
        node2_j = localbackground_node2_gll(2,inum)
    
        if (node1_i == node2_i .and. node1_i ==1)then
           localbackground_edges_type(inum) = 1  !left
         else if (node1_i == node2_i .and. node1_i == NGLLX) then
           localbackground_edges_type(inum) = 2  !right
         else if (node1_j == node2_j .and. node1_j == 1) then
           localbackground_edges_type(inum) = 3  !bottom
         else if (node1_j == node2_j .and. node1_j == NGLLZ) then
           localbackground_edges_type(inum) = 4  !top
         else 
            stop 'error, cannot detect what kind of edge for local/background element'
        endif
       enddo
  endif
  end subroutine check_nodesToGLL 


  subroutine write_info_localbackground_to_file(t_num)
  
  use specfem_par

  implicit none

  integer :: i, j, iglob, ispec, f_num
  integer, intent(in) :: t_num
  integer :: node1_i, node1_j, node2_i, node2_j
  character (len=80) :: fname

  if ( any_local_background_edges ) then
    do inum = 1, num_local_background_edges
       ispec = localbackground_local_ispec(inum) 
       node1_i = localbackground_node1_gll(1,inum)
       node1_j = localbackground_node1_gll(2,inum)
       node2_i = localbackground_node2_gll(1,inum)
       node2_j = localbackground_node2_gll(2,inum)
      ! edge_type = localbackground_edges_type(inum)
        
        if ( node1_i == node2_i ) then
           do j = 1,NGLLZ
             iglob = ibool(node1_i,j,ispec)
             dxd = veloc_elastic(1,iglob) 
             dyd = veloc_elastic(2,iglob)
             dzd = veloc_elastic(3,iglob)
             write(fname, "('./OUTPUT_FILES/&
                   &localboundaryinfo/elmnt',i8.8,'/',&
                   &i1.1,'_',i1.1)") &
                   ispec, node1_i,j
             !!!this generate an unique file number to every gll point of 
             !!!every element
             f_num = 777 + ispec * 100 + node1_i * 10 + j
             open(unit=f_num,file=trim(fname),status='unknown',&
                  position='append',iostat=ios)
             if( ios /= 0 ) stop 'error saving local/background boundary nodes info'
             write(f_num,*) t_num,dxd,dyd,dzd
             if(t_num == NSTEP) close(f_num)
           enddo
        endif 

        if ( node1_j == node2_j ) then   
           do i = 1,NGLLX
             iglob = ibool(i,node1_j,ispec)
             dxd = veloc_elastic(1,iglob) 
             dyd = veloc_elastic(2,iglob)
             dzd = veloc_elastic(3,iglob)
             write(fname, "('./OUTPUT_FILES/&
                   &localboundaryinfo/elmnt',i8.8,'/',&
                   &i1.1,'_',i1.1)") &
                   ispec, i,node1_j
             f_num = 777 + ispec * 100 + i * 10 + node1_j
             open(unit=f_num,file=trim(fname),status='unknown',&
                  position='append',iostat=ios)
             if( ios /= 0 ) stop 'error saving local/background boundary nodes info'
             write(f_num,*) t_num,dxd,dyd,dzd
             if(t_num == NSTEP) close(f_num)
           enddo
        endif 
    enddo
   endif

  end subroutine write_info_localbackground_to_file
