!this subroutine will be called to record the boundary (local\global) info

 subroutine record_bd_elemnt_prediction()
 !this subroutine is used to store prediction value of velocity and pot_dot for the
 !selected bd recording element, since the values involved in the acceleration calculation
 !is prediction value but not correction value
   
    use specfem_par, only: ibool,veloc_elastic,potential_dot_acoustic,& !original para
                           nspec_bd_elmt_elastic_pure,nspec_bd_elmt_acoustic_pure,&
                           vel_bd_elastic, pot_dot_bd_acoustic,&
                           ispec_bd_elmt_elastic_pure,ispec_bd_elmt_acoustic_pure,&
                           it,record_nt1,record_nt2,& !control time step for recording
                           coord
   
    implicit none
    include "constants.h"

    integer :: i,j,k,iglob,ispec
    !for test
    double precision :: x,z
     
    if (it < record_nt1 .or. it > record_nt2 ) return

    loop1:do k = 1,nspec_bd_elmt_elastic_pure
             ispec = ispec_bd_elmt_elastic_pure(k)
             !check whether we find the correct location
             x = dble(coord(1,ibool(2,2,ispec)))
             z = dble(coord(2,ibool(2,2,ispec)))
             
             ! print *,'k = ', k, 'ispec = ', ispec, ' x = ', x, ' z = ', z, ' from rank ', myrank,&
             !      'veloc_elastic = ', veloc_elastic(3,ibool(2,2,ispec))

             do i = 1, NGLLX
                do j = 1, NGLLZ
                   !ibool could have different meaning in serial model and parallal mode
                   iglob = ibool(i,j,ispec)
                   vel_bd_elastic(1,i,j,k) = veloc_elastic(1,iglob)
                   vel_bd_elastic(2,i,j,k) = veloc_elastic(2,iglob)
                   vel_bd_elastic(3,i,j,k) = veloc_elastic(3,iglob)
                enddo
             enddo
     enddo loop1
  
     loop2:do k =1,nspec_bd_elmt_acoustic_pure
              ispec = ispec_bd_elmt_acoustic_pure(k)
              do i = 1, NGLLX 
                 do j = 1, NGLLZ
                    iglob = ibool(i,j,ispec)
                    pot_dot_bd_acoustic(i,j,k) = potential_dot_acoustic(iglob)
                 enddo
              enddo
     enddo loop2
 end subroutine record_bd_elemnt_prediction

 subroutine record_bd_elmnt_elastic(ispec,i,j,&
            sigma_xx,sigma_xy,sigma_xz,sigma_zz,sigma_zy)

   use specfem_par, only: stress_bd_elastic,nspec_bd_elmt_elastic_pure,&
                          ispec_bd_elmt_elastic_pure,&
                          it,record_nt1,record_nt2 !control time step for recording

   implicit none
   include "constants.h"
   real(kind=CUSTOM_REAL), intent(in) :: sigma_xx,sigma_xy,sigma_xz,sigma_zz,sigma_zy
   integer, intent(in) :: ispec,i,j
   integer :: k   

   if (it < record_nt1 .or. it > record_nt2 ) return
 
   loop1:do k = 1,nspec_bd_elmt_elastic_pure
         if ( ispec_bd_elmt_elastic_pure(k) == ispec ) then
         !record the information for this element
         !  trac_bd_pnt_elastic(1,i,j,k) = tx_store  
         !  trac_bd_pnt_elastic(2,i,j,k) = ty_store
         !  trac_bd_pnt_elastic(3,i,j,k) = tz_store
       
         !is it a better way to store stress tensor, so that we can get rid of the 
         !possible confusion of normal vector(nx_pnt, nz_pnt)?
           stress_bd_elastic(1,i,j,k) = sigma_xx 
           stress_bd_elastic(2,i,j,k) = sigma_xy
           stress_bd_elastic(3,i,j,k) = sigma_xz
           stress_bd_elastic(4,i,j,k) = sigma_zz
           stress_bd_elastic(5,i,j,k) = sigma_zy

           exit loop1
         endif
   enddo loop1 
 end subroutine record_bd_elmnt_elastic


 subroutine record_bd_elmnt_acoustic(ispec,i,j,dux_dxl,dux_dzl)
  
   use specfem_par, only: grad_pot_bd_acoustic,nspec_bd_elmt_acoustic_pure,&
                          ispec_bd_elmt_acoustic_pure,&
                          it,record_nt1,record_nt2 !control time step for recording
   implicit none
   include "constants.h"
   real(kind=CUSTOM_REAL), intent(in) :: dux_dxl,dux_dzl
   integer, intent(in) :: ispec,i,j
   integer :: k

   if (it < record_nt1 .or. it > record_nt2 ) return
   loop2:do k = 1,nspec_bd_elmt_acoustic_pure
         if ( ispec_bd_elmt_acoustic_pure(k) == ispec ) then

            grad_pot_bd_acoustic(1,i,j,k) = dux_dxl
            grad_pot_bd_acoustic(2,i,j,k) = dux_dzl
               
            !test
            !print *,'time step: ', it, '  recording ispec i j :', ispec_bd_elmt_acoustic_pure(k),i,j
            !print *,'dux_dxl = ', dux_dxl, ' dux_dzl = ', dux_dzl

            exit loop2
         endif
   enddo loop2
 end subroutine record_bd_elmnt_acoustic


 subroutine write_bd_pnts()

   use mpi
   
   use specfem_par, only: elastic,acoustic,it,myrank,& !original para
                          npnt_local,&
                          nspec_bd_pnt_elastic_clt, nspec_bd_pnt_acoustic_clt,&
                          ispec_selected_bd_pnt,f_num,& !fname
                          ispec_bd_elmt_elastic_pure, ispec_bd_elmt_acoustic_pure,&
                          hxi_bd_store, hgammar_bd_store,&
                          stress_bd_elastic,vel_bd_elastic,grad_pot_bd_acoustic,pot_dot_bd_acoustic,&
                          stress_bd_pnt_elastic,vel_bd_pnt_elastic,trac_bd_pnt_elastic,&
                          nspec_bd_elmt_elastic_pure,nspec_bd_elmt_acoustic_pure,&
                          nspec_bd_pnt_elastic,nspec_bd_pnt_acoustic,&
                          grad_pot_bd_pnt_acoustic,pot_dot_bd_pnt_acoustic,&
                          nx_bd_pnt_elastic,nz_bd_pnt_elastic,&
                          record_nt1,record_nt2,& !control time step for recording
                          bg_record_acoustic, bg_record_elastic,&
                          num_step_output,step_count,write_time_count,&
                          o_rank_elastic,o_rank_acoustic,&
                          num_to_gather_elastic,data_gather_count_elastic,gather_offset_elastic,&
                          num_to_gather_acoustic,data_gather_count_acoustic,gather_offset_acoustic,&
                          data_assemble_elastic,one_time_slice_elastic,data_to_save_elastic,&
                          data_assemble_acoustic,one_time_slice_acoustic,data_to_save_acoustic

  implicit none
  include "constants.h"

  integer :: ispec_bd_pnt_elastic, ispec_bd_pnt_acoustic
  integer :: ipnt,ispec,k,kk,i,j
  double precision :: hlagrange
  !integer :: ios
  !MPI parameters
  !integer :: offset1, offset2
  integer (kind=MPI_OFFSET_KIND) :: offset
  integer :: size_elastic,size_acoustic,bd_info_type,ierror
  integer (kind=MPI_OFFSET_KIND) :: count
  !integer (kind=MPI_COUNT_KIND) :: count
  ! real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: temp
  !integer :: k,kk
  !logical :: switch = .false.

  if (it < record_nt1 .or. it > record_nt2 ) return


  step_count = step_count + 1

  ispec_bd_pnt_elastic = 0
  ispec_bd_pnt_acoustic = 0
  do ipnt = 1, npnt_local

     ispec =  ispec_selected_bd_pnt(ipnt)

     if ( elastic(ispec) ) then

        ispec_bd_pnt_elastic = ispec_bd_pnt_elastic + 1

        loop1: do k = 1,nspec_bd_elmt_elastic_pure
           !here we locate in which element the recording point locates
           if ( ispec_bd_elmt_elastic_pure(k) == ispec ) then

              stress_bd_pnt_elastic(:,ispec_bd_pnt_elastic) = 0.0
              vel_bd_pnt_elastic(:,ispec_bd_pnt_elastic) = 0.0
              !we use the lagrange interpolation to calculate
              !the value at the point based on all the GLL points
              !in that element
              do i = 1,NGLLX
                 do j = 1,NGLLZ

                    !print *,'size of hlagrange is ', hlagrange
                    !print *,'size of stress_bd_pnt_elastic(:,ispec_bd_pnt_elastic) is ',&
                    !         shape(stress_bd_pnt_elastic(:,ispec_bd_pnt_elastic))
                    !print *,'size of stress_bd_elastic(:,i,j,k) is ', shape(stress_bd_elastic(:,i,j,k)) 

                    hlagrange = hxi_bd_store(ipnt,i)*hgammar_bd_store(ipnt,j)

                    stress_bd_pnt_elastic(:,ispec_bd_pnt_elastic) = &
                         stress_bd_pnt_elastic(:,ispec_bd_pnt_elastic) + &
                         stress_bd_elastic(:,i,j,k)*hlagrange 
                    !print *,'stress_bd_pnt_elastic(:,ispec_bd_pnt_elastic) is'
                    !print *,stress_bd_pnt_elastic(:,ispec_bd_pnt_elastic)
                    !print *,'stress_bd_elastic(:,i,j,k) is'
                    !print *,stress_bd_elastic(:,i,j,k)
                    !stop

                    vel_bd_pnt_elastic(:,ispec_bd_pnt_elastic) = &
                         vel_bd_pnt_elastic(:,ispec_bd_pnt_elastic) + &
                         vel_bd_elastic(:,i,j,k)*hlagrange
                 enddo
              enddo

              exit loop1
           endif
        enddo loop1

     endif!elastic pnt

     if ( acoustic(ispec) ) then

        ispec_bd_pnt_acoustic = ispec_bd_pnt_acoustic + 1

        loop2: do kk = 1,nspec_bd_elmt_acoustic_pure
           if ( ispec_bd_elmt_acoustic_pure(kk) == ispec ) then

              grad_pot_bd_pnt_acoustic(:,ispec_bd_pnt_acoustic) = 0.0
              pot_dot_bd_pnt_acoustic(ispec_bd_pnt_acoustic) = 0.0

              do i = 1,NGLLX
                 do j = 1,NGLLZ
                    hlagrange = hxi_bd_store(ipnt,i)*hgammar_bd_store(ipnt,j)

                    grad_pot_bd_pnt_acoustic(:,ispec_bd_pnt_acoustic) = &
                         grad_pot_bd_pnt_acoustic(:,ispec_bd_pnt_acoustic) + &
                         grad_pot_bd_acoustic(:,i,j,kk)*hlagrange

                    pot_dot_bd_pnt_acoustic(ispec_bd_pnt_acoustic) = &
                         pot_dot_bd_pnt_acoustic(ispec_bd_pnt_acoustic) + &
                         pot_dot_bd_acoustic(i,j,kk)*hlagrange
                 enddo
              enddo

              exit loop2
           endif
        enddo loop2

     endif!acoustic pnt

  enddo

  !now for the boundary info recording, I think the economic way is to save info
  !of all steps into a large binary file for each partition.  An alternative way
  !is to save the info of all partition into a large binary file for each time
  !step.


  !for elastic 

  ! write(fname,"('./OUTPUT_FILES/bg_record/&
  !      &elastic_pnts/nt_',i6.6)")it


  if( nspec_bd_pnt_elastic /= 0 ) then

     !store valuse at recording points at every time step
     !we firstly compute the traction for those recording points
     !one tricky things here is the normal vector. Since the normal vectors we have
     !store in the local model is outer_pointing, for traction exerted on local boundary 
     !elements, I suppose we should transfer to the opposite direction
     !trac_x
     trac_bd_pnt_elastic(1,:) = nx_bd_pnt_elastic(:)*stress_bd_pnt_elastic(1,:) + &
          nz_bd_pnt_elastic(:)*stress_bd_pnt_elastic(3,:)
     !trac_z
     trac_bd_pnt_elastic(3,:) = nx_bd_pnt_elastic(:)*stress_bd_pnt_elastic(3,:) + &
          nz_bd_pnt_elastic(:)*stress_bd_pnt_elastic(4,:)
     !trac_y
     trac_bd_pnt_elastic(2,:) = nx_bd_pnt_elastic(:)*stress_bd_pnt_elastic(2,:) + &
          nz_bd_pnt_elastic(:)*stress_bd_pnt_elastic(5,:)



#ifdef USE_MPI

     data_assemble_elastic(1:3,:) = trac_bd_pnt_elastic
     data_assemble_elastic(4:6,:) = vel_bd_pnt_elastic
     
     
     call MPI_SIZEOF(trac_bd_pnt_elastic(1,1),size_elastic,ierror)
     call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,size_elastic,bd_info_type,ierror)

     !gather info from all processor every time step
     call MPI_GATHERV(data_assemble_elastic, num_to_gather_elastic(o_rank_elastic+1)*6, bd_info_type, one_time_slice_elastic,&
          data_gather_count_elastic*6, gather_offset_elastic*6, bd_info_type, 0, bg_record_elastic, ierror)

     if( o_rank_elastic == 0 )then
        data_to_save_elastic(:,:,step_count) = one_time_slice_elastic
     endif

#else
!!!this is for serial writting
!!!this is the recording length for unformatted recording
     inquire (iolength = length_unf_1) trac_bd_pnt_elastic(:,1),vel_bd_pnt_elastic(:,1)

     !formatted recording
     !open(unit=f_num,file=trim(fname),status='new',&
     !     action='write',iostat=ios) 
     f_num=113

     !unformatted recording
     open(unit=f_num,file=trim(fname),access='direct',status='new',&
          action='write',iostat=ios,recl=length_unf_1) 
     if( ios /= 0 ) stop 'error saving values at recording points'

     do k = 1, nspec_bd_pnt_elastic
        !write(f_num,111) trac_bd_pnt_elastic(:,k),vel_bd_pnt_elastic(:,k)
        write(f_num,rec=k) trac_bd_pnt_elastic(:,k),vel_bd_pnt_elastic(:,k)
     enddo

     close(f_num)
#endif

  endif

  !for acoustic

  if( nspec_bd_pnt_acoustic /= 0 ) then


#ifdef USE_MPI

     data_assemble_acoustic(1:2,:) = grad_pot_bd_pnt_acoustic
     data_assemble_acoustic(3,:) = pot_dot_bd_pnt_acoustic


     call MPI_SIZEOF(grad_pot_bd_pnt_acoustic(1,1),size_acoustic,ierror)
     call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,size_acoustic,bd_info_type,ierror)

     call MPI_GATHERV(data_assemble_acoustic, num_to_gather_acoustic(o_rank_acoustic+1)*3, bd_info_type, one_time_slice_acoustic, &
          data_gather_count_acoustic*3, gather_offset_acoustic*3, bd_info_type, 0, bg_record_acoustic, ierror)

     if( o_rank_acoustic == 0 )then
        data_to_save_acoustic(:,:,step_count) = one_time_slice_acoustic
     endif


#else
     inquire (iolength = length_unf_2) grad_pot_bd_pnt_acoustic(:,1),pot_dot_bd_pnt_acoustic(1)

     !formatted recording
     !open(unit=f_num,file=trim(fname),status='new',&
     !     action='write',iostat=ios) 


     f_num=113
     !unformatted recording
     open(unit=f_num,file=trim(fname),access='direct',status='new',&
          action='write',iostat=ios,recl=length_unf_2) 
     if( ios /= 0 ) stop 'error saving values at recording points'

     do kk = 1, nspec_bd_pnt_acoustic
        !write(f_num,112) grad_pot_bd_pnt_acoustic(:,kk),pot_dot_bd_pnt_acoustic(kk)
        write(f_num,rec=kk) grad_pot_bd_pnt_acoustic(:,kk),pot_dot_bd_pnt_acoustic(kk)
     enddo

     close(f_num)
#endif


  endif


#ifdef USE_MPI

  if( step_count == num_step_output .or. it == record_nt2 ) then

     !self_write
     if( nspec_bd_pnt_elastic /= 0 .and. o_rank_elastic == 0 )then

        call MPI_FILE_OPEN(MPI_COMM_SELF,'./OUTPUT_FILES/bg_record/elastic_pnts_data',&
             MPI_MODE_CREATE+MPI_MODE_WRONLY,&
             MPI_INFO_NULL,f_num,ierror)

        count = step_count * nspec_bd_pnt_elastic_clt * 6_8
        offset = write_time_count*num_step_output &
                 *nspec_bd_pnt_elastic_clt*6_8*size_elastic

         print *, 'step_count = ', step_count, ' count = ', count, &
               ' offset = ', offset, ' from processor ', myrank

        call MPI_FILE_WRITE_AT(f_num, offset, data_to_save_elastic, count,&
             bd_info_type, MPI_STATUS_IGNORE, ierror)

        call MPI_FILE_CLOSE(f_num,ierror)

     endif

     if( nspec_bd_pnt_acoustic /= 0 .and. o_rank_acoustic == 0 )then

        call MPI_FILE_OPEN(MPI_COMM_SELF,'./OUTPUT_FILES/bg_record/acoustic_pnts_data',&
             MPI_MODE_CREATE+MPI_MODE_WRONLY,&
             MPI_INFO_NULL,f_num,ierror)

        count = step_count * nspec_bd_pnt_acoustic_clt * 3_8
        offset = write_time_count*num_step_output &
                 *nspec_bd_pnt_acoustic_clt*3_8*size_acoustic

        call MPI_FILE_WRITE_AT(f_num, offset, data_to_save_acoustic, count,&
             bd_info_type, MPI_STATUS_IGNORE, ierror)

        call MPI_FILE_CLOSE(f_num,ierror)
        
     endif
     
     write_time_count = write_time_count + 1

     !step counter is reset to zero for a new round of recording
     step_count = 0 

  endif
#endif

  !111 format(6(es12.4,2x)) !112 column
  !112 format(3(es12.4,2x)) !36 column
 end subroutine write_bd_pnts
