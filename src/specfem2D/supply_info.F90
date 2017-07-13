!this subroutine is called to import the boundary information needed
!for local simulation. Every time step it will be called
subroutine supply_bd_pnt()

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par, only: it,read_nt1,read_nt2,myrank,& !original para
                         nspec_bd_pnt_elastic_supply,nspec_bd_pnt_acoustic_supply,&
                         x_final_bd_pnt_elastic_supply,z_final_bd_pnt_elastic_supply,&
                         trac_bd_pnt_elastic,vel_bd_pnt_elastic,&
                         x_final_bd_pnt_acoustic_supply,z_final_bd_pnt_acoustic_supply,&
                         grad_pot_bd_pnt_acoustic,pot_dot_bd_pnt_acoustic,&
                         f_num,&!,fname
                         num_step_input,&
                         data_to_supply_elastic,one_time_slice_elastic_1,one_time_slice_elastic_2,&
                         trac_bd_pnt_t1,trac_bd_pnt_t2, vel_bd_pnt_t1, vel_bd_pnt_t2,&
                         data_to_supply_acoustic,one_time_slice_acoustic_1,one_time_slice_acoustic_2,&
                         grad_pot_bd_pnt_t1,grad_pot_bd_pnt_t2, pot_dot_bd_pnt_t1,pot_dot_bd_pnt_t2


  implicit none
  include "constants.h"

  integer :: i,ios,temp_read
  integer :: ierror
  character(len=150) dummystring

  !caution: the usage of variables could be plausible here. Since we now
  !consider the local and global simulation separately when implementing MPI,
  !the meaning of the set of variables here could be
  !different. nspec_bd_pnt_elastic, nspec_bd_pnt_acoustic,
  !x/z_final_bd_pnt_elastic_supply, x/z_final_bd_pnt_acoustic_supply are actually defined
  !according to the bd_info profiles of global model if one tries to record the
  !bd_info preparing for the following global simulation

  !now we use nspec_bd_pnt_elastic_supply instead of nspec_bd_pnt_elastic,
  !nspec_bd_pnt_acoustic_supply instead of nspec_bd_pnt_acoustic but I think we
  !could still use x/z_final_bd_pnt_elastic_supply, x/z_final_bd_pnt_acoustic_supply since
  !these arrays have been deallocate in subroutine locate_recording_ponit() and
  !could be reused.
  if ( it == 1) then
     nspec_bd_pnt_elastic_supply = 0
     nspec_bd_pnt_acoustic_supply = 0

#ifdef USE_MPI

     if( myrank == 0 )then

        print *,'rank 0 is doing its work'
        !count the total boundary points for 
        open(unit=1,file='./OUTPUT_FILES/bg_record/elastic_pnts_profile_total',iostat=ios,status='old',action='read')

        do while(ios == 0)
           read(1,"(a)",iostat=ios) dummystring
           if(ios == 0) nspec_bd_pnt_elastic_supply = nspec_bd_pnt_elastic_supply + 1
        enddo

        close(1)

        open(unit=1,file='./OUTPUT_FILES/bg_record/acoustic_pnts_profile_total',iostat=ios,status='old',action='read')

        do while(ios == 0)
           read(1,"(a)",iostat=ios) dummystring
           if(ios == 0) nspec_bd_pnt_acoustic_supply = nspec_bd_pnt_acoustic_supply + 1
        enddo

        close(1)

     endif

     ! print *, 'start broadcasting '
     call MPI_BCAST(nspec_bd_pnt_elastic_supply, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
     call MPI_BCAST(nspec_bd_pnt_acoustic_supply, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

     ! print *, 'nspec_bd_pnt_elastic_supply = ', nspec_bd_pnt_elastic_supply, ' from rank ', myrank
     ! print *, 'nspec_bd_pnt_acoustic_supply = ', nspec_bd_pnt_acoustic_supply, ' from rank ', myrank
#else

     !count the total boundary points for 
     open(unit=1,file='./OUTPUT_FILES/bg_record/elastic_pnts_profile',iostat=ios,status='old',action='read')
     do while(ios == 0)
        read(1,"(a)",iostat=ios) dummystring
        if(ios == 0) nspec_bd_pnt_elastic_supply = nspec_bd_pnt_elastic_supply + 1
     enddo
     close(1)

     open(unit=1,file='./OUTPUT_FILES/bg_record/acoustic_pnts_profile',iostat=ios,status='old',action='read')
     do while(ios == 0)
        read(1,"(a)",iostat=ios) dummystring
        if(ios == 0) nspec_bd_pnt_acoustic_supply = nspec_bd_pnt_acoustic_supply + 1
     enddo
     close(1)

#endif


     !every processor allocates the variables

     !for elastic element
     if( nspec_bd_pnt_elastic_supply /= 0 )then

        !allocate the arrays needed at the first time step 
        
        allocate(x_final_bd_pnt_elastic_supply(nspec_bd_pnt_elastic_supply),&
             z_final_bd_pnt_elastic_supply(nspec_bd_pnt_elastic_supply))
        allocate(trac_bd_pnt_elastic(3,nspec_bd_pnt_elastic_supply))
        allocate(vel_bd_pnt_elastic(3,nspec_bd_pnt_elastic_supply))

        trac_bd_pnt_elastic = 0.0
        vel_bd_pnt_elastic = 0.0

     endif

     !for acoustic element
     if( nspec_bd_pnt_acoustic_supply /= 0 )then


        allocate(x_final_bd_pnt_acoustic_supply(nspec_bd_pnt_acoustic_supply),&
             z_final_bd_pnt_acoustic_supply(nspec_bd_pnt_acoustic_supply))
        allocate(grad_pot_bd_pnt_acoustic(2,nspec_bd_pnt_acoustic_supply))
        allocate(pot_dot_bd_pnt_acoustic(nspec_bd_pnt_acoustic_supply))  

        grad_pot_bd_pnt_acoustic = 0.0
        pot_dot_bd_pnt_acoustic = 0.0
     endif

  endif

  !read the coordinate of interpolation points
  if ( it == 1) then

#ifdef USE_MPI

     if( myrank == 0 )then

        ! call MPI_FILE_OPEN(MPI_COMM_SELF, './OUTPUT_FILES/bg_record/elastic_pnts_profile', &
        !      MPI_MODE_RDONLY, MPI_INFO_NULL, f_num, ierror)
        ! call MPI_FILE_READ()

        f_num = 111

        if( nspec_bd_pnt_elastic_supply /= 0 )then

           allocate(data_to_supply_elastic(6,nspec_bd_pnt_elastic_supply,num_step_input))

           !two time points for time interpolation
           allocate(one_time_slice_elastic_1(6,nspec_bd_pnt_elastic_supply),one_time_slice_elastic_2(6,nspec_bd_pnt_elastic_supply))
           allocate(trac_bd_pnt_t1(3,nspec_bd_pnt_elastic_supply),vel_bd_pnt_t1(3,nspec_bd_pnt_elastic_supply))
           allocate(trac_bd_pnt_t2(3,nspec_bd_pnt_elastic_supply),vel_bd_pnt_t2(3,nspec_bd_pnt_elastic_supply))


           open(f_num,file='./OUTPUT_FILES/bg_record/elastic_pnts_profile_total',iostat=ios,status='old',action='read')

           do i=1,nspec_bd_pnt_elastic_supply
              read(f_num,110) temp_read,x_final_bd_pnt_elastic_supply(i), z_final_bd_pnt_elastic_supply(i)
           enddo

           close(f_num)

        endif

        if( nspec_bd_pnt_acoustic_supply /= 0 )then

           allocate(data_to_supply_acoustic(3,nspec_bd_pnt_acoustic_supply,num_step_input))

           allocate(one_time_slice_acoustic_1(3,nspec_bd_pnt_acoustic_supply),&
                one_time_slice_acoustic_2(3,nspec_bd_pnt_acoustic_supply))
           allocate(grad_pot_bd_pnt_t1(2,nspec_bd_pnt_acoustic_supply),pot_dot_bd_pnt_t1(nspec_bd_pnt_acoustic_supply))
           allocate(grad_pot_bd_pnt_t2(2,nspec_bd_pnt_acoustic_supply),pot_dot_bd_pnt_t2(nspec_bd_pnt_acoustic_supply))

           open(f_num,file='./OUTPUT_FILES/bg_record/acoustic_pnts_profile_total',iostat=ios,status='old',action='read')

           do i=1,nspec_bd_pnt_acoustic_supply
              read(f_num,110) temp_read,x_final_bd_pnt_acoustic_supply(i), z_final_bd_pnt_acoustic_supply(i)
           enddo

           close(f_num)

        endif

     endif

     !broadcast the value
     if( nspec_bd_pnt_elastic_supply /= 0 )then
        call MPI_BCAST(x_final_bd_pnt_elastic_supply, nspec_bd_pnt_elastic_supply, &
             MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
        call MPI_BCAST(z_final_bd_pnt_elastic_supply, nspec_bd_pnt_elastic_supply, &
             MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
     endif

     if( nspec_bd_pnt_acoustic_supply /= 0 )then
        call MPI_BCAST(x_final_bd_pnt_acoustic_supply, nspec_bd_pnt_acoustic_supply, &
             MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
        call MPI_BCAST(z_final_bd_pnt_acoustic_supply, nspec_bd_pnt_acoustic_supply, &
             MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
     endif

#else

     f_num = 111

     if( nspec_bd_pnt_elastic_supply /= 0 )then
        open(f_num,file='./OUTPUT_FILES/bg_record/elastic_pnts_profile_total',iostat=ios,status='old',action='read')
        do i=1,nspec_bd_pnt_elastic_supply
           read(f_num,110) temp_read,x_final_bd_pnt_elastic_supply(i), z_final_bd_pnt_elastic_supply(i)
        enddo
        close(f_num)

        print *,'nspec_bd_pnt_elastic_supply = ', nspec_bd_pnt_elastic_supply
     endif

     if( nspec_bd_pnt_acoustic_supply /= 0 )then
        open(f_num,file='./OUTPUT_FILES/bg_record/acoustic_pnts_profile_total',iostat=ios,status='old',action='read')
        do i=1,nspec_bd_pnt_acoustic_supply
           read(f_num,110) temp_read,x_final_bd_pnt_acoustic_supply(i), z_final_bd_pnt_acoustic_supply(i)
        enddo
        close(f_num)
        print *,'nspec_bd_pnt_acoustic_supply = ', nspec_bd_pnt_acoustic_supply
     endif

#endif

  endif

110 format(i5,2(es12.4,2x))!consistent with format in 'locate_recording_point.F90'
  !read the stored boundary info
  if (it < read_nt1 .or. it > read_nt2 ) return

  !apply the time interpolation 
  call time_interplt_supply() 

end subroutine supply_bd_pnt

subroutine time_interplt_supply()

#ifdef USE_MPI
  use mpi
#endif
  
  use specfem_par, only: it,deltat_read,myrank,&
                         record_nt1,record_nt2, deltat_record,&
                         nspec_bd_pnt_elastic_supply,nspec_bd_pnt_acoustic_supply,&
                         trac_bd_pnt_elastic,vel_bd_pnt_elastic,&
                         grad_pot_bd_pnt_acoustic,pot_dot_bd_pnt_acoustic,&
                         num_step_input,&
                         data_to_supply_elastic,one_time_slice_elastic_1,one_time_slice_elastic_2,&
                         trac_bd_pnt_t1,trac_bd_pnt_t2, vel_bd_pnt_t1, vel_bd_pnt_t2,&
                         data_to_supply_acoustic,one_time_slice_acoustic_1,one_time_slice_acoustic_2,&
                         grad_pot_bd_pnt_t1,grad_pot_bd_pnt_t2, pot_dot_bd_pnt_t1,pot_dot_bd_pnt_t2


  implicit none
  include "constants.h"

  integer :: nt1_record,nt2_record
  integer :: f_num!,i,f_num_2
  integer :: slice_1,slice_2

  double precision :: diff_deltat
  !define the array for mpi read
  ! character(len=150) :: newfile
#ifdef USE_MPI
  !integer :: size,bd_info_type_elastic,bd_info_type_acoustic,ierror
  integer :: bd_info_type,ierror,size
  integer(kind=MPI_OFFSET_KIND) :: count
  integer(kind=MPI_OFFSET_KIND) :: offset
#else
  integer :: length_unf_1
  integer :: length_unf_2
  integer :: ios
#endif
  
#ifdef USE_MPI
  if( myrank == 0 )then

     !this doesn't give the expected result, which is weird
     !nt1_record = floor(it * deltat_read / deltat_record)
     nt1_record = floor(it * (deltat_read / deltat_record))

     if(nt1_record < record_nt1 ) nt1_record = record_nt1
     if(nt1_record >=  record_nt2 ) nt1_record = record_nt2 - 1
     nt2_record = nt1_record + 1
     !nt2_record = ceiling(it * deltat_read / deltat_record)
     diff_deltat = 0.5 * (deltat_record - deltat_read)

  endif
  
  call MPI_BCAST(nt1_record, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call MPI_BCAST(nt2_record, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call MPI_BCAST(diff_deltat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  
#else
  
  nt1_record = floor(it * deltat_read / deltat_record)
  if(nt1_record < record_nt1 ) nt1_record = record_nt1
  if(nt1_record >=  record_nt2 ) nt1_record = record_nt2 - 1
  nt2_record = nt1_record + 1
  !nt2_record = ceiling(it * deltat_read / deltat_record)
  diff_deltat = 0.5 * (deltat_record - deltat_read)
#endif
  

  !calculate the time interpolation
#ifdef USE_MPI
 
  
  if( myrank == 0) then

     if( nspec_bd_pnt_elastic_supply /= 0 ) then

        !during some specific time steps, read a pool of data
        if( mod((nt1_record-record_nt1),num_step_input) == 0 )then
           
           call MPI_SIZEOF(trac_bd_pnt_t1(1,1),size,ierror)
           call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,size,bd_info_type,ierror)

           if( (record_nt2 - nt1_record) < num_step_input )then
              count = record_nt2 - nt1_record + 1
           else
              count = num_step_input
           endif

           offset = (nt1_record-record_nt1)*nspec_bd_pnt_elastic_supply*6_8*size

           print *,'elastic: do data_reading at step, it = ', it, ' nt1_record = ', nt1_record,&
                ' offset = ', offset, ' count = ', count

           call MPI_FILE_OPEN(MPI_COMM_SELF, './OUTPUT_FILES/bg_record/elastic_pnts_data', &
                MPI_MODE_RDONLY, MPI_INFO_NULL, f_num, ierror)
           
           
           call MPI_FILE_READ_AT(f_num, offset, data_to_supply_elastic , 6*nspec_bd_pnt_elastic_supply*count,&
                bd_info_type, MPI_STATUS_IGNORE, ierror)

           call MPI_FILE_CLOSE(f_num,ierror)

        endif
        

        !read two time slices for interpolation
        if( mod((nt1_record - record_nt1), num_step_input) /= (num_step_input -1) )then

           slice_1 = mod((nt1_record - record_nt1), num_step_input) + 1
        else
           slice_1 = mod((nt1_record - record_nt1), num_step_input)
        endif
        
        slice_2 = slice_1 + 1

        one_time_slice_elastic_1 = data_to_supply_elastic(:,:,slice_1)
        one_time_slice_elastic_2 = data_to_supply_elastic(:,:,slice_2)

        trac_bd_pnt_t1 = one_time_slice_elastic_1(1:3,:)
        vel_bd_pnt_t1 = one_time_slice_elastic_1(4:6,:)

        trac_bd_pnt_t2 = one_time_slice_elastic_2(1:3,:)
        vel_bd_pnt_t2 = one_time_slice_elastic_2(4:6,:)


        ! allocate(trac_bd_pnt_t1(3,nspec_bd_pnt_elastic_supply),vel_bd_pnt_t1(3,nspec_bd_pnt_elastic_supply))
        ! allocate(trac_bd_pnt_t2(3,nspec_bd_pnt_elastic_supply),vel_bd_pnt_t2(3,nspec_bd_pnt_elastic_supply))

        ! call MPI_SIZEOF(trac_bd_pnt_t1(1,1),size,ierror)
        ! call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,size,bd_info_type,ierror)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         ! write(fname_1,"('./OUTPUT_FILES/bg_record/&
!         !      &elastic_pnts/nt_',i6.6)")nt1_record

!         ! print *,'t1 = ', nt1_record, ' t2 = ', nt2_record
        
!         call MPI_FILE_OPEN(MPI_COMM_SELF, './OUTPUT_FILES/bg_record/elastic_pnts_data', &
!              MPI_MODE_RDONLY, MPI_INFO_NULL, f_num_1, ierror)

!         offset_time = (nt1_record - record_nt1)*nspec_bd_pnt_elastic_supply*size*6_8
        
!         ! print *,'t1, read traction at ', offset_time
!         call MPI_FILE_READ_AT(f_num_1, offset_time, trac_bd_pnt_t1, 3*nspec_bd_pnt_elastic_supply,&
!              bd_info_type, MPI_STATUS_IGNORE, ierror)
!         ! call MPI_FILE_GET_POSITION(f_num_1, offset, ierror)
!         ! print *, '1st reading, offset: ', offset, ' from rank ', myrank

!         offset_time = offset_time + 3*nspec_bd_pnt_elastic_supply*size
        
!         ! print *,'t1, read velocity at ', offset_time
!         call MPI_FILE_READ_AT(f_num_1, offset_time, vel_bd_pnt_t1, 3*nspec_bd_pnt_elastic_supply,&
!              bd_info_type, MPI_STATUS_IGNORE, ierror)
!         ! call MPI_FILE_GET_POSITION(f_num_1, offset, ierror)
!         ! print *, '2nd reading, offset: ', offset, ' from rank ', myrank

!         ! call MPI_FILE_CLOSE(f_num_1,ierror)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         ! write(fname_2,"('./OUTPUT_FILES/bg_record/&
!         !      &elastic_pnts/nt_',i6.6)")nt2_record

!         ! call MPI_FILE_OPEN(MPI_COMM_SELF, './OUTPUT_FILES/bg_record/elastic_pnts_data', &
!         !      MPI_MODE_RDONLY, MPI_INFO_NULL, f_num_2, ierror)
        
!         offset_time = (nt2_record - record_nt1)*nspec_bd_pnt_elastic_supply*size*6

!         ! print *,'t2, read traction at ', offset_time
!         call MPI_FILE_READ_AT(f_num_1, offset_time, trac_bd_pnt_t2, 3*nspec_bd_pnt_elastic_supply,&
!              bd_info_type, MPI_STATUS_IGNORE, ierror)

!         offset_time = offset_time + 3*nspec_bd_pnt_elastic_supply*size

!         ! print *,'t2, read velocity at ', offset_time
!         call MPI_FILE_READ_AT(f_num_1, offset_time, vel_bd_pnt_t2, 3*nspec_bd_pnt_elastic_supply,&
!              bd_info_type, MPI_STATUS_IGNORE, ierror)

!         call MPI_FILE_CLOSE(f_num_1,ierror)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !using another faster way, loop over the component, but not points

!!!linear interpolation in time space
        trac_bd_pnt_elastic = (trac_bd_pnt_t2 - trac_bd_pnt_t1) * &
             (it*deltat_read - nt1_record*deltat_record)/deltat_record + &
             trac_bd_pnt_t1
        !note that velocity is recorded as prediction at half time step
        vel_bd_pnt_elastic = (vel_bd_pnt_t2 - vel_bd_pnt_t1) * &
             (it*deltat_read - nt1_record*deltat_record + diff_deltat)/deltat_record + &
             vel_bd_pnt_t1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     endif

     if( nspec_bd_pnt_acoustic_supply /= 0 ) then


        if( mod((nt1_record-record_nt1),num_step_input) == 0 )then
           
           call MPI_SIZEOF(grad_pot_bd_pnt_t1(1,1),size,ierror)
           call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,size,bd_info_type,ierror)

           if( (record_nt2 - nt1_record) < num_step_input )then
              count = record_nt2 - nt1_record + 1
           else
              count = num_step_input
           endif

           offset = (nt1_record-record_nt1)*nspec_bd_pnt_acoustic_supply*3_8*size

           print *,'acoustic: do data_reading at step, it = ', it, ' nt1_record = ', nt1_record,&
                ' offset = ', offset, ' count = ', count

           call MPI_FILE_OPEN(MPI_COMM_SELF, './OUTPUT_FILES/bg_record/acoustic_pnts_data', &
                MPI_MODE_RDONLY, MPI_INFO_NULL, f_num, ierror)

           call MPI_FILE_READ_AT(f_num, offset, data_to_supply_acoustic, 3*nspec_bd_pnt_acoustic_supply*count,&
                bd_info_type, MPI_STATUS_IGNORE, ierror)

           call MPI_FILE_CLOSE(f_num,ierror)
           
        endif
        
        !read two time slices for interpolation
        if( mod((nt1_record - record_nt1), num_step_input) /= (num_step_input -1) )then

           slice_1 = mod((nt1_record - record_nt1), num_step_input) + 1
        else
           slice_1 = mod((nt1_record - record_nt1), num_step_input)
        endif
        
        slice_2 = slice_1 + 1

        one_time_slice_acoustic_1 = data_to_supply_acoustic(:,:,slice_1)
        one_time_slice_acoustic_2 = data_to_supply_acoustic(:,:,slice_2)

        grad_pot_bd_pnt_t1 = one_time_slice_acoustic_1(1:2,:)
        pot_dot_bd_pnt_t1 = one_time_slice_acoustic_1(3,:)

        grad_pot_bd_pnt_t2 = one_time_slice_acoustic_2(1:2,:)
        pot_dot_bd_pnt_t2 = one_time_slice_acoustic_2(3,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         ! write(fname_1,"('./OUTPUT_FILES/bg_record/&
!         !      &acoustic_pnts/nt_',i6.6)")nt1_record

!         offset_time = (nt1_record - record_nt1)*nspec_bd_pnt_acoustic_supply*size*3_8

       
!         call MPI_FILE_READ_AT(f_num_1, offset_time, grad_pot_bd_pnt_t1, 2*nspec_bd_pnt_acoustic_supply,&
!              bd_info_type, MPI_STATUS_IGNORE, ierror)
        
!         offset_time = offset_time + nspec_bd_pnt_acoustic_supply*size*2
!         call MPI_FILE_READ_AT(f_num_1, offset_time, pot_dot_bd_pnt_t1, nspec_bd_pnt_acoustic_supply,&
!              bd_info_type, MPI_STATUS_IGNORE, ierror)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         ! write(fname_2,"('./OUTPUT_FILES/bg_record/&
!         !      &acoustic_pnts/nt_',i6.6)")nt2_record

!         ! call MPI_FILE_OPEN(MPI_COMM_SELF, fname_2, &
!         !      MPI_MODE_RDONLY, MPI_INFO_NULL, f_num_2, ierror)

!         offset_time = (nt2_record - record_nt1)*nspec_bd_pnt_acoustic_supply*size*3_8
!         call MPI_FILE_READ_AT(f_num_1, offset_time, grad_pot_bd_pnt_t2, 2*nspec_bd_pnt_acoustic_supply,&
!              bd_info_type, MPI_STATUS_IGNORE, ierror)

!         offset_time = offset_time + nspec_bd_pnt_acoustic_supply*size*2
!         call MPI_FILE_READ_AT(f_num_1, offset_time, pot_dot_bd_pnt_t2, nspec_bd_pnt_acoustic_supply,&
!              bd_info_type, MPI_STATUS_IGNORE, ierror)

!         call MPI_FILE_CLOSE(f_num_1,ierror)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! write(newfile,"('./OUTPUT_FILES/bg_record/&
        !      &acoustic_pnts/nt_',i6.6)")nt1_record
        
        ! open(unit=2,file=trim(newfile),status='unknown',action='write')
        
        ! do i= 1,nspec_bd_pnt_acoustic_supply

        !    ! write(2,111) grad_pot_bd_pnt_t1(:,i), pot_dot_bd_pnt_t1(i)
           
        !    grad_pot_bd_pnt_acoustic(:,i) = (grad_pot_bd_pnt_t2(:,i) - grad_pot_bd_pnt_t1(:,i)) * &
        !         (it*deltat_read - nt1_record*deltat_record)/deltat_record + &
        !         grad_pot_bd_pnt_t1(:,i)

        !    pot_dot_bd_pnt_acoustic(i) = (pot_dot_bd_pnt_t2(i) - pot_dot_bd_pnt_t1(i)) * &
        !         (it*deltat_read - nt1_record*deltat_record + diff_deltat)/deltat_record + &
        !         pot_dot_bd_pnt_t1(i)

        ! enddo

        ! close(2)
        ! 111 format(6(es12.4,2x)) !112 column
        !need to check whether this form could work (without loop)
        grad_pot_bd_pnt_acoustic = (grad_pot_bd_pnt_t2 - grad_pot_bd_pnt_t1) * &
             (it*deltat_read - nt1_record*deltat_record)/deltat_record + &
             grad_pot_bd_pnt_t1

        
        pot_dot_bd_pnt_acoustic(:) = (pot_dot_bd_pnt_t2(:) - pot_dot_bd_pnt_t1(:)) * &
             (it*deltat_read - nt1_record*deltat_record + diff_deltat)/deltat_record + &
             pot_dot_bd_pnt_t1(:)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     endif

     
  endif !endif for myrank == 0


  !broadcast the result
  if( nspec_bd_pnt_elastic_supply /= 0 ) then
     
     !call MPI_BCAST(bd_info_type_elastic, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

     call MPI_SIZEOF(trac_bd_pnt_elastic(1,1),size,ierror)
     call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,size,bd_info_type,ierror)

     call MPI_BCAST(trac_bd_pnt_elastic, 3*nspec_bd_pnt_elastic_supply, bd_info_type, 0, MPI_COMM_WORLD, ierror)
     call MPI_BCAST(vel_bd_pnt_elastic, 3*nspec_bd_pnt_elastic_supply, bd_info_type, 0, MPI_COMM_WORLD, ierror)

  endif
  
  if( nspec_bd_pnt_acoustic_supply /= 0 ) then


     !call MPI_BCAST(bd_info_type_acoustic, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

     call MPI_SIZEOF(grad_pot_bd_pnt_acoustic(1,1),size,ierror)
     call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,size,bd_info_type,ierror)
     
     call MPI_BCAST(grad_pot_bd_pnt_acoustic, 2*nspec_bd_pnt_acoustic_supply, bd_info_type, 0, MPI_COMM_WORLD, ierror)
     call MPI_BCAST(pot_dot_bd_pnt_acoustic, nspec_bd_pnt_acoustic_supply, bd_info_type, 0, MPI_COMM_WORLD, ierror)

  endif
  
#else
  !for the unformatted recording, now we use array reading instead element by
  !element, line by line. The sequence could be wrong, need to be carefully
  !checked against the way we save the file (without MPI implementation)
  
  !!!this is the recording length for unformatted recording
  !inquire (iolength = length_unf_1) trac_bd_pnt_elastic(:,1),vel_bd_pnt_elastic(:,1)
  !inquire (iolength = length_unf_2) grad_pot_bd_pnt_acoustic(:,1),pot_dot_bd_pnt_acoustic(1)
  if( nspec_bd_pnt_elastic_supply /= 0 ) then

     allocate(trac_bd_pnt_t1(3,nspec_bd_pnt_elastic_supply),vel_bd_pnt_t1(3,nspec_bd_pnt_elastic_supply))
     allocate(trac_bd_pnt_t2(3,nspec_bd_pnt_elastic_supply),vel_bd_pnt_t2(3,nspec_bd_pnt_elastic_supply))

    inquire (iolength = length_unf_1) trac_bd_pnt_t1(:,1),vel_bd_pnt_t1(:,1)

    !elstic elements
    f_num_1=113
    write(fname_1,"('./OUTPUT_FILES/bg_record/&
          &elastic_pnts/nt_',i6.6)")nt1_record
    

    !unformatted reading
    open(unit=f_num_1,file=trim(fname_1),access='direct',status='old',&
         action='read',iostat=ios,recl=length_unf_1)

    if( ios /= 0 ) stop 'error reading values at profile points' 
    
    f_num_2=114
    write(fname_2,"('./OUTPUT_FILES/bg_record/&
          &elastic_pnts/nt_',i6.6)")nt2_record

    !unformatted reading
    open(unit=f_num_2,file=trim(fname_2),access='direct',status='old',&
         action='read',iostat=ios,recl=length_unf_1)

    if( ios /= 0 ) stop 'error reading values at profile points' 


    do i=1,nspec_bd_pnt_elastic_supply
       read(f_num_1,rec=i) trac_bd_pnt_t1(:,i),vel_bd_pnt_t1(:,i)
       read(f_num_2,rec=i) trac_bd_pnt_t2(:,i),vel_bd_pnt_t2(:,i)

       !!!linear interpolation in time space
       trac_bd_pnt_elastic(:,i) = (trac_bd_pnt_t2(:,i) - trac_bd_pnt_t1(:,i)) * &
                                  (it*deltat_read - nt1_record*deltat_record)/deltat_record + &
                                  trac_bd_pnt_t1(:,i)
       !note that velocity is recorded as prediction at half time step
       vel_bd_pnt_elastic(:,i) = (vel_bd_pnt_t2(:,i) - vel_bd_pnt_t1(:,i)) * &
                                 (it*deltat_read - nt1_record*deltat_record + diff_deltat)/deltat_record + &
                                 vel_bd_pnt_t1(:,i)
       
    enddo
 
    close(f_num_1)
    close(f_num_2)

    deallocate(trac_bd_pnt_t1,vel_bd_pnt_t1)
    deallocate(trac_bd_pnt_t2,vel_bd_pnt_t2)
  endif


  !acoustic elements
  if( nspec_bd_pnt_acoustic_supply /= 0 ) then 

     allocate(grad_pot_bd_pnt_t1(2,nspec_bd_pnt_acoustic_supply),pot_dot_bd_pnt_t1(nspec_bd_pnt_acoustic_supply))
     allocate(grad_pot_bd_pnt_t2(2,nspec_bd_pnt_acoustic_supply),pot_dot_bd_pnt_t2(nspec_bd_pnt_acoustic_supply))

     inquire (iolength = length_unf_2) grad_pot_bd_pnt_t1(:,1),pot_dot_bd_pnt_t1(1)

     write(fname_1,"('./OUTPUT_FILES/bg_record/&
           &acoustic_pnts/nt_',i6.6)")nt1_record

     !unformatted reading
     open(unit=f_num_1,file=trim(fname_1),access='direct',status='old',&
          action='read',iostat=ios,recl=length_unf_2)

     if( ios /= 0 ) stop 'error reading values at profile points' 
     
     write(fname_2,"('./OUTPUT_FILES/bg_record/&
           &acoustic_pnts/nt_',i6.6)")nt2_record

     !unformatted reading
     open(unit=f_num_2,file=trim(fname_2),access='direct',status='old',&
          action='read',iostat=ios,recl=length_unf_2)

     if( ios /= 0 ) stop 'error reading values at profile points' 


     do i= 1,nspec_bd_pnt_acoustic_supply
        read(f_num_1,rec=i) grad_pot_bd_pnt_t1(:,i),pot_dot_bd_pnt_t1(i)
        read(f_num_2,rec=i) grad_pot_bd_pnt_t2(:,i),pot_dot_bd_pnt_t2(i)

        grad_pot_bd_pnt_acoustic(:,i) = (grad_pot_bd_pnt_t2(:,i) - grad_pot_bd_pnt_t1(:,i)) * &
                                        (it*deltat_read - nt1_record*deltat_record)/deltat_record + &
                                        grad_pot_bd_pnt_t1(:,i)

        pot_dot_bd_pnt_acoustic(i) = (pot_dot_bd_pnt_t2(i) - pot_dot_bd_pnt_t1(i)) * &
                                     (it*deltat_read - nt1_record*deltat_record + diff_deltat)/deltat_record + &
                                     pot_dot_bd_pnt_t1(i)

     enddo
     
     close(f_num_1)
     close(f_num_2)
  endif

  deallocate(grad_pot_bd_pnt_t1,pot_dot_bd_pnt_t1)
  deallocate(grad_pot_bd_pnt_t2,pot_dot_bd_pnt_t2)
#endif

end subroutine time_interplt_supply


subroutine supply_pnt_reconst()

  use specfem_par, only: it, read_nt1_reconst, read_nt2_reconst
  

  implicit none
  include "constants.h"
  
   
  if (it < read_nt1_reconst .or. it > read_nt2_reconst ) return

  !apply the time interpolation 
  call time_interplt_supply_reconst() 

end subroutine supply_pnt_reconst



subroutine time_interplt_supply_reconst()

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par, only: it,p_sv,deltat_read_reconst,&
                         record_nt1_reconst,record_nt2_reconst, deltat_record_reconst,&
                         supply_reconst_elastic,supply_reconst_acoustic,&
                         o_rank_elastic,o_rank_acoustic,&
                         nspec_bd_pnt_elastic_supply,nspec_bd_pnt_acoustic_supply, &
                         trac_f,m_xx,m_xz,m_zz,m_yx,m_yz,&
                         Grad_pot,Pot_x,Pot_z,&
                         
                         nspec_bd_pnt_elastic_supply_total,nspec_bd_pnt_acoustic_supply_total,&
                         num_step_input,&
                         data_to_supply_elastic,&
                         one_time_slice_elastic_1,one_time_slice_elastic_2,&
                         trac_f_t1,trac_f_t2,trac_f_total,&
                         m_xx_t1,m_xz_t1,m_zz_t1,m_yx_t1,m_yz_t1,&
                         m_xx_t2,m_xz_t2,m_zz_t2,m_yx_t2,m_yz_t2,&
                         m_xx_total,m_xz_total,m_zz_total,m_yx_total,m_yz_total,&
                         data_to_supply_acoustic,&
                         one_time_slice_acoustic_1,one_time_slice_acoustic_2,&
                         Grad_pot_t1,Pot_x_t1,Pot_z_t1,&
                         Grad_pot_t2,Pot_x_t2,Pot_z_t2,&
                         Grad_pot_total,Pot_x_total,Pot_z_total,&
                         booking_reconst_elastic,booking_reconst_acoustic
                         

  implicit none
  include "constants.h"

  integer :: nt1_record_reconst,nt2_record_reconst
  integer :: f_num
 
  integer (kind=MPI_OFFSET_KIND) :: offset
  integer (kind=MPI_OFFSET_KIND) :: count, count_type
  integer :: bd_info_type,ierror,size
  integer :: slice_1,slice_2
  
#ifdef USE_MPI
#else
  integer :: i,ios 
  integer :: length_unf_1
  integer :: length_unf_2
#endif
  
  ! trac_f_t1 = 0.0
  ! trac_f_t2 = 0.0

  
  nt1_record_reconst = floor(it * (deltat_read_reconst / deltat_record_reconst))
  if(nt1_record_reconst < record_nt1_reconst ) nt1_record_reconst = record_nt1_reconst
  if(nt1_record_reconst >=  record_nt2_reconst ) nt1_record_reconst = record_nt2_reconst - 1
  nt2_record_reconst = nt1_record_reconst + 1
  !nt2_record_reconst = ceiling(it * deltat_read_reconst / deltat_record_reconst)

  ! print *,'it = ', it, ' nt1_record_reconst = ', nt1_record_reconst
  !calculate the time interpolation
#ifdef USE_MPI
     
  if( nspec_bd_pnt_elastic_supply /= 0 .and. o_rank_elastic == 0 ) then

     !do the data input
     if( mod((nt1_record_reconst - record_nt1_reconst), num_step_input) == 0 )then

        call MPI_SIZEOF(trac_f_t1(1,1),size,ierror)
        call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,size,bd_info_type,ierror)

        if( (record_nt2_reconst - nt1_record_reconst ) < num_step_input )then
           count = record_nt2_reconst - nt1_record_reconst
        else
           count = num_step_input
        endif

        if( p_sv )then

           offset = (nt1_record_reconst - record_nt1_reconst)*nspec_bd_pnt_elastic_supply_total*5_8*size
           count_type = 5_8*nspec_bd_pnt_elastic_supply_total*count

        else

           offset = (nt1_record_reconst - record_nt1_reconst)*nspec_bd_pnt_elastic_supply_total*3_8*size
           count_type = 3_8*nspec_bd_pnt_elastic_supply_total*count

        endif

        ! print *,'offset = ', offset, ' count_type = ', count_type

        call MPI_FILE_OPEN(MPI_COMM_SELF, './OUTPUT_FILES/reconst_record/elastic_pnts_data', &
             MPI_MODE_RDONLY, MPI_INFO_NULL, f_num, ierror)

        call MPI_FILE_READ_AT(f_num, offset, data_to_supply_elastic, count_type,&
             bd_info_type, MPI_STATUS_IGNORE, ierror)

        call MPI_FILE_CLOSE(f_num,ierror)

     endif

     !read two time slices for interpolation
     if( mod((nt1_record_reconst - record_nt1_reconst), num_step_input) /= (num_step_input -1) )then

        slice_1 = mod((nt1_record_reconst - record_nt1_reconst), num_step_input) + 1
     else
        slice_1 = mod((nt1_record_reconst - record_nt1_reconst), num_step_input)
     endif

     slice_2 = slice_1 + 1

     one_time_slice_elastic_1 = data_to_supply_elastic(:,:,slice_1)
     one_time_slice_elastic_2 = data_to_supply_elastic(:,:,slice_2)

     if( p_sv )then

        trac_f_t1(1,:) = one_time_slice_elastic_1(1,:) 
        trac_f_t1(3,:) = one_time_slice_elastic_1(2,:) 
        m_xx_t1 = one_time_slice_elastic_1(3,:) 
        m_xz_t1 = one_time_slice_elastic_1(4,:) 
        m_zz_t1 = one_time_slice_elastic_1(5,:)

        trac_f_t2(1,:) = one_time_slice_elastic_2(1,:) 
        trac_f_t2(3,:) = one_time_slice_elastic_2(2,:) 
        m_xx_t2 = one_time_slice_elastic_2(3,:) 
        m_xz_t2 = one_time_slice_elastic_2(4,:) 
        m_zz_t2 = one_time_slice_elastic_2(5,:)  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!linear interpolation in time domain
        trac_f_total(1,:) = (trac_f_t2(1,:) - trac_f_t1(1,:)) * &
             (it*deltat_read_reconst - nt1_record_reconst*deltat_record_reconst)/deltat_record_reconst + &
             trac_f_t1(1,:)

        trac_f_total(3,:) = (trac_f_t2(3,:) - trac_f_t1(3,:)) * &
             (it*deltat_read_reconst - nt1_record_reconst*deltat_record_reconst)/deltat_record_reconst + &
             trac_f_t1(3,:)

        m_xx_total = (m_xx_t2 - m_xx_t1) * &
             (it*deltat_read_reconst - nt1_record_reconst*deltat_record_reconst)/deltat_record_reconst + &
             m_xx_t1
        m_xz_total = (m_xz_t2 - m_xz_t1) * &
             (it*deltat_read_reconst - nt1_record_reconst*deltat_record_reconst)/deltat_record_reconst + &
             m_xz_t1
        m_zz_total = (m_zz_t2 - m_zz_t1) * &
             (it*deltat_read_reconst - nt1_record_reconst*deltat_record_reconst)/deltat_record_reconst + &
             m_zz_t1

     else

        trac_f_t1(2,:) = one_time_slice_elastic_1(1,:) 
        m_yx_t1 = one_time_slice_elastic_1(2,:) 
        m_yz_t1 = one_time_slice_elastic_1(3,:) 

        trac_f_t2(2,:) = one_time_slice_elastic_2(1,:) 
        m_yx_t2 = one_time_slice_elastic_2(2,:) 
        m_yz_t2 = one_time_slice_elastic_2(3,:) 

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!linear interpolation in time domain
        !only y component
        trac_f_total(2,:) = (trac_f_t2(2,:) - trac_f_t1(2,:)) * &
             (it*deltat_read_reconst - nt1_record_reconst*deltat_record_reconst)/deltat_record_reconst + &
             trac_f_t1(2,:)

        m_yx_total = (m_yx_t2 - m_yx_t1) * &
             (it*deltat_read_reconst - nt1_record_reconst*deltat_record_reconst)/deltat_record_reconst + &
             m_yx_t1
        m_yz_total = (m_yz_t2 - m_yz_t1) * &
             (it*deltat_read_reconst - nt1_record_reconst*deltat_record_reconst)/deltat_record_reconst + &
             m_yz_t1
        
        
     endif
     

  endif

     
  if( nspec_bd_pnt_acoustic_supply /= 0 .and. p_sv .and. o_rank_acoustic == 0 )then

     !do the data input
     
     if( mod((nt1_record_reconst - record_nt1_reconst), num_step_input) == 0 )then
        call MPI_SIZEOF(Grad_pot_t1(1),size,ierror)
        call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,size,bd_info_type,ierror)

        if( (record_nt2_reconst - nt1_record_reconst ) < num_step_input )then
           count = record_nt2_reconst - nt1_record_reconst
        else

           count = num_step_input

        endif

        offset = (nt1_record_reconst - record_nt1_reconst)*nspec_bd_pnt_acoustic_supply_total*3_8*size
        count_type = 3_8*nspec_bd_pnt_acoustic_supply_total*count
        
        call MPI_FILE_OPEN(MPI_COMM_SELF, './OUTPUT_FILES/reconst_record/acoustic_pnts_data', &
             MPI_MODE_RDONLY, MPI_INFO_NULL, f_num, ierror)

        call MPI_FILE_READ_AT(f_num, offset, data_to_supply_acoustic, count_type,&
             bd_info_type, MPI_STATUS_IGNORE, ierror)

        call MPI_FILE_CLOSE(f_num,ierror)

     endif

     !read two time slices for interpolation
     if( mod((nt1_record_reconst - record_nt1_reconst), num_step_input) /= (num_step_input -1) )then

        slice_1 = mod((nt1_record_reconst - record_nt1_reconst), num_step_input) + 1
     else
        slice_1 = mod((nt1_record_reconst - record_nt1_reconst), num_step_input)
     endif

     slice_2 = slice_1 + 1

     one_time_slice_acoustic_1 = data_to_supply_acoustic(:,:,slice_1)
     one_time_slice_acoustic_2 = data_to_supply_acoustic(:,:,slice_2)

     
     Grad_pot_t1 = one_time_slice_acoustic_1(1,:)
     Pot_x_t1 = one_time_slice_acoustic_1(2,:)
     Pot_z_t1 = one_time_slice_acoustic_1(3,:)


     Grad_pot_t2 = one_time_slice_acoustic_2(1,:)
     Pot_x_t2 = one_time_slice_acoustic_2(2,:)
     Pot_z_t2 = one_time_slice_acoustic_2(3,:)

     !do the time interpolation
     Grad_pot_total = (Grad_pot_t2 - Grad_pot_t1) * &
          (it*deltat_read_reconst - nt1_record_reconst*deltat_record_reconst)/deltat_record_reconst + &
          Grad_pot_t1
     Pot_x_total = (Pot_x_t2 - Pot_x_t1) * &
          (it*deltat_read_reconst - nt1_record_reconst*deltat_record_reconst)/deltat_record_reconst + &
          Pot_x_t1
     Pot_z_total = (Pot_z_t2 - Pot_z_t1) * &
          (it*deltat_read_reconst - nt1_record_reconst*deltat_record_reconst)/deltat_record_reconst + &
          Pot_z_t1


  endif


  !broadcast the result
  if( nspec_bd_pnt_elastic_supply /= 0 ) then

     call MPI_SIZEOF(trac_f_total(1,1),size,ierror)
     call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,size,bd_info_type,ierror)

     if( p_sv )then

        call MPI_BCAST(trac_f_total,3*nspec_bd_pnt_elastic_supply_total,bd_info_type,0,supply_reconst_elastic,ierror)
        call MPI_BCAST(m_xx_total,nspec_bd_pnt_elastic_supply_total,bd_info_type,0,supply_reconst_elastic,ierror)
        call MPI_BCAST(m_xz_total,nspec_bd_pnt_elastic_supply_total,bd_info_type,0,supply_reconst_elastic,ierror)
        call MPI_BCAST(m_zz_total,nspec_bd_pnt_elastic_supply_total,bd_info_type,0,supply_reconst_elastic,ierror)

        trac_f((/1,3/),:) = trac_f_total((/1,3/),booking_reconst_elastic)  
        m_xx   = m_xx_total(booking_reconst_elastic)
        m_xz   = m_xz_total(booking_reconst_elastic)
        m_zz   = m_zz_total(booking_reconst_elastic)
        
     else

        call MPI_BCAST(trac_f_total,3*nspec_bd_pnt_elastic_supply_total,bd_info_type,0,supply_reconst_elastic,ierror)
        call MPI_BCAST(m_yx_total,nspec_bd_pnt_elastic_supply_total,bd_info_type,0,supply_reconst_elastic,ierror)
        call MPI_BCAST(m_yz_total,nspec_bd_pnt_elastic_supply_total,bd_info_type,0,supply_reconst_elastic,ierror)

        trac_f = trac_f_total(:,booking_reconst_elastic)  
        m_yx   = m_yx_total(booking_reconst_elastic)
        m_yz   = m_yz_total(booking_reconst_elastic)

     endif

  endif

  if( nspec_bd_pnt_acoustic_supply /= 0 .and. p_sv )then
     
     call MPI_SIZEOF(Grad_pot_total(1),size,ierror)
     call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,size,bd_info_type,ierror)

     call MPI_BCAST(Grad_pot,nspec_bd_pnt_acoustic_supply_total,bd_info_type,0,supply_reconst_acoustic,ierror)
     call MPI_BCAST(Pot_x_total,nspec_bd_pnt_acoustic_supply_total,bd_info_type,0,supply_reconst_acoustic,ierror)
     call MPI_BCAST(Pot_z_total,nspec_bd_pnt_acoustic_supply_total,bd_info_type,0,supply_reconst_acoustic,ierror)

     Grad_pot  =   Grad_pot_total(booking_reconst_acoustic)
     Pot_x     =   Pot_x_total(booking_reconst_acoustic)
     Pot_z     =   Pot_z_total(booking_reconst_acoustic)

  endif
  
#endif


end subroutine time_interplt_supply_reconst

