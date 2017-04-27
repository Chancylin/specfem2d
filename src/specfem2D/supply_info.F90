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
                         f_num!,fname


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

     print *, 'start broadcasting '
     call MPI_BCAST(nspec_bd_pnt_elastic_supply, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
     call MPI_BCAST(nspec_bd_pnt_acoustic_supply, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

     print *, 'nspec_bd_pnt_elastic_supply = ', nspec_bd_pnt_elastic_supply, ' from rank ', myrank
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

           open(f_num,file='./OUTPUT_FILES/bg_record/elastic_pnts_profile_total',iostat=ios,status='old',action='read')

           do i=1,nspec_bd_pnt_elastic_supply
              read(f_num,110) temp_read,x_final_bd_pnt_elastic_supply(i), z_final_bd_pnt_elastic_supply(i)
           enddo

           close(f_num)

        endif

        if( nspec_bd_pnt_acoustic_supply /= 0 )then

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
                         grad_pot_bd_pnt_acoustic,pot_dot_bd_pnt_acoustic

  implicit none
  include "constants.h"

  integer :: nt1_record,nt2_record
  integer :: i,f_num_1,f_num_2
  character(len=150) :: fname_1,fname_2

  ! real(kind=CUSTOM_REAL), dimension(3) :: trac_bd_pnt_t1,vel_bd_pnt_t1,trac_bd_pnt_t2,vel_bd_pnt_t2
  ! real(kind=CUSTOM_REAL), dimension(2) :: grad_pot_bd_pnt_t1, grad_pot_bd_pnt_t2
  ! real(kind=CUSTOM_REAL) :: pot_dot_bd_pnt_t1, pot_dot_bd_pnt_t2 
  double precision :: diff_deltat
  !define the array for mpi read
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: trac_bd_pnt_t1,vel_bd_pnt_t1,trac_bd_pnt_t2,vel_bd_pnt_t2
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: grad_pot_bd_pnt_t1, grad_pot_bd_pnt_t2
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: pot_dot_bd_pnt_t1, pot_dot_bd_pnt_t2 
#ifdef USE_MPI
  !integer :: size,bd_info_type_elastic,bd_info_type_acoustic,ierror
  integer :: size,bd_info_type,ierror
  ! integer(kind=MPI_OFFSET_KIND) :: offset
#else
  integer :: length_unf_1
  integer :: length_unf_2
  integer :: ios
#endif
  
#ifdef USE_MPI
  if( myrank == 0 )then

     nt1_record = floor(it * deltat_read / deltat_record)
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

        allocate(trac_bd_pnt_t1(3,nspec_bd_pnt_elastic_supply),vel_bd_pnt_t1(3,nspec_bd_pnt_elastic_supply))
        allocate(trac_bd_pnt_t2(3,nspec_bd_pnt_elastic_supply),vel_bd_pnt_t2(3,nspec_bd_pnt_elastic_supply))

        call MPI_SIZEOF(trac_bd_pnt_t1(1,1),size,ierror)
        call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,size,bd_info_type,ierror)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(fname_1,"('./OUTPUT_FILES/bg_record/&
             &elastic_pnts/nt_',i6.6)")nt1_record

        call MPI_FILE_OPEN(MPI_COMM_SELF, fname_1, &
             MPI_MODE_RDONLY, MPI_INFO_NULL, f_num_1, ierror)
        call MPI_FILE_READ(f_num_1, trac_bd_pnt_t1, 3*nspec_bd_pnt_elastic_supply,&
             bd_info_type, MPI_STATUS_IGNORE, ierror)
        ! call MPI_FILE_GET_POSITION(f_num_1, offset, ierror)
        ! print *, '1st reading, offset: ', offset, ' from rank ', myrank

        call MPI_FILE_READ(f_num_1, vel_bd_pnt_t1, 3*nspec_bd_pnt_elastic_supply,&
             bd_info_type, MPI_STATUS_IGNORE, ierror)
        ! call MPI_FILE_GET_POSITION(f_num_1, offset, ierror)
        ! print *, '2nd reading, offset: ', offset, ' from rank ', myrank

        call MPI_FILE_CLOSE(f_num_1,ierror)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(fname_2,"('./OUTPUT_FILES/bg_record/&
             &elastic_pnts/nt_',i6.6)")nt2_record

        call MPI_FILE_OPEN(MPI_COMM_SELF, fname_2, &
             MPI_MODE_RDONLY, MPI_INFO_NULL, f_num_2, ierror)
        call MPI_FILE_READ(f_num_2, trac_bd_pnt_t2, 3*nspec_bd_pnt_elastic_supply,&
             bd_info_type, MPI_STATUS_IGNORE, ierror)
        call MPI_FILE_READ(f_num_2, vel_bd_pnt_t2, 3*nspec_bd_pnt_elastic_supply,&
             bd_info_type, MPI_STATUS_IGNORE, ierror)

        call MPI_FILE_CLOSE(f_num_2,ierror)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i=1,nspec_bd_pnt_elastic_supply
!!!linear interpolation in time space
           trac_bd_pnt_elastic(:,i) = (trac_bd_pnt_t2(:,i) - trac_bd_pnt_t1(:,i)) * &
                (it*deltat_read - nt1_record*deltat_record)/deltat_record + &
                trac_bd_pnt_t1(:,i)
           !note that velocity is recorded as prediction at half time step
           vel_bd_pnt_elastic(:,i) = (vel_bd_pnt_t2(:,i) - vel_bd_pnt_t1(:,i)) * &
                (it*deltat_read - nt1_record*deltat_record + diff_deltat)/deltat_record + &
                vel_bd_pnt_t1(:,i)

        enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        deallocate(trac_bd_pnt_t1,vel_bd_pnt_t1)
        deallocate(trac_bd_pnt_t2,vel_bd_pnt_t2)

     endif

     if( nspec_bd_pnt_acoustic_supply /= 0 ) then

        allocate(grad_pot_bd_pnt_t1(2,nspec_bd_pnt_acoustic_supply),pot_dot_bd_pnt_t1(nspec_bd_pnt_acoustic_supply))
        allocate(grad_pot_bd_pnt_t2(2,nspec_bd_pnt_acoustic_supply),pot_dot_bd_pnt_t2(nspec_bd_pnt_acoustic_supply))

        call MPI_SIZEOF(grad_pot_bd_pnt_t1(1,1),size,ierror)
        call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,size,bd_info_type,ierror)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(fname_1,"('./OUTPUT_FILES/bg_record/&
             &acoustic_pnts/nt_',i6.6)")nt1_record

        call MPI_FILE_OPEN(MPI_COMM_SELF, fname_1, &
             MPI_MODE_RDONLY, MPI_INFO_NULL, f_num_1, ierror)
        call MPI_FILE_READ(f_num_1, grad_pot_bd_pnt_t1, 2*nspec_bd_pnt_acoustic_supply,&
             bd_info_type, MPI_STATUS_IGNORE, ierror)
        
        call MPI_FILE_READ(f_num_1, pot_dot_bd_pnt_t1, nspec_bd_pnt_acoustic_supply,&
             bd_info_type, MPI_STATUS_IGNORE, ierror)

        call MPI_FILE_CLOSE(f_num_1,ierror)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(fname_2,"('./OUTPUT_FILES/bg_record/&
             &acoustic_pnts/nt_',i6.6)")nt2_record

        call MPI_FILE_OPEN(MPI_COMM_SELF, fname_2, &
             MPI_MODE_RDONLY, MPI_INFO_NULL, f_num_2, ierror)
        call MPI_FILE_READ(f_num_2, grad_pot_bd_pnt_t2, 2*nspec_bd_pnt_acoustic_supply,&
             bd_info_type, MPI_STATUS_IGNORE, ierror)
        call MPI_FILE_READ(f_num_2, pot_dot_bd_pnt_t2, nspec_bd_pnt_acoustic_supply,&
             bd_info_type, MPI_STATUS_IGNORE, ierror)

        call MPI_FILE_CLOSE(f_num_2,ierror)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i= 1,nspec_bd_pnt_acoustic_supply

           grad_pot_bd_pnt_acoustic(:,i) = (grad_pot_bd_pnt_t2(:,i) - grad_pot_bd_pnt_t1(:,i)) * &
                (it*deltat_read - nt1_record*deltat_record)/deltat_record + &
                grad_pot_bd_pnt_t1(:,i)

           pot_dot_bd_pnt_acoustic(i) = (pot_dot_bd_pnt_t2(i) - pot_dot_bd_pnt_t1(i)) * &
                (it*deltat_read - nt1_record*deltat_record + diff_deltat)/deltat_record + &
                pot_dot_bd_pnt_t1(i)

        enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        deallocate(grad_pot_bd_pnt_t1,pot_dot_bd_pnt_t1)
        deallocate(grad_pot_bd_pnt_t2,pot_dot_bd_pnt_t2)

     endif


     ! if( nspec_bd_pnt_elastic_supply /= 0 ) then
     !    call MPI_SIZEOF(trac_bd_pnt_elastic(1,1),size,ierror)
     !    print *, 'size = ', size, ' from rank ', myrank
     !    call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,size,bd_info_type_elastic,ierror)
     ! endif

     ! if( nspec_bd_pnt_acoustic_supply /= 0 ) then
     !    call MPI_SIZEOF(grad_pot_bd_pnt_acoustic(1,1),size,ierror)
     !    print *, 'size = ', size, ' from rank ', myrank
     !    call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL,size,bd_info_type_acoustic,ierror)
     ! endif
     
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
  
  !!!this is the recording length for unformatted recording
  !inquire (iolength = length_unf_1) trac_bd_pnt_elastic(:,1),vel_bd_pnt_elastic(:,1)
  !inquire (iolength = length_unf_2) grad_pot_bd_pnt_acoustic(:,1),pot_dot_bd_pnt_acoustic(1)
  if( nspec_bd_pnt_elastic_supply /= 0 ) then

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

  endif  


  !acoustic elements
  if( nspec_bd_pnt_acoustic_supply /= 0 ) then 

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

#endif

end subroutine time_interplt_supply


subroutine supply_pnt_reconst()

  implicit none
  include "constants.h"
 
  !read the stored boundary info
  !if (it < read_nt1_reconst .or. it > read_nt2_reconst ) return

  !apply the time interpolation 
  call time_interplt_supply_reconst() 

end subroutine supply_pnt_reconst



subroutine time_interplt_supply_reconst()

  use specfem_par, only: it,p_sv,deltat_read_reconst,&
                         record_nt1_reconst,record_nt2_reconst, deltat_record_reconst,&
                         nspec_bd_pnt_elastic_supply,nspec_bd_pnt_acoustic_supply, &
                         trac_f,m_xx,m_xz,m_zz,m_yx,m_yz,&
                         Grad_pot,Pot_x,Pot_z

  implicit none
  include "constants.h"

  integer :: nt1_record_reconst,nt2_record_reconst
  integer :: i,f_num_1,f_num_2,ios
  character(len=150) :: fname_1,fname_2
  integer :: length_unf_1
  integer :: length_unf_2

  
  real(kind=CUSTOM_REAL), dimension(3) :: trac_f_t1,trac_f_t2
  real(kind=CUSTOM_REAL) :: m_xx_t1,m_xz_t1,m_zz_t1,m_yx_t1,m_yz_t1
  real(kind=CUSTOM_REAL) :: m_xx_t2,m_xz_t2,m_zz_t2,m_yx_t2,m_yz_t2
  real(kind=CUSTOM_REAL) :: Grad_pot_t1,Pot_x_t1,Pot_z_t1 
  real(kind=CUSTOM_REAL) :: Grad_pot_t2,Pot_x_t2,Pot_z_t2 
  trac_f_t1 = 0.0
  trac_f_t2 = 0.0

  nt1_record_reconst = floor(it * deltat_read_reconst / deltat_record_reconst)
  if(nt1_record_reconst < record_nt1_reconst ) nt1_record_reconst = record_nt1_reconst
  if(nt1_record_reconst >=  record_nt2_reconst ) nt1_record_reconst = record_nt2_reconst - 1
  nt2_record_reconst = nt1_record_reconst + 1
  !nt2_record_reconst = ceiling(it * deltat_read_reconst / deltat_record_reconst)

  if( nspec_bd_pnt_elastic_supply /= 0 ) then

     if( p_sv ) then
        
        !inquire (iolength = length_unf_1) trac_f_t1(:),m_xx_t1,m_xz_t1,m_zz_t1
        inquire (iolength = length_unf_1) trac_f_t1(1),trac_f_t1(3),m_xx_t1,m_xz_t1,m_zz_t1

        !elstic elements
        f_num_1=113
        write(fname_1,"('./OUTPUT_FILES/reconst_record/&
             &elastic_pnts/nt_',i6.6)")nt1_record_reconst


        !unformatted reading
        open(unit=f_num_1,file=trim(fname_1),access='direct',status='old',&
             action='read',iostat=ios,recl=length_unf_1)

        if( ios /= 0 ) stop 'error reading values at profile points' 

        f_num_2=114
        write(fname_2,"('./OUTPUT_FILES/reconst_record/&
             &elastic_pnts/nt_',i6.6)")nt2_record_reconst

        !unformatted reading
        open(unit=f_num_2,file=trim(fname_2),access='direct',status='old',&
             action='read',iostat=ios,recl=length_unf_1)

        if( ios /= 0 ) stop 'error reading values at profile points' 


        do i=1,nspec_bd_pnt_elastic_supply
           read(f_num_1,rec=i) trac_f_t1(1),trac_f_t1(3),m_xx_t1,m_xz_t1,m_zz_t1
           read(f_num_2,rec=i) trac_f_t2(1),trac_f_t2(3),m_xx_t2,m_xz_t2,m_zz_t2

!!!linear interpolation in time domain
           trac_f(:,i) = (trac_f_t2(:) - trac_f_t1(:)) * &
                (it*deltat_read_reconst - nt1_record_reconst*deltat_record_reconst)/deltat_record_reconst + &
                trac_f_t1(:)
           m_xx(i) = (m_xx_t2 - m_xx_t1) * &
                (it*deltat_read_reconst - nt1_record_reconst*deltat_record_reconst)/deltat_record_reconst + &
                m_xx_t1
           m_xz(i) = (m_xz_t2 - m_xz_t1) * &
                (it*deltat_read_reconst - nt1_record_reconst*deltat_record_reconst)/deltat_record_reconst + &
                m_xz_t1
           m_zz(i) = (m_zz_t2 - m_zz_t1) * &
                (it*deltat_read_reconst - nt1_record_reconst*deltat_record_reconst)/deltat_record_reconst + &
                m_zz_t1

        enddo

        close(f_num_1)
        close(f_num_2)

     else !SH case
        !only the y-comp traction is needed
        inquire (iolength = length_unf_1) trac_f_t1(2),m_yx_t1,m_yz_t1

        !elstic elements
        f_num_1=113
        write(fname_1,"('./OUTPUT_FILES/reconst_record/&
             &elastic_pnts/nt_',i6.6)")nt1_record_reconst


        !unformatted reading
        open(unit=f_num_1,file=trim(fname_1),access='direct',status='old',&
             action='read',iostat=ios,recl=length_unf_1)

        if( ios /= 0 ) stop 'error reading values at profile points' 

        f_num_2=114
        write(fname_2,"('./OUTPUT_FILES/reconst_record/&
             &elastic_pnts/nt_',i6.6)")nt2_record_reconst

        !unformatted reading
        open(unit=f_num_2,file=trim(fname_2),access='direct',status='old',&
             action='read',iostat=ios,recl=length_unf_1)

        if( ios /= 0 ) stop 'error reading values at profile points' 


        do i=1,nspec_bd_pnt_elastic_supply
           read(f_num_1,rec=i) trac_f_t1(2),m_yx_t1,m_yz_t1
           read(f_num_2,rec=i) trac_f_t2(2),m_yx_t2,m_yz_t2

!!!linear interpolation in time domain
           trac_f(:,i) = (trac_f_t2(:) - trac_f_t1(:)) * &
                (it*deltat_read_reconst - nt1_record_reconst*deltat_record_reconst)/deltat_record_reconst + &
                trac_f_t1(:)
           m_yx(i) = (m_yx_t2 - m_yx_t1) * &
                (it*deltat_read_reconst - nt1_record_reconst*deltat_record_reconst)/deltat_record_reconst + &
                m_yx_t1
           m_yz(i) = (m_yz_t2 - m_yz_t1) * &
                (it*deltat_read_reconst - nt1_record_reconst*deltat_record_reconst)/deltat_record_reconst + &
                m_yz_t1

        enddo

        close(f_num_1)
        close(f_num_2)
     endif
     

  endif

  if( nspec_bd_pnt_acoustic_supply /= 0 .and. p_sv )then

    inquire (iolength = length_unf_2) Grad_pot_t1,Pot_x_t1,Pot_z_t1 
    
    f_num_1=113
    write(fname_1,"('./OUTPUT_FILES/reconst_record/&
          &acoustic_pnts/nt_',i6.6)")nt1_record_reconst
    

    !unformatted reading
    open(unit=f_num_1,file=trim(fname_1),access='direct',status='old',&
         action='read',iostat=ios,recl=length_unf_2)

    if( ios /= 0 ) stop 'error reading values at profile points' 
    
    f_num_2=114
    write(fname_2,"('./OUTPUT_FILES/reconst_record/&
          &acoustic_pnts/nt_',i6.6)")nt2_record_reconst

    !unformatted reading
    open(unit=f_num_2,file=trim(fname_2),access='direct',status='old',&
         action='read',iostat=ios,recl=length_unf_2)

    if( ios /= 0 ) stop 'error reading values at profile points' 

    do i=1,nspec_bd_pnt_acoustic_supply

       read(f_num_1,rec=i) Grad_pot_t1,Pot_x_t1,Pot_z_t1
       read(f_num_2,rec=i) Grad_pot_t2,Pot_x_t2,Pot_z_t2
       
       !!linear interpolation in time domain
       
       Grad_pot(i) = (Grad_pot_t2 - Grad_pot_t1) * &
                                  (it*deltat_read_reconst - nt1_record_reconst*deltat_record_reconst)/deltat_record_reconst + &
                                  Grad_pot_t1
       Pot_x(i) = (Pot_x_t2 - Pot_x_t1) * &
                                  (it*deltat_read_reconst - nt1_record_reconst*deltat_record_reconst)/deltat_record_reconst + &
                                  Pot_x_t1
       Pot_z(i) = (Pot_z_t2 - Pot_z_t1) * &
                                  (it*deltat_read_reconst - nt1_record_reconst*deltat_record_reconst)/deltat_record_reconst + &
                                  Pot_z_t1

    enddo

    close(f_num_1)
    close(f_num_2)

  endif

end subroutine time_interplt_supply_reconst

