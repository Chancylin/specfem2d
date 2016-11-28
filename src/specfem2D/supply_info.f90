!this subroutine is called to import the boundary information needed
!for local simulation. Every time step it will be called
subroutine supply_bd_pnt()

  use specfem_par, only: it,read_nt1,read_nt2,& !original para
                         nspec_bd_pnt_elastic,nspec_bd_pnt_acoustic,&
                         x_final_bd_pnt_elastic,z_final_bd_pnt_elastic,&
                         trac_bd_pnt_elastic,vel_bd_pnt_elastic,&
                         x_final_bd_pnt_acoustic,z_final_bd_pnt_acoustic,&
                         grad_pot_bd_pnt_acoustic,pot_dot_bd_pnt_acoustic,&
                         f_num!,fname
                         

  implicit none
  include "constants.h"

  integer :: i,ios,temp_read                         
  character(len=150) dummystring

  if ( it == 1) then
     nspec_bd_pnt_elastic = 0
     nspec_bd_pnt_acoustic = 0
     !count the total boundary points for 
     open(unit=1,file='./OUTPUT_FILES/bg_record/elastic_pnts_profile',iostat=ios,status='old',action='read')
     do while(ios == 0)
        read(1,"(a)",iostat=ios) dummystring
        if(ios == 0) nspec_bd_pnt_elastic = nspec_bd_pnt_elastic + 1
     enddo
     close(1)

     open(unit=1,file='./OUTPUT_FILES/bg_record/acoustic_pnts_profile',iostat=ios,status='old',action='read')
     do while(ios == 0)
        read(1,"(a)",iostat=ios) dummystring
        if(ios == 0) nspec_bd_pnt_acoustic = nspec_bd_pnt_acoustic + 1
     enddo
     close(1)
 
     !for elastic element
     if( nspec_bd_pnt_elastic /= 0 )then
       
       !allocate the arrays needed at the first time step 
       allocate(x_final_bd_pnt_elastic(nspec_bd_pnt_elastic),z_final_bd_pnt_elastic(nspec_bd_pnt_elastic))
       allocate(trac_bd_pnt_elastic(3,nspec_bd_pnt_elastic))
       allocate(vel_bd_pnt_elastic(3,nspec_bd_pnt_elastic))

       trac_bd_pnt_elastic = 0.0
       vel_bd_pnt_elastic = 0.0

     endif

     !for acoustic element
     if( nspec_bd_pnt_acoustic /= 0 )then
       allocate(x_final_bd_pnt_acoustic(nspec_bd_pnt_acoustic),z_final_bd_pnt_acoustic(nspec_bd_pnt_acoustic))
       allocate(grad_pot_bd_pnt_acoustic(2,nspec_bd_pnt_acoustic))
       allocate(pot_dot_bd_pnt_acoustic(nspec_bd_pnt_acoustic))  

       grad_pot_bd_pnt_acoustic = 0.0
       pot_dot_bd_pnt_acoustic = 0.0
     endif

  endif

  !read the coordinate of interpolation points
  if ( it == 1) then
     f_num = 111

     if( nspec_bd_pnt_elastic /= 0 )then
       open(f_num,file='./OUTPUT_FILES/bg_record/elastic_pnts_profile',iostat=ios,status='old',action='read')
       do i=1,nspec_bd_pnt_elastic
          read(f_num,110) temp_read, x_final_bd_pnt_elastic(i), z_final_bd_pnt_elastic(i)
       enddo
       close(f_num)
     endif

     if( nspec_bd_pnt_acoustic /= 0 )then
       open(f_num,file='./OUTPUT_FILES/bg_record/acoustic_pnts_profile',iostat=ios,status='old',action='read')
       do i=1,nspec_bd_pnt_acoustic
          read(f_num,110) temp_read, x_final_bd_pnt_acoustic(i), z_final_bd_pnt_acoustic(i)
       enddo
       close(f_num)
     endif 

  endif

  110 format(i5,2(es12.4,2x))!consistent with format in 'locate_recording_point.F90'
 
  !read the stored boundary info
  if (it < read_nt1 .or. it > read_nt2 ) return

  !apply the time interpolation 
  call time_interplt_supply() 

  !f_num=113
  !write(fname,"('./OUTPUT_FILES/bg_record/&
  !      &elastic_pnts/nt_',i6.6)")it
  !
  !!formatted reading
  !!open(unit=f_num,file=trim(fname),status='old',&
  !!     action='read',iostat=ios)

  !!unformatted reading
  !open(unit=f_num,file=trim(fname),access='direct',status='old',&
  !     action='read',iostat=ios,recl=length_unf_1)

  !if( ios /= 0 ) stop 'error reading values at profile points' 
  !
  !do i=1,nspec_bd_pnt_elastic
  !   !read(f_num,111) trac_bd_pnt_elastic(:,i),vel_bd_pnt_elastic(:,i)
  !   read(f_num,rec=i) trac_bd_pnt_elastic(:,i),vel_bd_pnt_elastic(:,i)
  !enddo  
 
  !close(f_num)

  !write(fname,"('./OUTPUT_FILES/bg_record/&
  !      &acoustic_pnts/nt_',i6.6)")it
  !!formatted reading
  !!open(unit=f_num,file=trim(fname),status='old',&
  !!     action='read',iostat=ios)

  !!unformatted reading
  !open(unit=f_num,file=trim(fname),access='direct',status='old',&
  !     action='read',iostat=ios,recl=length_unf_2)

  !if( ios /= 0 ) stop 'error reading values at profile points' 
  !
  !do i= 1,nspec_bd_pnt_acoustic
  !   !read(f_num,112) grad_pot_bd_pnt_acoustic(:,i),pot_dot_bd_pnt_acoustic(i)
  !   read(f_num,rec=i) grad_pot_bd_pnt_acoustic(:,i),pot_dot_bd_pnt_acoustic(i)
  !enddo
  !
  !close(f_num)
  !format need to be consistent with format in 'record_bd_pnt.F90'
  !111 format(6(es12.4,2x)) !112 column
  !112 format(3(es12.4,2x)) !36 column

end subroutine supply_bd_pnt

subroutine time_interplt_supply()

  use specfem_par, only: it,deltat_read,&
                         record_nt1,record_nt2, deltat_record,&
                         nspec_bd_pnt_elastic,nspec_bd_pnt_acoustic,&
                         trac_bd_pnt_elastic,vel_bd_pnt_elastic,&
                         grad_pot_bd_pnt_acoustic,pot_dot_bd_pnt_acoustic

  implicit none
  include "constants.h"

  integer :: nt1_record,nt2_record
  integer :: i,f_num_1,f_num_2,ios
  character(len=150) :: fname_1,fname_2
  integer :: length_unf_1
  integer :: length_unf_2

  real(kind=CUSTOM_REAL), dimension(3) :: trac_bd_pnt_t1,vel_bd_pnt_t1,trac_bd_pnt_t2,vel_bd_pnt_t2
  real(kind=CUSTOM_REAL), dimension(2) :: grad_pot_bd_pnt_t1, grad_pot_bd_pnt_t2
  real(kind=CUSTOM_REAL) :: pot_dot_bd_pnt_t1, pot_dot_bd_pnt_t2 
  double precision :: diff_deltat

  nt1_record = floor(it * deltat_read / deltat_record)
  if(nt1_record < record_nt1 ) nt1_record = record_nt1
  if(nt1_record >=  record_nt2 ) nt1_record = record_nt2 - 1
  nt2_record = nt1_record + 1
  !nt2_record = ceiling(it * deltat_read / deltat_record)
  diff_deltat = 0.5 * (deltat_record - deltat_read)


  !!!this is the recording length for unformatted recording
  !inquire (iolength = length_unf_1) trac_bd_pnt_elastic(:,1),vel_bd_pnt_elastic(:,1)
  !inquire (iolength = length_unf_2) grad_pot_bd_pnt_acoustic(:,1),pot_dot_bd_pnt_acoustic(1)
  if( nspec_bd_pnt_elastic /= 0 ) then

    inquire (iolength = length_unf_1) trac_bd_pnt_t1(:),vel_bd_pnt_t1(:)

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


    do i=1,nspec_bd_pnt_elastic
       read(f_num_1,rec=i) trac_bd_pnt_t1(:),vel_bd_pnt_t1(:)
       read(f_num_2,rec=i) trac_bd_pnt_t2(:),vel_bd_pnt_t2(:)

       !!!linear interpolation in time space
       trac_bd_pnt_elastic(:,i) = (trac_bd_pnt_t2(:) - trac_bd_pnt_t1(:)) * &
                                  (it*deltat_read - nt1_record*deltat_record)/deltat_record + &
                                  trac_bd_pnt_t1(:)
       !note that velocity is recorded as prediction at half time step
       vel_bd_pnt_elastic(:,i) = (vel_bd_pnt_t2(:) - vel_bd_pnt_t1(:)) * &
                                 (it*deltat_read - nt1_record*deltat_record + diff_deltat)/deltat_record + &
                                 vel_bd_pnt_t1(:)
       
    enddo  
 
    close(f_num_1)
    close(f_num_2)

  endif  


  !acoustic elements
  if( nspec_bd_pnt_acoustic /= 0 ) then 

     inquire (iolength = length_unf_2) grad_pot_bd_pnt_t1(:),pot_dot_bd_pnt_t1

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


     do i= 1,nspec_bd_pnt_acoustic
        read(f_num_1,rec=i) grad_pot_bd_pnt_t1(:),pot_dot_bd_pnt_t1
        read(f_num_2,rec=i) grad_pot_bd_pnt_t2(:),pot_dot_bd_pnt_t2

        grad_pot_bd_pnt_acoustic(:,i) = (grad_pot_bd_pnt_t2(:) - grad_pot_bd_pnt_t1(:)) * &
                                        (it*deltat_read - nt1_record*deltat_record)/deltat_record + &
                                        grad_pot_bd_pnt_t1(:)

        pot_dot_bd_pnt_acoustic(i) = (pot_dot_bd_pnt_t2 - pot_dot_bd_pnt_t1) * &
                                     (it*deltat_read - nt1_record*deltat_record + diff_deltat)/deltat_record + &
                                     pot_dot_bd_pnt_t1

     enddo
     
     close(f_num_1)
     close(f_num_2)
  endif

  

end subroutine time_interplt_supply


subroutine supply_pnt_reconst()

  !use specfem_par, only: it,read_nt1_reconst,read_nt2 !original para
                         !nspec_bd_pnt_elastic,nspec_bd_pnt_acoustic,&
                         !x_final_bd_pnt_elastic,z_final_bd_pnt_elastic,&
                         !trac_f,&
                         !x_final_bd_pnt_acoustic,z_final_bd_pnt_acoustic
                        !grad_pot_bd_pnt_acoustic,pot_dot_bd_pnt_acoustic,&
                         

  implicit none
  include "constants.h"
 
  !read the stored boundary info
  !if (it < read_nt1_reconst .or. it > read_nt2_reconst ) return

  !apply the time interpolation 
  call time_interplt_supply_reconst() 

end subroutine supply_pnt_reconst



subroutine time_interplt_supply_reconst()

  use specfem_par, only: it,deltat_read_reconst,&
                         record_nt1_reconst,record_nt2_reconst, deltat_record_reconst,&
                         nspec_bd_pnt_elastic,nspec_bd_pnt_acoustic, &
                         trac_f,m_xx,m_xz,m_zz,&
                         Grad_pot,Pot_x,Pot_z

  implicit none
  include "constants.h"

  integer :: nt1_record_reconst,nt2_record_reconst
  integer :: i,f_num_1,f_num_2,ios
  character(len=150) :: fname_1,fname_2
  integer :: length_unf_1
  integer :: length_unf_2

  
  real(kind=CUSTOM_REAL), dimension(3) :: trac_f_t1,trac_f_t2
  real(kind=CUSTOM_REAL) :: m_xx_t1,m_xz_t1,m_zz_t1
  real(kind=CUSTOM_REAL) :: m_xx_t2,m_xz_t2,m_zz_t2
  real(kind=CUSTOM_REAL) :: Grad_pot_t1,Pot_x_t1,Pot_z_t1 
  real(kind=CUSTOM_REAL) :: Grad_pot_t2,Pot_x_t2,Pot_z_t2 
  !real(kind=CUSTOM_REAL), dimension(3) :: trac_bd_pnt_t1,vel_bd_pnt_t1,trac_bd_pnt_t2,vel_bd_pnt_t2
  !real(kind=CUSTOM_REAL), dimension(2) :: grad_pot_bd_pnt_t1, grad_pot_bd_pnt_t2
  !real(kind=CUSTOM_REAL) :: pot_dot_bd_pnt_t1, pot_dot_bd_pnt_t2 

  nt1_record_reconst = floor(it * deltat_read_reconst / deltat_record_reconst)
  if(nt1_record_reconst < record_nt1_reconst ) nt1_record_reconst = record_nt1_reconst
  if(nt1_record_reconst >=  record_nt2_reconst ) nt1_record_reconst = record_nt2_reconst - 1
  nt2_record_reconst = nt1_record_reconst + 1
  !nt2_record_reconst = ceiling(it * deltat_read_reconst / deltat_record_reconst)

  if( nspec_bd_pnt_elastic /= 0 ) then

    !you may change the iolength once you take the moment density tensor into account
    inquire (iolength = length_unf_1) trac_f_t1(:),m_xx_t1,m_xz_t1,m_zz_t1

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


    do i=1,nspec_bd_pnt_elastic
       read(f_num_1,rec=i) trac_f_t1(:),m_xx_t1,m_xz_t1,m_zz_t1
       read(f_num_2,rec=i) trac_f_t2(:),m_xx_t2,m_xz_t2,m_zz_t2

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

  endif  

  if( nspec_bd_pnt_acoustic /= 0 )then

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

    do i=1,nspec_bd_pnt_acoustic

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

