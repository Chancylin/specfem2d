!this subroutine is called to import the boundary information needed
!for local simulation. Every time step it will be called
subroutine supply_bd_pnt()

  use specfem_par, only: it,& !original para
                         nspec_bd_pnt_elastic,nspec_bd_pnt_acoustic,&
                         x_final_bd_pnt_elastic,z_final_bd_pnt_elastic,&
                         trac_bd_pnt_elastic,vel_bd_pnt_elastic,&
                         x_final_bd_pnt_acoustic,z_final_bd_pnt_acoustic,&
                         grad_pot_bd_pnt_acoustic,pot_dot_bd_pnt_acoustic,&
                         f_num,fname
                         

  implicit none
  include "constants.h"

  integer :: i,ios,temp_read                         
  character(len=150) dummystring
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
 
  !allocate the arrays needed at the first time step 
  if ( it == 1) then
     allocate(x_final_bd_pnt_elastic(nspec_bd_pnt_elastic),z_final_bd_pnt_elastic(nspec_bd_pnt_elastic))
     allocate(trac_bd_pnt_elastic(3,nspec_bd_pnt_elastic))
     allocate(vel_bd_pnt_elastic(3,nspec_bd_pnt_elastic))

     allocate(x_final_bd_pnt_acoustic(nspec_bd_pnt_acoustic),z_final_bd_pnt_acoustic(nspec_bd_pnt_acoustic))
     allocate(grad_pot_bd_pnt_acoustic(2,nspec_bd_pnt_acoustic))
     allocate(pot_dot_bd_pnt_acoustic(nspec_bd_pnt_acoustic))  
  endif

  !read the coordinate of interpolation points
  if ( it == 1) then
     f_num = 111
     open(f_num,file='./OUTPUT_FILES/bg_record/elastic_pnts_profile',iostat=ios,status='old',action='read')
     do i=1,nspec_bd_pnt_elastic
        read(f_num,110) temp_read, x_final_bd_pnt_elastic(i), z_final_bd_pnt_elastic(i)
     enddo
     close(f_num)

     
     open(f_num,file='./OUTPUT_FILES/bg_record/acoustic_pnts_profile',iostat=ios,status='old',action='read')
     do i=1,nspec_bd_pnt_acoustic
        read(f_num,110) temp_read, x_final_bd_pnt_acoustic(i), z_final_bd_pnt_acoustic(i)
     enddo
     close(f_num)
  endif

  110 format(i5,2(es12.4,2x))!consistent with format in 'locate_recording_point.F90'
 
  !read the stored boundary info
  f_num=113
  write(fname,"('./OUTPUT_FILES/bg_record/&
        &elastic_pnts/nt_',i6.6)")it
  open(unit=f_num,file=trim(fname),status='old',&
       action='read',iostat=ios)
  if( ios /= 0 ) stop 'error reading values at profile points' 
  
  do i=1,nspec_bd_pnt_elastic
     read(f_num,111) trac_bd_pnt_elastic(:,i),vel_bd_pnt_elastic(:,i)
  enddo  
 
  close(f_num)

  write(fname,"('./OUTPUT_FILES/bg_record/&
        &acoustic_pnts/nt_',i6.6)")it
  open(unit=f_num,file=trim(fname),status='old',&
       action='read',iostat=ios)
  if( ios /= 0 ) stop 'error reading values at profile points' 
  
  do i= 1,nspec_bd_pnt_acoustic
     read(f_num,112) grad_pot_bd_pnt_acoustic(:,i),pot_dot_bd_pnt_acoustic(i)
  enddo
  
  close(f_num)
  !format need to be consistent with format in 'record_bd_pnt.F90'
  111 format(6(es12.4,2x)) !112 column
  112 format(3(es12.4,2x)) !36 column

end subroutine supply_bd_pnt
