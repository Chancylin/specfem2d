!this subroutine is to construct the wavefield in global model by taking the information provide from local simulation
subroutine reconst_field_solid()
!Oct. 24, 2016
!at the current stage, we suppose the the part golbal model out of the local model
!share the some inteface with the local model (i.e., same mesh around the interface)

!the traction term
!---left absorbing boundary
  if( type_side == 'L' ) then
    i = 1
    do j = 1,NGLLZ
       iglob = ibool(i,j,ispec)

       ! external velocity model
       if( assign_external_model ) then
         cpl = vpext(i,j,ispec)
         csl = vsext(i,j,ispec)
         rhol = rhoext(i,j,ispec)
       endif

       rho_vp = rhol*cpl
       rho_vs = rhol*csl

       xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
       zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
       jacobian1D = sqrt(xgamma**2 + zgamma**2)
       nx = - zgamma / jacobian1D
       nz = + xgamma / jacobian1D

       weight = jacobian1D * wzgll(j)

       call pnt_info_interpl_solid(iglob,veloc_x_store,veloc_y_store,veloc_z_store,&
                             tx_store,ty_store,tz_store)

       vx = veloc_elastic(1,iglob) - veloc_x_store
       vy = veloc_elastic(2,iglob) - veloc_y_store
       vz = veloc_elastic(3,iglob) - veloc_z_store

       vn = nx*vx+nz*vz

       tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
       ty = rho_vs*vy
       tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)
 
       if( (codeabs_corner(1,ispecabs) .and. j == 1) .or. (codeabs_corner(3,ispecabs) .and. j == NGLLZ) ) then
       !apply the averaging to deal with the boundary corner: Left-Bottom and Left-Top
         accel_elastic(1,iglob) = accel_elastic(1,iglob) - 0.5*(tx - tx_store)*weight
         accel_elastic(2,iglob) = accel_elastic(2,iglob) - 0.5*(ty - ty_store)*weight
         accel_elastic(3,iglob) = accel_elastic(3,iglob) - 0.5*(tz - tz_store)*weight
 
         else

       !!confirm that tx_store should have plus sign as the contribution to the RHS of the equation
           accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tx - tx_store)*weight
           accel_elastic(2,iglob) = accel_elastic(2,iglob) - (ty - ty_store)*weight
           accel_elastic(3,iglob) = accel_elastic(3,iglob) - (tz - tz_store)*weight
       endif

    enddo
  endif



!the moment density term

end subroutine reconst_field_solid


subroutine compute_add_trac_f_viscoelastic(accel_elastic,it)
       
  use specfem_par, only: p_sv,elastic,nglob_elastic,&
                         nspec_bd_pnt_elastic,&
                         trac_f,&
                         ispec_selected_elastic_source,&  !this should be located by a subroutine
                         hxis_trac_f_store,hgammas_trac_f_store,ibool
  implicit none
  include "constants.h"
  
  real(kind=CUSTOM_REAL), dimension(3,nglob_elastic) :: accel_elastic
  integer :: it
  integer :: i_f_source,i,j,iglob
  double precision :: hlagrange 

  !there must be a subroutine to assign the trac_f for the time step 'it'
  !the time interpolation may be applied
  if( nspec_bd_pnt_elastic /= 0 ) then


      do i_f_source=1,nspec_bd_pnt_elastic

         if( p_sv ) then ! P-SV calculation

           do j = 1,NGLLZ
             do i = 1,NGLLX
               !there must be step to locate these force source. write another subroutine
               iglob = ibool(i,j,ispec_selected_elastic_source(i_f_source))
               hlagrange = hxis_trac_f_store(i_f_source,i) &
                           * hgammas_trac_f_store(i_f_source,j)
               accel_elastic(1,iglob) = accel_elastic(1,iglob) & 
                                        + trac_f(1,i_f_source)*hlagrange
               accel_elastic(3,iglob) = accel_elastic(3,iglob) &
                                        + trac_f(3,i_f_source)*hlagrange
             enddo
           enddo
         else    ! SH (membrane) calculation
           do j = 1,NGLLZ
             do i = 1,NGLLX
               iglob = ibool(i,j,ispec_selected_elastic_source(i_f_source))
               hlagrange = hxis_trac_f_store(i_f_source,i) &
                           * hgammas_trac_f_store(i_f_source,j)
               accel_elastic(2,iglob) = accel_elastic(2,iglob) &
                                        + trac_f(2,i_f_source)*hlagrange
             enddo
           enddo

         endif

  enddo
  endif
end subroutine compute_add_trac_f_viscoelastic

subroutine setup_trac_f_sources()
 
  use specfem_par, only: nspec_bd_pnt_elastic,nspec_bd_pnt_acoustic, &
                         coord,ibool,nglob,nspec,elastic, &
                         hxis_trac_f,hgammas_trac_f,&
                         hpxis_trac_f,hpgammas_trac_f,hxis_trac_f_store,hgammas_trac_f_store,&
                         x_final_bd_pnt_elastic,z_final_bd_pnt_elastic,ispec_selected_elastic_source, &
                         is_proc_trac_f_source,nb_proc_trac_f_source, & !does this really matter?
                         xigll,zigll,npgeo, &
                         nproc,myrank,xi_trac_f,gamma_trac_f,coorg,knods,ngnod &

  implicit none
 
  include "constants.h"
 
  ! Local variables
  integer :: i,ios,temp_read,temp1_read,temp2_read,f_num
  character(len=150) dummystring
  integer :: f_num 
  integer :: i_f_source,ispec
  integer :: iglob_trac_f_source !local or global variable?


  !calculate the total sources number
  if ( it == 1) then
   nspec_bd_pnt_elastic = 0
   nspec_bd_pnt_acoustic = 0
   !count the total boundary points for 
   open(unit=1,file='./OUTPUT_FILES/reconst_record/elastic_pnts_profile',iostat=ios,status='old',action='read')
   do while(ios == 0)
      read(1,"(a)",iostat=ios) dummystring
      if(ios == 0) nspec_bd_pnt_elastic = nspec_bd_pnt_elastic + 1
   enddo 
   close(1)

   open(unit=1,file='./OUTPUT_FILES/reconst_record/acoustic_pnts_profile',iostat=ios,status='old',action='read')
   do while(ios == 0)
      read(1,"(a)",iostat=ios) dummystring
      if(ios == 0) nspec_bd_pnt_acoustic = nspec_bd_pnt_acoustic + 1
   enddo
   close(1)

   endif

  !read the coordinates of the source points. The coordinate will be the key information
  f_num = 111                                                                 
                                                                              
  if( nspec_bd_pnt_elastic /= 0 ) then
     open(f_num,file='./OUTPUT_FILES/reconst_record/elastic_pnts_profile',iostat=ios,status='old',action='read')
     do i=1,nspec_bd_pnt_elastic
        read(f_num,111) temp_read,temp1_read,temp2_read, x_final_bd_pnt_elastic(i), z_final_bd_pnt_elastic(i)
     enddo                                                                    
     close(f_num)                                                             
  endif
  if( nspec_bd_pnt_acoustic /= 0 )then 
     open(f_num,file='./OUTPUT_FILES/reconst_record/acoustic_pnts_profile',iostat=ios,status='old',action='read')
     do i=1,nspec_bd_pnt_acoustic
        read(f_num,111) temp_read,temp1_read,temp2_read, x_final_bd_pnt_acoustic(i), z_final_bd_pnt_acoustic(i)
     enddo
     close(f_num)
  endif
  
  !consistent with format in 'locate_recording_point.F90'
  111 format(i5,i1,i1,2(es12.4,2x))!

  allocate(hxis_trac_f(NGLLX))
  allocate(hgammas_trac_f(NGLLZ))
  allocate(hpxis_trac_f(NGLLX))
  allocate(hpgammas_trac_f(NGLLZ))
  allocate(hxis_trac_f_store(nspec_bd_pnt_elastic,NGLLX))
  allocate(hgammas_trac_f_store(nspec_bd_pnt_elastic,NGLLZ))

  print *, 'For wavefield reconstruction: here we locate the traction force sources'
  !note: because here we locate the tranction force source according to their coordinates, following the build-in way locating
  !the sources and receviers. Therefore, I think the spatial interploation is not necessary?
 
   do i_f_source=1,nspec_bd_pnt_elastic


     ! collocated force source: here we just take advantange of the available subroutine 'locate_source_force'
     ! the main purpose here is to calculate ispec_selected_elastic_source,xi_trac_f(i_f_source),
     ! gamma_trac_f
     call locate_source_force(ibool,coord,nspec,nglob,xigll,zigll,x_final_bd_pnt_elastic(i_f_source),z_final_bd_pnt_elastic(i_f_source), &
         ispec_selected_elastic_source(i_f_source),is_proc_source(i_f_source),nb_proc_source(i_f_source), &
         nproc,myrank,xi_trac_f(i_f_source),gamma_trac_f(i_f_source),coorg,knods,ngnod,npgeo, &
         iglob_trac_f_source)

   enddo


   ! define and store Lagrange interpolators at all the sources
   !compute hxis_trac_f_store, hgammas_trac_f_store, which will be used in 'compute_add_trac_f_viscoelastic'
   do i_f_source=1,nspec_bd_pnt_elastic

      call lagrange_any(xi_trac_f(i_f_source),NGLLX,xigll,hxis_trac_f,hpxis_trac_f)
      call lagrange_any(gamma_trac_f(i_f_source),NGLLZ,zigll,hgammas_trac_f,hpgammas_trac_f)

      hxis_trac_f_store(i_f_source,:) = hxis_trac_f(:)
      hgammas_trac_f_store(i_f_source,:) = hgammas_trac_f(:)

  enddo

end subroutine setup_trac_f_sources





