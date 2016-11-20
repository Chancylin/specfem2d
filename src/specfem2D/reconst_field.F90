!this subroutine is to construct the wavefield in global model by taking the information provide from local simulation
!Note: the basic idea here is to provide the traction force and moment force as the sources to excite
!the wavefield, i.e., reconstructing.
!Nov. 4, 2016
!Since the excitation terms will be treated as the input sources in the global model, the local mesh could
!be different from the global mesh

subroutine compute_add_trac_f_viscoelastic(accel_elastic,it)
       
  use specfem_par, only: p_sv,nglob_elastic,&
                         nspec_bd_pnt_elastic,&
                         trac_f, m_xx,m_xz,m_zz,&
                         ispec_selected_elastic_source_reconst,&  !this should be located by a subroutine
                         hxis_trac_f_store,hgammas_trac_f_store,ibool,&
                         hpxis_trac_f_store,hpgammas_trac_f_store, &
                         xix,xiz,gammax,gammaz, &
                         read_nt1_reconst,read_nt2_reconst
 
  implicit none
  include "constants.h"
  
  real(kind=CUSTOM_REAL), dimension(3,nglob_elastic) :: accel_elastic
  integer :: it
  integer :: i_f_source,i,j,iglob, a,b,k,iv,ir
  double precision :: hlagrange 
  double precision :: xixd,xizd,gammaxd,gammazd
  double precision, dimension(NGLLX,NGLLZ) :: G_11,G_13,G_31,G_33
  double precision, dimension(2,NGLLX,NGLLZ) :: mid_temp
  logical :: switch
  
  switch = .TRUE.
  if (it < read_nt1_reconst .or. it > read_nt2_reconst ) return
  !there must be a subroutine to assign the trac_f for the time step 'it'
  !the time interpolation may be applied
  call supply_pnt_reconst()

  if( nspec_bd_pnt_elastic /= 0 ) then


      do i_f_source=1,nspec_bd_pnt_elastic

         if( p_sv ) then ! P-SV calculation

           !traction force point sources
           do j = 1,NGLLZ
             do i = 1,NGLLX
               !there must be step to locate these force source. write another subroutine
               iglob = ibool(i,j,ispec_selected_elastic_source_reconst(i_f_source))
               hlagrange = hxis_trac_f_store(i_f_source,i) &
                           * hgammas_trac_f_store(i_f_source,j)
               accel_elastic(1,iglob) = accel_elastic(1,iglob) & 
                                        - trac_f(1,i_f_source)*hlagrange
               accel_elastic(3,iglob) = accel_elastic(3,iglob) &
                                        - trac_f(3,i_f_source)*hlagrange
             enddo
           enddo
         else    ! SH (membrane) calculation
           do j = 1,NGLLZ
             do i = 1,NGLLX
               iglob = ibool(i,j,ispec_selected_elastic_source_reconst(i_f_source))
               hlagrange = hxis_trac_f_store(i_f_source,i) &
                           * hgammas_trac_f_store(i_f_source,j)
               accel_elastic(2,iglob) = accel_elastic(2,iglob) &
                                        - trac_f(2,i_f_source)*hlagrange
             enddo
           enddo

         endif

        if( switch ) then

         if( p_sv ) then !P-SV calculation
           !calcualte G_ik
           do k = 1,NGLLZ
              do i = 1,NGLLX
          
                 xixd    = xix(i,k,ispec_selected_elastic_source_reconst(i_f_source))
                 xizd    = xiz(i,k,ispec_selected_elastic_source_reconst(i_f_source))
                 gammaxd = gammax(i,k,ispec_selected_elastic_source_reconst(i_f_source))
                 gammazd = gammaz(i,k,ispec_selected_elastic_source_reconst(i_f_source))

                 G_11(i,k) = m_xx(i_f_source)*xixd  + m_xz(i_f_source)*xizd
                 G_13(i,k) = m_xx(i_f_source)*gammaxd  + m_xz(i_f_source)*gammazd
                 G_31(i,k) = m_xz(i_f_source)*xixd  + m_zz(i_f_source)*xizd
                 G_33(i,k) = m_xz(i_f_source)*gammaxd  + m_zz(i_f_source)*gammazd
             enddo
         enddo

         do a = 1,NGLLX
            do b = 1,NGLLZ
               mid_temp(:,a,b) = ZERO

               do iv=1,NGLLZ
                  do ir=1,NGLLX

                     mid_temp(1,a,b) = mid_temp(1,a,b) + hxis_trac_f_store(i_f_source,ir) &
                                       *hgammas_trac_f_store(i_f_source,iv) &
                                       *( G_11(ir,iv)*hpxis_trac_f_store(i_f_source,a) &
                                         *hgammas_trac_f_store(i_f_source,b) &
                                        + G_13(ir,iv)*hxis_trac_f_store(i_f_source,a) &
                                         *hpgammas_trac_f_store(i_f_source,b))

                     mid_temp(2,a,b) = mid_temp(2,a,b) + hxis_trac_f_store(i_f_source,ir) &
                                       *hgammas_trac_f_store(i_f_source,iv) &
                                       *( G_31(ir,iv)*hpxis_trac_f_store(i_f_source,a) &
                                         *hgammas_trac_f_store(i_f_source,b) &
                                        + G_33(ir,iv)*hxis_trac_f_store(i_f_source,a) &
                                         *hpgammas_trac_f_store(i_f_source,b))

                  enddo
               enddo
               
            enddo
         enddo

        do i = 1,NGLLX
           do j = 1,NGLLZ
                 !the moment density tensor point sources are coinsiding with the traction point sources
                 iglob = ibool(i,j,ispec_selected_elastic_source_reconst(i_f_source))

                 accel_elastic(1,iglob) = accel_elastic(1,iglob) &
                                          + mid_temp(1,i,j) 
                 accel_elastic(3,iglob) = accel_elastic(3,iglob) &
                                          + mid_temp(2,i,j) 
           enddo
       enddo

         else 
           stop 'SH case not supported for moment density tensor so far'  

         endif

        endif!whether to add the moment tensor term

    enddo

  endif
end subroutine compute_add_trac_f_viscoelastic

subroutine setup_trac_f_sources()
 
  use specfem_par, only: nspec_bd_pnt_elastic,nspec_bd_pnt_acoustic, &
                         coord,ibool,nglob,nspec, &
                         trac_f, m_xx,m_xz,m_zz, &
                         hxis_trac_f,hgammas_trac_f,hxis_trac_f_store,hgammas_trac_f_store,&
                         hpxis_trac_f,hpgammas_trac_f,hpxis_trac_f_store,hpgammas_trac_f_store,&
                         x_final_bd_pnt_elastic,z_final_bd_pnt_elastic,ispec_selected_elastic_source_reconst, &
                         x_final_bd_pnt_acoustic,z_final_bd_pnt_acoustic,&
                         !is_proc_trac_f_source,nb_proc_trac_f_source, & !does this really matter?
                         xigll,zigll,npgeo, &
                         nproc,myrank,coorg,knods,ngnod 

  implicit none
 
  include "constants.h"
 
  ! Local variables
  integer :: i,ios,temp_read,temp1_read,temp2_read,f_num
  character(len=150) dummystring
  integer :: i_f_source
  integer :: iglob_trac_f_source !local or global variable?
  double precision, dimension(:), allocatable :: xi_trac_f,gamma_trac_f
  integer, dimension(:), allocatable :: is_proc_source,nb_proc_source 

  !calculate the total sources number
   nspec_bd_pnt_elastic = 0
   nspec_bd_pnt_acoustic = 0
   !count the total boundary points for 
   open(unit=1,file='./OUTPUT_FILES/reconst_record/elastic_pnts_profile',iostat=ios,status='old',action='read')
   if( ios /= 0 ) stop 'error reading elastic points profile'
   do while(ios == 0)
      read(1,"(a)",iostat=ios) dummystring
      if(ios == 0) nspec_bd_pnt_elastic = nspec_bd_pnt_elastic + 1
   enddo 
   print *,"number of recording points in elastic region: ", nspec_bd_pnt_elastic
   close(1)

   open(unit=1,file='./OUTPUT_FILES/reconst_record/acoustic_pnts_profile',iostat=ios,status='old',action='read')
   if( ios /= 0 ) stop 'error reading acoustic points profile'
   do while(ios == 0)
      read(1,"(a)",iostat=ios) dummystring
      if(ios == 0) nspec_bd_pnt_acoustic = nspec_bd_pnt_acoustic + 1
   enddo
   close(1)

   allocate(x_final_bd_pnt_elastic(nspec_bd_pnt_elastic),z_final_bd_pnt_elastic(nspec_bd_pnt_elastic))
   allocate(x_final_bd_pnt_acoustic(nspec_bd_pnt_acoustic),z_final_bd_pnt_acoustic(nspec_bd_pnt_acoustic))
   allocate(trac_f(3,nspec_bd_pnt_elastic))
   !allocate(m_f(nspec_bd_pnt_elastic,3))
   allocate(m_xx(nspec_bd_pnt_elastic),m_xz(nspec_bd_pnt_elastic),m_zz(nspec_bd_pnt_elastic))

  !read the coordinates of the source points. The coordinate will be the key information
  f_num = 111                                                                 
                                                                              
  if( nspec_bd_pnt_elastic /= 0 ) then
     open(f_num,file='./OUTPUT_FILES/reconst_record/elastic_pnts_profile',iostat=ios,status='old',action='read')
     do i=1,nspec_bd_pnt_elastic
        read(f_num,111) temp_read,temp1_read,temp2_read, x_final_bd_pnt_elastic(i),z_final_bd_pnt_elastic(i)
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
  111 format(i5,2x,i1,2x,i1,2x,2(es12.4,2x))

  allocate(hxis_trac_f(NGLLX))
  allocate(hgammas_trac_f(NGLLZ))
  allocate(hpxis_trac_f(NGLLX))
  allocate(hpgammas_trac_f(NGLLZ))
  allocate(hxis_trac_f_store(nspec_bd_pnt_elastic,NGLLX))
  allocate(hgammas_trac_f_store(nspec_bd_pnt_elastic,NGLLZ))
  allocate(hpxis_trac_f_store(nspec_bd_pnt_elastic,NGLLX))
  allocate(hpgammas_trac_f_store(nspec_bd_pnt_elastic,NGLLZ))
  allocate(ispec_selected_elastic_source_reconst(nspec_bd_pnt_elastic))
  allocate(xi_trac_f(nspec_bd_pnt_elastic),gamma_trac_f(nspec_bd_pnt_elastic))
  allocate(is_proc_source(nspec_bd_pnt_elastic),nb_proc_source(nspec_bd_pnt_elastic))

  print *, 'For wavefield reconstruction: here we locate the traction force sources'
  !note: because here we locate the tranction force source according to their coordinates, following the build-in way locating
  !the sources and receviers. Therefore, I think the spatial interploation is not necessary?
 
   do i_f_source=1,nspec_bd_pnt_elastic


     ! collocated force source: here we just take advantange of the available subroutine 'locate_source_force'
     ! the main purpose here is to calculate ispec_selected_elastic_source_reconst,xi_trac_f(i_f_source),
     ! gamma_trac_f
     call locate_source_force(ibool,coord,nspec,nglob,xigll,zigll,x_final_bd_pnt_elastic(i_f_source),&
          z_final_bd_pnt_elastic(i_f_source),ispec_selected_elastic_source_reconst(i_f_source),&
          is_proc_source(i_f_source),nb_proc_source(i_f_source), &
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
      
      hpxis_trac_f_store(i_f_source,:) = hpxis_trac_f(:)
      hpgammas_trac_f_store(i_f_source,:) = hpgammas_trac_f(:)
  enddo

  deallocate(xi_trac_f,gamma_trac_f)
  deallocate(is_proc_source,nb_proc_source)

end subroutine setup_trac_f_sources


