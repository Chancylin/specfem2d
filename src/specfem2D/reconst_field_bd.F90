!this subroutine is to construct the wavefield in global model by
!taking the information provide from local simulation Note: in this
!implementation, in contrast to treating the exciations as virtual
!sources, we consider them as boundary condition, which will be
!supplied exactly to the corresponding GLL points.
!Jan 12, 2017

subroutine compute_add_trac_f_viscoelastic_bd(accel_elastic,it)
       
  use specfem_par, only: p_sv,nglob_elastic,hprime_xx,hprime_zz,&
       ibool,xix,xiz,gammax,gammaz,&
       nspec_bd_pnt_elastic,&
                                !important para
       ispec_bd_elmt_elastic,bd_pnt_i_bg_elastic,bd_pnt_j_bg_elastic,&
       trac_f, m_xx,m_xz,m_zz,m_yx,m_yz,&
       read_nt1_reconst,read_nt2_reconst

  implicit none
  include "constants.h"

  real(kind=CUSTOM_REAL), dimension(3,nglob_elastic) :: accel_elastic
  integer :: it
  integer :: i_f_source,i,k,iglob
  integer :: a,b
  integer :: delta_func
  double precision :: xixd,xizd,gammaxd,gammazd
  double precision :: G_11,G_13,G_31,G_33,G_21,G_23
  double precision :: mid_temp_1,mid_temp_2,mid_temp_3
  !logical :: switch1,switch2

  !integer :: ios,f_num  

  !switch1 = .FALSE.
  !switch2 = .TRUE.

  if (it < read_nt1_reconst .or. it > read_nt2_reconst ) return
  !there must be a subroutine to assign the trac_f for the time step 'it'
  !the time interpolation may be applied
  !call supply_pnt_reconst()

  if( nspec_bd_pnt_elastic /= 0 ) then

     ! !this is used to export the excitation points to check they are
     ! !the same as the recording point (i.e., exact mesh grids shared
     ! !by local and global domain)
     !  if( it == 1)then
     !    f_num = 117
     !    open(unit=f_num,file='./OUTPUT_FILES/elastic_excitation_info',status='new',&
     !         action='write',iostat=ios)
     !    if( ios /= 0 ) stop 'error saving elastic point profile'
     !  endif

      do i_f_source=1,nspec_bd_pnt_elastic
 
         iglob = ibool(bd_pnt_i_bg_elastic(i_f_source), bd_pnt_j_bg_elastic(i_f_source),&
                 ispec_bd_elmt_elastic(i_f_source))

         !if ( switch1 ) then

         if( p_sv ) then ! P-SV calculation

           !supply the traction force to the corresponding GLL points

           accel_elastic(1,iglob) = accel_elastic(1,iglob) - trac_f(1,i_f_source)
           accel_elastic(3,iglob) = accel_elastic(3,iglob) - trac_f(3,i_f_source) 
         else    ! SH (membrane) calculation

           accel_elastic(2,iglob) = accel_elastic(2,iglob) - trac_f(2,i_f_source)

         endif
         !endif !switch for this term

        !if( switch2 ) then

         if( p_sv ) then !P-SV calculation
        !the procedure for moment tensor here should be checked carefully

           !calcualte G_ik at the GLL point, which is just the excitation point
          
           i = bd_pnt_i_bg_elastic(i_f_source)
           k = bd_pnt_j_bg_elastic(i_f_source)

           xixd    = xix(i,k,ispec_bd_elmt_elastic(i_f_source))
           xizd    = xiz(i,k,ispec_bd_elmt_elastic(i_f_source))
           gammaxd = gammax(i,k,ispec_bd_elmt_elastic(i_f_source))
           gammazd = gammaz(i,k,ispec_bd_elmt_elastic(i_f_source))

           G_11 = m_xx(i_f_source)*xixd  + m_xz(i_f_source)*xizd
           G_13 = m_xx(i_f_source)*gammaxd  + m_xz(i_f_source)*gammazd
           G_31 = m_xz(i_f_source)*xixd  + m_zz(i_f_source)*xizd
           G_33 = m_xz(i_f_source)*gammaxd  + m_zz(i_f_source)*gammazd

           do a = 1,NGLLX
              do b = 1,NGLLZ

                 mid_temp_1 =  G_11*hprime_xx(i,a)*delta_func(k,b) &
                      + G_13*hprime_zz(k,b)*delta_func(i,a) 

                 mid_temp_2 =  G_31*hprime_xx(i,a)*delta_func(k,b) & 
                      + G_33*hprime_zz(k,b)*delta_func(i,a) 

                 iglob = ibool(a,b,ispec_bd_elmt_elastic(i_f_source))

                 accel_elastic(1,iglob) = accel_elastic(1,iglob) &
                      + mid_temp_1
                 accel_elastic(3,iglob) = accel_elastic(3,iglob) &
                      + mid_temp_2 
!!!!note: I guess the idea here is: we already calculate the effective discrete moment tensor point sources
!!!!in the local simulation. Here every discrete moment tensor point source concides with one of GLL points,
!!!!thus we can follow Komatitsch & Tromp, 1999 (A11) to code. Alternative approach is to consider is as 
!!!!finite fault plane and then follow (A12)

               enddo
           enddo

         else 
           !stop 'SH case not supported for moment density tensor so far'  
            i = bd_pnt_i_bg_elastic(i_f_source)
            k = bd_pnt_j_bg_elastic(i_f_source)

            xixd    = xix(i,k,ispec_bd_elmt_elastic(i_f_source))
            xizd    = xiz(i,k,ispec_bd_elmt_elastic(i_f_source))
            gammaxd = gammax(i,k,ispec_bd_elmt_elastic(i_f_source))
            gammazd = gammaz(i,k,ispec_bd_elmt_elastic(i_f_source))

            G_21 = m_yx(i_f_source)*xixd + m_yz(i_f_source)*xizd
            G_23 = m_yx(i_f_source)*gammaxd + m_yz(i_f_source)*gammazd

            do a = 1,NGLLX 
               do b= 1,NGLLZ
                  
                  mid_temp_3 =  G_21*hprime_xx(i,a)*delta_func(k,b) &
                       + G_23*hprime_zz(k,b)*delta_func(i,a) 

                  iglob = ibool(a,b,ispec_bd_elmt_elastic(i_f_source))

                  accel_elastic(2,iglob) = accel_elastic(2,iglob) &
                       + mid_temp_3
                  
               enddo
            enddo

         endif

        !endif!whether to add the moment tensor term
         
         ! if( it == 1 )then
         ! write(f_num,113) ispec_bd_elmt_elastic(i_f_source), &
         !                  bd_pnt_i_bg_elastic(i_f_source), bd_pnt_j_bg_elastic(i_f_source)
         ! endif
    enddo
         ! if( it == 1 )then
         !   close(f_num)
         !   !stop 'obtain the excitation location in the mesh'
         ! endif

  endif

  !113 format(i5,2x,i1,2x,i1)

end subroutine compute_add_trac_f_viscoelastic_bd

subroutine compute_add_pot_f_acoustic_bd(potential_dot_dot_acoustic,it)

  use specfem_par, only: p_sv,nglob_acoustic,hprime_xx,hprime_zz,&
                         ibool,xix,xiz,gammax,gammaz, &
                         nspec_bd_pnt_acoustic,&
                         !acoustic para
                         ispec_bd_elmt_acoustic,bd_pnt_i_bg_acoustic,bd_pnt_j_bg_acoustic,&
                         Grad_pot,Pot_x,Pot_z,&
                         rhostore,&
                         read_nt1_reconst,read_nt2_reconst

  implicit none
  include "constants.h"

  real(kind=CUSTOM_REAL), dimension(nglob_acoustic) :: potential_dot_dot_acoustic
  integer :: it
  integer :: i_pot_source,iglob,i,j,iv,ir,a,b
  integer :: delta_func
  double precision :: xixd,xizd,gammaxd,gammazd
  double precision :: B_1,B_2,mid_temp_acoustic


  if (it < read_nt1_reconst .or. it > read_nt2_reconst ) return

  if( nspec_bd_pnt_acoustic /= 0 ) then

    !stop 'acoustic element detected!'
    do i_pot_source = 1, nspec_bd_pnt_acoustic
       
       if( p_sv ) then

       !this is the term due to gradient of potential

         i = bd_pnt_i_bg_acoustic(i_pot_source)
         j = bd_pnt_j_bg_acoustic(i_pot_source)

         iglob = ibool(i,j,ispec_bd_elmt_acoustic(i_pot_source))

         potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                             - Grad_pot(i_pot_source) &
                                             / rhostore(i,j,ispec_bd_elmt_acoustic(i_pot_source))
                                          !!note how rhol is obtained here


       !this term due to potential itself

        !compute the intermediate term B_i
        iv = j
        ir = i

        xixd    = xix(ir,iv,ispec_bd_elmt_acoustic(i_pot_source))
        xizd    = xiz(ir,iv,ispec_bd_elmt_acoustic(i_pot_source))
        gammaxd = gammax(ir,iv,ispec_bd_elmt_acoustic(i_pot_source))
        gammazd = gammaz(ir,iv,ispec_bd_elmt_acoustic(i_pot_source))

        B_1 = Pot_x(i_pot_source)*xixd + Pot_z(i_pot_source)*xizd
        B_2 = Pot_x(i_pot_source)*gammaxd + Pot_z(i_pot_source)*gammazd

        do a = 1,NGLLX 
           do b = 1,NGLLZ

           mid_temp_acoustic =  B_1*hprime_xx(i,a)*delta_func(j,b) & 
                              + B_2*hprime_zz(j,b)*delta_func(i,a)
           iglob = ibool(a,b,ispec_bd_elmt_acoustic(i_pot_source))

           potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                           + mid_temp_acoustic &
                                           / rhostore(i,j,ispec_bd_elmt_acoustic(i_pot_source))

           enddo
         enddo

       !else !basically nothing needs to do for SH case
       ! stop 'not designed for SH case yet'
       endif

    enddo

  endif

end subroutine compute_add_pot_f_acoustic_bd


subroutine locate_virtual_bd_pnt()

  use specfem_par, only: p_sv,nspec_bd_pnt_elastic,ispec_bd_elmt_elastic,&
                         bd_pnt_i_bg_elastic,bd_pnt_j_bg_elastic,&
                         trac_f,m_xx,m_xz,m_zz,m_yx,m_yz, &
                         nspec_bd_pnt_acoustic,ispec_bd_elmt_acoustic,&
                         bd_pnt_i_bg_acoustic,bd_pnt_j_bg_acoustic,&
                         Grad_pot,Pot_x,Pot_z
  implicit none

  include "constants.h"

  ! Local variables
  integer :: i,ios,f_num
  character(len=150) dummystring
  double precision :: temp1_read,temp2_read

  !calculate the total sources number
  nspec_bd_pnt_elastic = 0
  !count the total boundary points for 
  open(unit=1,file='./OUTPUT_FILES/bg_record/elastic_pnts_profile',iostat=ios,status='old',action='read')
  if( ios /= 0 ) stop 'error reading elastic points profile'
  do while(ios == 0)
     read(1,"(a)",iostat=ios) dummystring
     if(ios == 0) nspec_bd_pnt_elastic = nspec_bd_pnt_elastic + 1
  enddo 
  close(1)

  if( nspec_bd_pnt_elastic /= 0 ) then

     allocate(ispec_bd_elmt_elastic(nspec_bd_pnt_elastic))
     allocate(bd_pnt_i_bg_elastic(nspec_bd_pnt_elastic),bd_pnt_j_bg_elastic(nspec_bd_pnt_elastic))

     allocate(trac_f(3,nspec_bd_pnt_elastic))
     trac_f = 0.0
     
     if( p_sv ) then
        allocate(m_xx(nspec_bd_pnt_elastic),m_xz(nspec_bd_pnt_elastic),m_zz(nspec_bd_pnt_elastic))
     else
        allocate(m_yx(nspec_bd_pnt_elastic),m_yz(nspec_bd_pnt_elastic))
     endif
     

     !read the coordinates of the source points. The coordinate will be the key information
     f_num = 111                                                                 

     open(f_num,file='./OUTPUT_FILES/bg_record/elastic_pnts_profile',iostat=ios,status='old',action='read')
     do i=1,nspec_bd_pnt_elastic
        read(f_num,111) ispec_bd_elmt_elastic(i), bd_pnt_i_bg_elastic(i), bd_pnt_j_bg_elastic(i),&
                        temp1_read,temp2_read
     enddo                                                                    
     close(f_num)                                                             

  endif

  nspec_bd_pnt_acoustic = 0
  open(unit=1,file='./OUTPUT_FILES/bg_record/acoustic_pnts_profile',iostat=ios,status='old',action='read')
  if( ios /= 0 ) stop 'error reading acoustic points profile'
  do while(ios == 0)
     read(1,"(a)",iostat=ios) dummystring
     if(ios == 0) nspec_bd_pnt_acoustic = nspec_bd_pnt_acoustic + 1
  enddo

  close(1)

  if( nspec_bd_pnt_acoustic /= 0 ) then

     allocate(ispec_bd_elmt_acoustic(nspec_bd_pnt_acoustic))
     allocate(bd_pnt_i_bg_acoustic(nspec_bd_pnt_acoustic),bd_pnt_j_bg_acoustic(nspec_bd_pnt_acoustic))


     !only P-SV case has excitation in acoustic domain
     if( p_sv ) then
        allocate(Grad_pot(nspec_bd_pnt_acoustic))
        allocate(Pot_x(nspec_bd_pnt_acoustic),Pot_z(nspec_bd_pnt_acoustic))
     endif

     f_num = 111
     open(f_num,file='./OUTPUT_FILES/bg_record/acoustic_pnts_profile',iostat=ios,status='old',action='read')
     do i=1,nspec_bd_pnt_acoustic
        read(f_num,111) ispec_bd_elmt_acoustic(i), bd_pnt_i_bg_acoustic(i), bd_pnt_j_bg_acoustic(i),&
                        temp1_read,temp2_read
     enddo                                                                    
     close(f_num)                                                             

  endif

  !consistent with format in 'locate_recording_point.F90'
  111 format(i5,2x,i1,2x,i1,2x,2(es12.4,2x))

end subroutine locate_virtual_bd_pnt


function delta_func(a,b)

  integer,intent(in) :: a,b
  integer :: delta_func
  if ( a == b )then
       delta_func=1
     else
       delta_func=0
  endif

end function delta_func
