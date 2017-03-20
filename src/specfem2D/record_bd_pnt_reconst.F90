!!!lcx:this file is following 'record_bd_pnt.F90', recording the information along the local model boundary
!!(i.e., traction T_i and f_i = \partial_k m_{ki}, here m_{kl}=u_i n_j C_{ijkl} is the moment density tensor)


!!this subroutine is to record the stress tensor along the local model boundary. 
!!And It seems we can just keep 'record_bd_elmnt_elastic', but not need a new one.
 subroutine record_bd_elmnt_elastic_reconst_f(ispec,i,j,&
            sigma_xx,sigma_xy,sigma_xz,sigma_zz,sigma_zy)

   use specfem_par, only: ispec_bd_elmt_elastic_i,ispec_bd_elmt_elastic_j,ispec_bd_elmt_elastic,&
                          trac_bd_pnt_elastic_reconst,trac_f,nx_bd_pnt_elastic,nz_bd_pnt_elastic,&
                          it,record_nt1_reconst,record_nt2_reconst,& !control time step for recording
                          side_type_elastic,nspec_bd_pnt_elastic,wzgll,wxgll,&
                          xiz,xix,gammaz,gammax,jacobian

   implicit none
   include "constants.h"
   real(kind=CUSTOM_REAL), intent(in) :: sigma_xx,sigma_xy,sigma_xz,sigma_zz,sigma_zy
   integer, intent(in) :: ispec,i,j
   integer :: ispec_bd_pnt_elastic   
   real(kind=CUSTOM_REAL) :: weight,xxi,zxi,xgamma,zgamma,jacobian1D

   if (it < record_nt1_reconst .or. it > record_nt2_reconst ) return

   ispec_bd_pnt_elastic = 0 

   loop1:do ispec_bd_pnt_elastic = 1, nspec_bd_pnt_elastic

      !locate the corresponding recording point
         if ( ispec_bd_elmt_elastic(ispec_bd_pnt_elastic) == ispec .and. ispec_bd_elmt_elastic_i(ispec_bd_pnt_elastic) == i &
            .and. ispec_bd_elmt_elastic_j(ispec_bd_pnt_elastic) == j ) then
            !trac_x
            trac_bd_pnt_elastic_reconst(1,ispec_bd_pnt_elastic) = nx_bd_pnt_elastic(ispec_bd_pnt_elastic)*sigma_xx + &
                 nz_bd_pnt_elastic(ispec_bd_pnt_elastic)*sigma_xz
            !trac_z
            trac_bd_pnt_elastic_reconst(3,ispec_bd_pnt_elastic) = nx_bd_pnt_elastic(ispec_bd_pnt_elastic)*sigma_xz + &
                 nz_bd_pnt_elastic(ispec_bd_pnt_elastic)*sigma_zz
            !trac_y
            trac_bd_pnt_elastic_reconst(2,ispec_bd_pnt_elastic) = nx_bd_pnt_elastic(ispec_bd_pnt_elastic)*sigma_xy + &
                 nz_bd_pnt_elastic(ispec_bd_pnt_elastic)*sigma_zy

            !multiply the weigth coefficients, depending on different sides
            !1. the sign may need to be changed because now these traction force
            !are like external to the interested region
            !2. may merge the 'right and left' and 'bottom and top' case
            if( side_type_elastic(ispec_bd_pnt_elastic) == 'L' )then  !Left

               xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
               zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
               jacobian1D = sqrt(xgamma**2 + zgamma**2)
               weight = jacobian1D * wzgll(j)

               trac_f(:,ispec_bd_pnt_elastic) = trac_bd_pnt_elastic_reconst(:,ispec_bd_pnt_elastic)*weight

            else if( side_type_elastic(ispec_bd_pnt_elastic) == 'R' ) then !Right

               xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
               zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
               jacobian1D = sqrt(xgamma**2 + zgamma**2)
               weight = jacobian1D * wzgll(j)

               trac_f(:,ispec_bd_pnt_elastic) = trac_bd_pnt_elastic_reconst(:,ispec_bd_pnt_elastic)*weight

            else if( side_type_elastic(ispec_bd_pnt_elastic) == 'B' ) then !Bottom

               xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
               zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
               jacobian1D = sqrt(xxi**2 + zxi**2)
               weight = jacobian1D * wxgll(i)

               trac_f(:,ispec_bd_pnt_elastic) = trac_bd_pnt_elastic_reconst(:,ispec_bd_pnt_elastic)*weight

            else if( side_type_elastic(ispec_bd_pnt_elastic) == 'T' ) then !Top

               xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
               zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
               jacobian1D = sqrt(xxi**2 + zxi**2)
               weight = jacobian1D * wxgll(i)

               trac_f(:,ispec_bd_pnt_elastic) = trac_bd_pnt_elastic_reconst(:,ispec_bd_pnt_elastic)*weight

            else

               stop 'type of side is unknown'

            endif

           exit loop1
        
        endif !finding corresponding recording point

   enddo loop1 

 end subroutine record_bd_elmnt_elastic_reconst_f

!a new subroutine to obtain the moment tensor in the local simulation, 
!by calculate the integral of moment density tensor, and save the discrete terms at each GLL point
 subroutine record_bd_elmnt_elastic_reconst_m(ispec,i,j,displ_elastic,&
                           lambdaplus2mu_unrelaxed_elastic,lambdal_unrelaxed_elastic,mul_unrelaxed_elastic)

   use specfem_par, only: p_sv,ispec_bd_elmt_elastic_i,ispec_bd_elmt_elastic_j,ispec_bd_elmt_elastic,&
                          m_xx,m_xz,m_zz,m_zx,&
                          m_xx_reconst,m_xz_reconst,m_zz_reconst,m_zx_reconst,&
                          m_yx,m_yz,m_yx_reconst,m_yz_reconst,&
                          nx_bd_pnt_elastic,nz_bd_pnt_elastic,&
                          it,record_nt1_reconst,record_nt2_reconst,& !control time step for recording
                          side_type_elastic,nspec_bd_pnt_elastic,wzgll,wxgll,&
                          xiz,xix,gammaz,gammax,jacobian

   implicit none
   include "constants.h"

   integer, intent(in) :: ispec,i,j
   real(kind=CUSTOM_REAL), dimension(3), intent(in) :: displ_elastic 
   !!!do we need the external velocity model to calculate he lame paramteters?
   real(kind=CUSTOM_REAL), intent(in) :: mul_unrelaxed_elastic,lambdal_unrelaxed_elastic, &
                             lambdaplus2mu_unrelaxed_elastic
   real(kind=CUSTOM_REAL) :: nx,nz,jacobian1D,xxi,zxi,xgamma,zgamma,weight
   integer :: ispec_bd_pnt_elastic


   if (it < record_nt1_reconst .or. it > record_nt2_reconst ) return
   
   !print *,'ispec = ', ispec
   !print *,lambdaplus2mu_unrelaxed_elastic,lambdal_unrelaxed_elastic,mul_unrelaxed_elastic
   !stop
   ispec_bd_pnt_elastic = 0

   loop1:do ispec_bd_pnt_elastic = 1, nspec_bd_pnt_elastic
         
         !locate the corresponding recording point
         if ( ispec_bd_elmt_elastic(ispec_bd_pnt_elastic) == ispec .and. ispec_bd_elmt_elastic_i(ispec_bd_pnt_elastic) == i &
            .and. ispec_bd_elmt_elastic_j(ispec_bd_pnt_elastic) == j ) then

            nx = nx_bd_pnt_elastic(ispec_bd_pnt_elastic)
            nz = nz_bd_pnt_elastic(ispec_bd_pnt_elastic)

            if( side_type_elastic(ispec_bd_pnt_elastic) == 'L' ) then !left

               xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
               zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
               jacobian1D = sqrt(xgamma**2 + zgamma**2)
               weight = jacobian1D * wzgll(j)

            else if( side_type_elastic(ispec_bd_pnt_elastic) == 'R' ) then !right
                 xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
                 zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
                 jacobian1D = sqrt(xgamma**2 + zgamma**2)
                 weight = jacobian1D * wzgll(j)

            else if( side_type_elastic(ispec_bd_pnt_elastic) == 'B' ) then ! bottom
                 xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                 zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
                 jacobian1D = sqrt(xxi**2 + zxi**2)
                 weight = jacobian1D * wxgll(i)

            else if( side_type_elastic(ispec_bd_pnt_elastic) == 'T' ) then !top

                 xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                 zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
                 jacobian1D = sqrt(xxi**2 + zxi**2)
                 weight = jacobian1D * wxgll(i)

            else

              stop 'type of side is unknown' 
            endif

            !m_kl = u_i * n_j * C_ijkl
            !here we simply use the Lame coefficients for C_ijkl, because we deal with the elastic case.
            !If full anisotropic material is dealed with, then need to use a full tensor C_ijkl
            if( p_sv ) then
               m_xx_reconst(ispec_bd_pnt_elastic) = displ_elastic(1)*nx &
                    *lambdaplus2mu_unrelaxed_elastic + displ_elastic(3)*nz*lambdal_unrelaxed_elastic
               m_xz_reconst(ispec_bd_pnt_elastic) = displ_elastic(1)*nz &
                    *mul_unrelaxed_elastic + displ_elastic(3)*nx* mul_unrelaxed_elastic 
               m_zz_reconst(ispec_bd_pnt_elastic) = displ_elastic(1)*nx &
                    *lambdal_unrelaxed_elastic + displ_elastic(3)*nz*lambdaplus2mu_unrelaxed_elastic 
               m_zx_reconst(ispec_bd_pnt_elastic) = m_xz_reconst(ispec_bd_pnt_elastic) 

               !calculate the linear integral to obtain the discrete moment tensor at each GLL points
               
               m_xx(ispec_bd_pnt_elastic) = m_xx_reconst(ispec_bd_pnt_elastic)*weight
               m_xz(ispec_bd_pnt_elastic) = m_xz_reconst(ispec_bd_pnt_elastic)*weight
               m_zz(ispec_bd_pnt_elastic) = m_zz_reconst(ispec_bd_pnt_elastic)*weight
               m_zx(ispec_bd_pnt_elastic) = m_xz(ispec_bd_pnt_elastic) 

            else
               !calculate m_yx,m_yy,m_yz for SH case
               m_yx_reconst(ispec_bd_pnt_elastic) = displ_elastic(2)*nx &
                    *mul_unrelaxed_elastic
               m_yz_reconst(ispec_bd_pnt_elastic) = displ_elastic(2)*nz &
                    *mul_unrelaxed_elastic

               !calculate the linear integral
               m_yx(ispec_bd_pnt_elastic) = m_yx_reconst(ispec_bd_pnt_elastic)*weight
               m_yz(ispec_bd_pnt_elastic) = m_yz_reconst(ispec_bd_pnt_elastic)*weight

            endif
            
            exit loop1

       endif !locate the corresponding recording point
           
   enddo loop1  
     

 end subroutine record_bd_elmnt_elastic_reconst_m

 subroutine record_bd_elmnt_acoustic_reconst_Grad_pot(ispec,i,j,&
            dux_dxl,dux_dzl)

   use specfem_par, only: p_sv,ispec_bd_elmt_acoustic_i,ispec_bd_elmt_acoustic_j,ispec_bd_elmt_acoustic,&
                          grad_pot_x_reconst,grad_pot_z_reconst,Grad_pot,&
                          nx_bd_pnt_acoustic,nz_bd_pnt_acoustic, &
                          it,record_nt1_reconst,record_nt2_reconst,&
                          side_type_acoustic,nspec_bd_pnt_acoustic,wzgll,wxgll,&
                          xiz,xix,gammaz,gammax,jacobian

   implicit none
   include "constants.h"


   real(kind=CUSTOM_REAL), intent(in) :: dux_dxl,dux_dzl
   integer, intent(in) :: ispec,i,j
   integer :: ispec_bd_pnt_acoustic   
   real(kind=CUSTOM_REAL) :: weight,xxi,zxi,xgamma,zgamma,jacobian1D
   
   if (it < record_nt1_reconst .or. it > record_nt2_reconst ) return

   if ( .not. p_sv ) return !not any recording needed if SH case
   
  loop1:do ispec_bd_pnt_acoustic = 1, nspec_bd_pnt_acoustic

      !locate the corresponding recording point
       if ( ispec_bd_elmt_acoustic(ispec_bd_pnt_acoustic) == ispec .and. ispec_bd_elmt_acoustic_i(ispec_bd_pnt_acoustic) == i &
            .and. ispec_bd_elmt_acoustic_j(ispec_bd_pnt_acoustic) == j ) then
       
         grad_pot_x_reconst(ispec_bd_pnt_acoustic) = dux_dxl
         grad_pot_z_reconst(ispec_bd_pnt_acoustic) = dux_dzl

        ! Grad_pot(ispec_bd_pnt_acoustic) = grad_pot_x_reconst(ispec_bd_pnt_acoustic)&
        !         *nx_bd_pnt_acoustic(ispec_bd_pnt_acoustic)&
        !         + grad_pot_z_reconst(ispec_bd_pnt_acoustic) &
        !         *nz_bd_pnt_acoustic(ispec_bd_pnt_acoustic)

         Grad_pot(ispec_bd_pnt_acoustic) = dux_dxl &
                 *nx_bd_pnt_acoustic(ispec_bd_pnt_acoustic)&
                 + dux_dzl&
                 *nz_bd_pnt_acoustic(ispec_bd_pnt_acoustic)

         if( side_type_acoustic(ispec_bd_pnt_acoustic) == 'L' )then

             xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
             zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
             jacobian1D = sqrt(xgamma**2 + zgamma**2)
             weight = jacobian1D * wzgll(j)

             Grad_pot(ispec_bd_pnt_acoustic) = Grad_pot(ispec_bd_pnt_acoustic)*weight

         else if( side_type_acoustic(ispec_bd_pnt_acoustic) == 'R' ) then
             xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
             zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
             jacobian1D = sqrt(xgamma**2 + zgamma**2)
             weight = jacobian1D * wzgll(j)

             Grad_pot(ispec_bd_pnt_acoustic) = Grad_pot(ispec_bd_pnt_acoustic)*weight

         else if( side_type_acoustic(ispec_bd_pnt_acoustic) == 'B' ) then
             xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
             zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
             jacobian1D = sqrt(xxi**2 + zxi**2)
             weight = jacobian1D * wxgll(i)

             Grad_pot(ispec_bd_pnt_acoustic) = Grad_pot(ispec_bd_pnt_acoustic)*weight

         else if( side_type_acoustic(ispec_bd_pnt_acoustic) == 'T' ) then
             xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
             zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
             jacobian1D = sqrt(xxi**2 + zxi**2)
             weight = jacobian1D * wxgll(i)

             Grad_pot(ispec_bd_pnt_acoustic) = Grad_pot(ispec_bd_pnt_acoustic)*weight

         else

             stop 'type of side is unknown'
         endif
         exit loop1

       endif

   enddo loop1

 end subroutine record_bd_elmnt_acoustic_reconst_Grad_pot

 subroutine record_bd_elmnt_acoustic_reconst_Pot(ispec,i,j,potential_acoustic)
 
   use specfem_par, only: p_sv,ispec_bd_elmt_acoustic_i,ispec_bd_elmt_acoustic_j,ispec_bd_elmt_acoustic,&
                          Pot_x,Pot_z,&
                          nx_bd_pnt_acoustic,nz_bd_pnt_acoustic, &
                          it,record_nt1_reconst,record_nt2_reconst,&
                          side_type_acoustic,nspec_bd_pnt_acoustic,wzgll,wxgll,&
                          xiz,xix,gammaz,gammax,jacobian

   implicit none
   include "constants.h"


   real(kind=CUSTOM_REAL), intent(in) :: potential_acoustic
   real(kind=CUSTOM_REAL) :: pot_x_temp,pot_z_temp
   integer, intent(in) :: ispec,i,j
   integer :: ispec_bd_pnt_acoustic   
   real(kind=CUSTOM_REAL) :: weight,xxi,zxi,xgamma,zgamma,jacobian1D

   if (it < record_nt1_reconst .or. it > record_nt2_reconst ) return
   if ( .not. p_sv ) return !not any recording needed if SH case

   loop1:do ispec_bd_pnt_acoustic = 1, nspec_bd_pnt_acoustic

      if ( ispec_bd_elmt_acoustic(ispec_bd_pnt_acoustic) == ispec .and. ispec_bd_elmt_acoustic_i(ispec_bd_pnt_acoustic) == i &
            .and. ispec_bd_elmt_acoustic_j(ispec_bd_pnt_acoustic) == j ) then

         pot_x_temp = potential_acoustic*nx_bd_pnt_acoustic(ispec_bd_pnt_acoustic)

         pot_z_temp = potential_acoustic*nz_bd_pnt_acoustic(ispec_bd_pnt_acoustic)   
      
         if( side_type_acoustic(ispec_bd_pnt_acoustic) == 'L' )then

             xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
             zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
             jacobian1D = sqrt(xgamma**2 + zgamma**2)
             weight = jacobian1D * wzgll(j)

             Pot_x(ispec_bd_pnt_acoustic) = pot_x_temp*weight
             Pot_z(ispec_bd_pnt_acoustic) = pot_z_temp*weight

         else if( side_type_acoustic(ispec_bd_pnt_acoustic) == 'R' ) then
             xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
             zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
             jacobian1D = sqrt(xgamma**2 + zgamma**2)
             weight = jacobian1D * wzgll(j)

             Pot_x(ispec_bd_pnt_acoustic) = pot_x_temp*weight
             Pot_z(ispec_bd_pnt_acoustic) = pot_z_temp*weight

         else if( side_type_acoustic(ispec_bd_pnt_acoustic) == 'B' ) then
             xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
             zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
             jacobian1D = sqrt(xxi**2 + zxi**2)
             weight = jacobian1D * wxgll(i)

             Pot_x(ispec_bd_pnt_acoustic) = pot_x_temp*weight
             Pot_z(ispec_bd_pnt_acoustic) = pot_z_temp*weight

         else if( side_type_acoustic(ispec_bd_pnt_acoustic) == 'T' ) then
             xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
             zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
             jacobian1D = sqrt(xxi**2 + zxi**2)
             weight = jacobian1D * wxgll(i)

             Pot_x(ispec_bd_pnt_acoustic) = pot_x_temp*weight
             Pot_z(ispec_bd_pnt_acoustic) = pot_z_temp*weight

         else

             stop 'type of side is unknown'
         endif
         exit loop1
     endif 

   enddo loop1
 end subroutine record_bd_elmnt_acoustic_reconst_Pot


 subroutine write_bd_pnts_reconst()

  use specfem_par, only: it,p_sv,& !original para
                         fname,f_num,&
                         nspec_bd_pnt_elastic,nspec_bd_pnt_acoustic,&
                         trac_f,m_xx,m_xz,m_zz,m_yx,m_yz, &
                         Grad_pot,Pot_x,Pot_z, &
                         record_nt1_reconst,record_nt2_reconst !control time step for recording
 
 
  implicit none
  include "constants.h"

  integer :: k
  integer :: ios
  integer :: length_unf_1
  integer :: length_unf_2

  if (it < record_nt1_reconst .or. it > record_nt2_reconst ) return


  !for elastic 
  if( nspec_bd_pnt_elastic /= 0 )then

    f_num=113
    write(fname,"('./OUTPUT_FILES/reconst_record/&
          &elastic_pnts/nt_',i6.6)")it

    if( p_sv ) then
       
       !!!this is the recording length for unformatted recording
       inquire (iolength = length_unf_1) trac_f(1,1),trac_f(1,1),m_xx(1),m_xx(1),m_xx(1)
       !unformatted recording
       open(unit=f_num,file=trim(fname),access='direct',status='new',&
            action='write',iostat=ios,recl=length_unf_1) 
       if( ios /= 0 ) stop 'error saving values at recording points'

       do k = 1, nspec_bd_pnt_elastic
          !save traction_x, traction_z, m_xx, m_xz, m_zz
          write(f_num,rec=k) trac_f(1,k),trac_f(3,k),m_xx(k),m_xz(k),m_zz(k)
       enddo
       
    else
       
       inquire (iolength = length_unf_1) trac_f(1,1),m_yx(1),m_yz(1)
       !unformatted recording
       open(unit=f_num,file=trim(fname),access='direct',status='new',&
            action='write',iostat=ios,recl=length_unf_1) 
       if( ios /= 0 ) stop 'error saving values at recording points'

       do k = 1, nspec_bd_pnt_elastic
          !save traction_y, m_yx, m_yz
          write(f_num,rec=k) trac_f(2,k),m_yx(k),m_yz(k)
       enddo
       
    endif

    close(f_num)
  
  endif 

  !for acoustic
  if( nspec_bd_pnt_acoustic /= 0 )then

     if( p_sv ) then !only need to do the recording if P-SV case

        inquire (iolength = length_unf_2) Grad_pot(1),Pot_x(1),Pot_z(1)

        f_num=114
        write(fname,"('./OUTPUT_FILES/reconst_record/&
             &acoustic_pnts/nt_',i6.6)")it

        !unformatted recording
        open(unit=f_num,file=trim(fname),access='direct',status='new',&
             action='write',iostat=ios,recl=length_unf_2)

        if( ios /= 0 ) stop 'error saving values at recording points'

        do k = 1, nspec_bd_pnt_acoustic
           write(f_num,rec=k) Grad_pot(k),Pot_x(k),Pot_z(k)
        enddo

        close(f_num)

     endif

  endif

 end subroutine write_bd_pnts_reconst
