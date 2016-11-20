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
      !1. the sign may need to be changed because now these traction force are like external to the interested region
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

   use specfem_par, only: ispec_bd_elmt_elastic_i,ispec_bd_elmt_elastic_j,ispec_bd_elmt_elastic,&
                          m_xx,m_xz,m_zz,m_zx,&
                          m_xx_reconst,m_xz_reconst,m_zz_reconst,m_zx_reconst,&
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
            !m_kl = u_i * n_j * C_ijkl
            m_xx_reconst(ispec_bd_pnt_elastic) = displ_elastic(1)*nx &
                 *lambdaplus2mu_unrelaxed_elastic + displ_elastic(3)*nz*lambdal_unrelaxed_elastic
            m_xz_reconst(ispec_bd_pnt_elastic) = displ_elastic(1)*nz &
                 *mul_unrelaxed_elastic + displ_elastic(3)*nx* mul_unrelaxed_elastic 
            m_zz_reconst(ispec_bd_pnt_elastic) = displ_elastic(1)*nx &
                 *lambdal_unrelaxed_elastic + displ_elastic(3)*nz*lambdaplus2mu_unrelaxed_elastic 
            m_zx_reconst(ispec_bd_pnt_elastic) = m_xz_reconst(ispec_bd_pnt_elastic) 

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
            !calculate the linear integral to obtain the discrete moment tensor at each GLL points

            m_xx(ispec_bd_pnt_elastic) = m_xx_reconst(ispec_bd_pnt_elastic)*weight
            m_xz(ispec_bd_pnt_elastic) = m_xz_reconst(ispec_bd_pnt_elastic)*weight
            m_zz(ispec_bd_pnt_elastic) = m_zz_reconst(ispec_bd_pnt_elastic)*weight
            m_zx(ispec_bd_pnt_elastic) = m_xz(ispec_bd_pnt_elastic) 
            exit loop1

       endif !locate the corresponding recording point
           
   enddo loop1  
     

 end subroutine record_bd_elmnt_elastic_reconst_m
!------------------------------------------------------------------------
!note: these two subroutines together are to compute the moment density tensor integral term in SEM (weak form)
!but the implementation could be wrong
! subroutine record_bd_elmnt_elastic_reconst_m(ispec,ispecabs,displ_elastic,&
!                           lambdaplus2mu_unrelaxed_elastic,lambdal_unrelaxed_elastic,mul_unrelaxed_elastic)
!
!   use specfem_par, only: it,& !original para
!                          record_nt1_reconst,record_nt2_reconst,& !control time step for recording
!                          ispec_bd_elmt_elastic_pure_edge,ibool,&
!                          nspec_bd_elmt_elastic_pure_edge,codeabs,&
!                          xix,xiz,gammax,gammaz,jacobian,&
!                          m_xx,m_xz,m_zz,m_zx,nglob
!
!   implicit none
!   include "constants.h"
!
!   integer, intent(in) :: ispec,ispecabs
!   real(kind=CUSTOM_REAL), dimension(3,nglob) :: displ_elastic 
!   !!!do we need the external velocity model to calculate he lame paramteters?
!   real(kind=CUSTOM_REAL), intent(in) :: mul_unrelaxed_elastic,lambdal_unrelaxed_elastic, &
!                             lambdaplus2mu_unrelaxed_elastic
!   real(kind=CUSTOM_REAL) :: nx,nz,jacobian1D,xxi,zxi,xgamma,zgamma
!   integer :: ispec_bd_edge_elastic
!
!   integer :: i,j,iglob
!
!   if (it < record_nt1_reconst .or. it > record_nt2_reconst ) return
!   
!   !print *,'ispec = ', ispec
!   !print *,lambdaplus2mu_unrelaxed_elastic,lambdal_unrelaxed_elastic,mul_unrelaxed_elastic
!   !stop
!   loop1:do ispec_bd_edge_elastic = 1, nspec_bd_elmt_elastic_pure_edge
!         
!         !!!!locate which element in ispec_bd_elmt_elastic_pure_edge(:) will be recorded
!         !!!!note that here we use codeabs() to check the typeside
!         if ( ispec_bd_elmt_elastic_pure_edge(ispec_bd_edge_elastic) == ispec ) then
!            !!left absorbing boundary
!            if( codeabs(IEDGE4,ispecabs) ) then
!              i = 1
!              do j = 1,NGLLZ
!                 iglob = ibool(i,j,ispec)
!                 xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
!                 zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
!                 jacobian1D = sqrt(xgamma**2 + zgamma**2)
!                 nx = - zgamma / jacobian1D
!                 nz = + xgamma / jacobian1D
!
!                 m_xx(ispec_bd_edge_elastic,j) = displ_elastic(1,iglob)*nx &
!                      *lambdaplus2mu_unrelaxed_elastic + displ_elastic(3,iglob)*nz*lambdal_unrelaxed_elastic
!                 m_xz(ispec_bd_edge_elastic,j) = displ_elastic(1,iglob)*nz &
!                      *mul_unrelaxed_elastic + displ_elastic(3,iglob)*nx* mul_unrelaxed_elastic 
!                 m_zz(ispec_bd_edge_elastic,j) = displ_elastic(1,iglob)*nx &
!                      *lambdal_unrelaxed_elastic + displ_elastic(3,iglob)*nz*lambdaplus2mu_unrelaxed_elastic 
!                 m_zx(ispec_bd_edge_elastic,j) = m_xz(ispec_bd_edge_elastic,j) 
!              enddo
!            endif
!            !!right absorbing boundary
!            if( codeabs(IEDGE2,ispecabs) ) then
!              i = NGLLX
!              do j = 1,NGLLZ
!                 iglob = ibool(i,j,ispec)
!                 xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
!                 zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
!                 jacobian1D = sqrt(xgamma**2 + zgamma**2)
!                 nx = + zgamma / jacobian1D
!                 nz = - xgamma / jacobian1D
!
!                 m_xx(ispec_bd_edge_elastic,j) = displ_elastic(1,iglob)*nx &
!                      *lambdaplus2mu_unrelaxed_elastic + displ_elastic(3,iglob)*nz*lambdal_unrelaxed_elastic
!                 m_xz(ispec_bd_edge_elastic,j) = displ_elastic(1,iglob)*nz &
!                      *mul_unrelaxed_elastic + displ_elastic(3,iglob)*nx* mul_unrelaxed_elastic 
!                 m_zz(ispec_bd_edge_elastic,j) = displ_elastic(1,iglob)*nx &
!                      *lambdal_unrelaxed_elastic + displ_elastic(3,iglob)*nz*lambdaplus2mu_unrelaxed_elastic 
!                 m_zx(ispec_bd_edge_elastic,j) = m_xz(ispec_bd_edge_elastic,j) 
!              enddo
!            endif
!            !!bottom absorbing boundary
!            if( codeabs(IEDGE1,ispecabs) ) then
!              j = 1
!              do i = 1,NGLLX
!                 iglob = ibool(i,j,ispec)
!                 xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
!                 zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
!                 jacobian1D = sqrt(xxi**2 + zxi**2)
!                 nx = + zxi / jacobian1D
!                 nz = - xxi / jacobian1D
!
!                 m_xx(ispec_bd_edge_elastic,i) = displ_elastic(1,iglob)*nx &
!                      *lambdaplus2mu_unrelaxed_elastic + displ_elastic(3,iglob)*nz*lambdal_unrelaxed_elastic
!                 m_xz(ispec_bd_edge_elastic,i) = displ_elastic(1,iglob)*nz &
!                      *mul_unrelaxed_elastic + displ_elastic(3,iglob)*nx* mul_unrelaxed_elastic 
!                 m_zz(ispec_bd_edge_elastic,i) = displ_elastic(1,iglob)*nx &
!                      *lambdal_unrelaxed_elastic + displ_elastic(3,iglob)*nz*lambdaplus2mu_unrelaxed_elastic 
!                 m_zx(ispec_bd_edge_elastic,i) = m_xz(ispec_bd_edge_elastic,i) 
!              enddo
!            endif
!            !!top absorbing boundary
!            if( codeabs(IEDGE3,ispecabs) ) then
!              j = NGLLZ
!              do i = 1,NGLLX
!                 iglob = ibool(i,j,ispec)
!                 xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
!                 zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
!                 jacobian1D = sqrt(xxi**2 + zxi**2)
!                 nx = - zxi / jacobian1D
!                 nz = + xxi / jacobian1D
!
!                 m_xx(ispec_bd_edge_elastic,i) = displ_elastic(1,iglob)*nx &
!                      *lambdaplus2mu_unrelaxed_elastic + displ_elastic(3,iglob)*nz*lambdal_unrelaxed_elastic
!                 m_xz(ispec_bd_edge_elastic,i) = displ_elastic(1,iglob)*nz &
!                      *mul_unrelaxed_elastic + displ_elastic(3,iglob)*nx* mul_unrelaxed_elastic 
!                 m_zz(ispec_bd_edge_elastic,i) = displ_elastic(1,iglob)*nx &
!                      *lambdal_unrelaxed_elastic + displ_elastic(3,iglob)*nz*lambdaplus2mu_unrelaxed_elastic 
!                 m_zx(ispec_bd_edge_elastic,i) = m_xz(ispec_bd_edge_elastic,i) 
!              enddo
!            endif
!        ! ! get unrelaxed elastic parameters of current spectral element
!        ! lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
!        ! mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec))
!        ! lambdaplus2mu_unrelaxed_elastic = poroelastcoef(3,1,kmato(ispec))
!        ! lambdalplusmul_unrelaxed_elastic = lambdal_unrelaxed_elastic + mul_unrelaxed_elastic
!
!        ! if( assign_external_model ) then
!        !   cpl = vpext(i_pnt,j_pnt,ispec)
!        !   csl = vsext(i_pnt,j_pnt,ispec)
!        !   rhol = rhoext(i_pnt,j_pnt,ispec)
!        !   mul_unrelaxed_elastic = rhol*csl*csl
!        !   lambdal_unrelaxed_elastic = rhol*cpl*cpl - TWO*mul_unrelaxed_elastic
!        !   lambdaplus2mu_unrelaxed_elastic = lambdal_unrelaxed_elastic + TWO*mul_unrelaxed_elastic
!        !   lambdalplusmul_unrelaxed_elastic = lambdal_unrelaxed_elastic + mul_unrelaxed_elastic
!        ! endif
!
!        ! do k = 1,NGLJ
!        !    dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j_pnt,ispec))*hprime_xx(i_pnt,k)
!        !    duy_dxi = duy_dxi + displ_elastic(2,ibool(k,j_pnt,ispec))*hprime_xx(i_pnt,k)
!        !    duz_dxi = duz_dxi + displ_elastic(3,ibool(k,j_pnt,ispec))*hprime_xx(i_pnt,k)
!        !    dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i_pnt,k,ispec))*hprime_zz(j_pnt,k)
!        !    duy_dgamma = duy_dgamma + displ_elastic(2,ibool(i_pnt,k,ispec))*hprime_zz(j_pnt,k)
!        !    duz_dgamma = duz_dgamma + displ_elastic(3,ibool(i_pnt,k,ispec))*hprime_zz(j_pnt,k)
!        ! enddo
!
!        ! xixl = xix(i_pnt,j_pnt,ispec)
!        ! xizl = xiz(i_pnt,j_pnt,ispec)
!        ! gammaxl = gammax(i_pnt,j_pnt,ispec)
!        ! gammazl = gammaz(i_pnt,j_pnt,ispec)
! 
!        ! ! derivatives of displacement
!        ! dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
!        ! dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl
! 
!        ! duy_dxl = duy_dxi*xixl + duy_dgamma*gammaxl
!        ! duy_dzl = duy_dxi*xizl + duy_dgamma*gammazl
! 
!        ! duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
!        ! duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl
!!!!     !mathematic for calculating moment tensor point source
!         !1.compute m_{kl} calculation
!         !!for isotropic material
!         !!m_xx = u_1*n_1*(lambda+2mu) + u_3*n_3*(lamba)
!         !!m_zz = u_3*n_3*(lambda+2mu) + u_1*n_1*(lamba)
!         !!m_xz = u_1*n_3*mu + u_3*n_1*mu
!         !!m_zx = m_xz
!         !!m_xy = u_2*n_1*mu !!should we consider m_{y?}
!         !!m_zy = u_2*n_3*mu
!!!!     tthat: f = - \nabla \cdot m. If we consider that
!!!!     )he material properties are continuous (or we don't inculde the ponts at the interface);
!!!!     )he geometry is smooth enough (i.e., curvature changes very slowly);
!!!!     e the partial derivative only effects on the displacement?
!!!!      n simplify the calculation (need to further concern)
!         !!f_1 = -{ \partial_1 m_11 + \partial_3 m_31}
!         !!    = -{(\partial_1 u_1)*n_1*(lambda+2mu) + (\partial_1 u_3)*n_3*(lamba)
!         !!        + (\partial_3 u_1)*n_3*mu + (\partial_3 u_3)*n_1*mu}
!
!         !!f_3 = -{ \partial_1 m_31 + \partial_3 m_33}
!         !!    = -{(\partial_1 u_1)*n_3*mu + (\partial_1 u_3)*n_1*mu
!         !!        + (\partial_3 u_3)*n_3*(lambda+2mu) + (\partial_3 u_1)*n_1*(lamba)}
!         !!f_2 = -{(\partial_1 u_2)*n_1*mu
!         !!        (\partial_3 u_2)*n_3*mu}
!
!  
!!!i     lentation, code
!         !m_xx(ipnt) = displ_elastic(1,iglob)*nx_bd_pnt_elastic(ipnt)*lambdaplus2mu_unrelaxed_elastic &
!         !+ displ_elastic(3,iglob)*nz_bd_pnt_elastic(ipnt)*lambdal_unrelaxed_elastic
!         !m_xz(ipnt) = displ_elastic(1,iglob)*nz_bd_pnt_elastic(ipnt)*mul_unrelaxed_elastic &
!         !+ displ_elastic(3,iglob)*nx_bd_pnt_elastic(ipnt)* mul_unrelaxed_elastic
!         !m_zz(ipnt) = displ_elastic(1,iglob)*nx_bd_pnt_elastic(ipnt)*lambdal_unrelaxed_elastic &
!         !+ displ_elastic(3,iglob)*nz_bd_pnt_elastic(ipnt)*lambdaplus2mu_unrelaxed_elastic
!         !!m_zx(ipnt) = m_xz(ipnt)
!
!         !m_f_bd_pnt_elastic(1,ispec_bd_pnt_elastic) = &
!         !   - (dux_dxl*nx_bd_pnt_elastic(ispec_bd_pnt_elastic)*lambdaplus2mu_unrelaxed_elastic &
!         !   + duz_dxl*nz_bd_pnt_elastic(ispec_bd_pnt_elastic)*lambdal_unrelaxed_elastic &
!         !   + dux_dzl*nz_bd_pnt_elastic(ispec_bd_pnt_elastic)*mul_unrelaxed_elastic &
!         !   + duz_dzl*nx_bd_pnt_elastic(ispec_bd_pnt_elastic)*mul_unrelaxed_elastic)
!         !
!         !m_f_bd_pnt_elastic(3,ispec_bd_pnt_elastic) = &
!         !   - (duz_dzl*nz_bd_pnt_elastic(ispec_bd_pnt_elastic)*lambdaplus2mu_unrelaxed_elastic & 
!         !   + dux_dzl*nx_bd_pnt_elastic(ispec_bd_pnt_elastic)*lambdal_unrelaxed_elastic &
!         !   + dux_dxl*nz_bd_pnt_elastic(ispec_bd_pnt_elastic)*mul_unrelaxed_elastic &
!         !   + duz_dxl*nx_bd_pnt_elastic(ispec_bd_pnt_elastic)*mul_unrelaxed_elastic)
!
!         !m_f_bd_pnt_elastic(2,ispec_bd_pnt_elastic) = &
!         !   - (duy_dxl*nx_bd_pnt_elastic(ispec_bd_pnt_elastic)*mul_unrelaxed_elastic &
!         !   + duy_dzl*nz_bd_pnt_elastic(ispec_bd_pnt_elastic)*mul_unrelaxed_elastic)
!
!        
!         exit loop1
!
!       endif !locate the corresponding recording point
!           
!   enddo loop1  
!     
!
! end subroutine record_bd_elmnt_elastic_reconst_m
!
! subroutine calculate_bd_elastic_reconst_m_f()
!
!   use specfem_par, only: it,xix,xiz,gammax,gammaz,jacobian,wxgll,wzgll,&
!                          ispec_bd_elmt_elastic,ispec_bd_elmt_elastic_i,ispec_bd_elmt_elastic_j,&
!                          nspec_bd_pnt_elastic,side_type_elastic,&
!                          ispec_bd_elmt_elastic_pure_edge,nspec_bd_elmt_elastic_pure_edge, & 
!                          ispec_bd_elmt_elastic_pure_side,&
!                          hprime_zz,hprime_xx,&
!                          m_xx,m_xz,m_zz,m_f
!   implicit none
!   include "constants.h"
!   
!   double precision, dimension(NGLLX) :: G11,G13,G31,G33
!   integer :: i,j,k,ispec_selected_m_f,ispec_bd_edge_elastic
!   double precision :: xixd,xizd,gammaxd,gammazd !for G_ik calculation
!   real(kind=CUSTOM_REAL) :: jacobian1D,weight,xxi,zxi,xgamma,zgamma !for integeral weight calculation
!
!   if( nspec_bd_pnt_elastic /= 0 ) then
!     do k = 1, nspec_bd_pnt_elastic
!
!   loop1: do ispec_bd_edge_elastic = 1, nspec_bd_elmt_elastic_pure_edge  
!        if ( ispec_bd_elmt_elastic(k) == ispec_bd_elmt_elastic_pure_edge(ispec_bd_edge_elastic) &
!           .and. side_type_elastic(k) == ispec_bd_elmt_elastic_pure_side(ispec_bd_edge_elastic) ) then
!
!              ispec_selected_m_f = ispec_bd_elmt_elastic(k)
!
!           if(it == 1) print *,'ispec = ', ispec_selected_m_f, '  side type: ', side_type_elastic(k)
!
!           if( side_type_elastic(k) == 'L') then
!              m_f(k,1) = 0.0
!              m_f(k,3) = 0.0
!              !!calculate the G_ik for each GLL point along the edge
!              do j = 1,NGLLZ
!                xixd    = xix(1,j,ispec_selected_m_f)
!                xizd    = xiz(1,j,ispec_selected_m_f)
!                gammaxd = gammax(1,j,ispec_selected_m_f)
!                gammazd = gammaz(1,j,ispec_selected_m_f)
!  
!                !they are G11(1,j) indeed
!                G11(j) = m_xx(ispec_bd_edge_elastic,j)*xixd+m_xz(ispec_bd_edge_elastic,j)*xizd 
!                G13(j) = m_xx(ispec_bd_edge_elastic,j)*gammaxd+m_xz(ispec_bd_edge_elastic,j)*gammazd
!                G31(j) = m_xz(ispec_bd_edge_elastic,j)*xixd+m_zz(ispec_bd_edge_elastic,j)*xizd
!                G33(j) = m_xz(ispec_bd_edge_elastic,j)*gammaxd+m_zz(ispec_bd_edge_elastic,j)*gammazd
!              enddo
!             i =1
!             do j = 1,NGLLZ
!                xgamma = - xiz(i,j,ispec_selected_m_f) * jacobian(i,j,ispec_selected_m_f)
!                zgamma = + xix(i,j,ispec_selected_m_f) * jacobian(i,j,ispec_selected_m_f)
!                jacobian1D = sqrt(xgamma**2 + zgamma**2)
!                weight = jacobian1D * wzgll(j)
!                !hprime_zz(j,i) is the derivative of ith Lagrange interpolant at jth GLL point
!                m_f(k,1) = m_f(k,1) + weight*G11(j)*hprime_zz(j,ispec_bd_elmt_elastic_j(k))
!                m_f(k,3) = m_f(k,3) + weight*G31(j)*hprime_zz(j,ispec_bd_elmt_elastic_j(k))
!                
!             enddo
!
!             j = ispec_bd_elmt_elastic_j(k)
!             xgamma = - xiz(i,j,ispec_selected_m_f) * jacobian(i,j,ispec_selected_m_f)
!             zgamma = + xix(i,j,ispec_selected_m_f) * jacobian(i,j,ispec_selected_m_f)
!             jacobian1D = sqrt(xgamma**2 + zgamma**2)
!             weight = jacobian1D * wzgll(j)
!             m_f(k,1) = m_f(k,1) + weight*G13(j)*hprime_xx(i,ispec_bd_elmt_elastic_i(k))
!             m_f(k,3) = m_f(k,3) + weight*G33(j)*hprime_xx(i,ispec_bd_elmt_elastic_i(k))
!             !m_f(k,1) =  
!             !m_f(k,3) = 
!             !m_f(k,2) = 
!           
!             exit loop1 
!           endif !left
!
!           if( side_type_elastic(k) == 'R' ) then
!              m_f(k,1) = 0.0
!              m_f(k,3) = 0.0
!              !!calculate the G_ik for each GLL point along the edge
!              i = NGLLX
!
!              do j = 1,NGLLZ
!                xixd    = xix(i,j,ispec_selected_m_f)
!                xizd    = xiz(i,j,ispec_selected_m_f)
!                gammaxd = gammax(i,j,ispec_selected_m_f)
!                gammazd = gammaz(i,j,ispec_selected_m_f)
!  
!                !they are G11(NGLLX,j) indeed
!                G11(j) = m_xx(ispec_bd_edge_elastic,j)*xixd+m_xz(ispec_bd_edge_elastic,j)*xizd 
!                G13(j) = m_xx(ispec_bd_edge_elastic,j)*gammaxd+m_xz(ispec_bd_edge_elastic,j)*gammazd
!                G31(j) = m_xz(ispec_bd_edge_elastic,j)*xixd+m_zz(ispec_bd_edge_elastic,j)*xizd
!                G33(j) = m_xz(ispec_bd_edge_elastic,j)*gammaxd+m_zz(ispec_bd_edge_elastic,j)*gammazd
!              enddo
!
!              do j = 1,NGLLZ
!                 xgamma = - xiz(i,j,ispec_selected_m_f) * jacobian(i,j,ispec_selected_m_f)
!                 zgamma = + xix(i,j,ispec_selected_m_f) * jacobian(i,j,ispec_selected_m_f)
!                 jacobian1D = sqrt(xgamma**2 + zgamma**2)
!                 weight = jacobian1D * wzgll(j)
!                 !hprime_zz(j,i) is the derivative of ith Lagrange interpolant at jth GLL point
!                 m_f(k,1) = m_f(k,1) + weight*G11(j)*hprime_zz(j,ispec_bd_elmt_elastic_j(k))
!                 m_f(k,3) = m_f(k,3) + weight*G31(j)*hprime_zz(j,ispec_bd_elmt_elastic_j(k))
!                 
!              enddo
!
!              j = ispec_bd_elmt_elastic_j(k)
!              xgamma = - xiz(i,j,ispec_selected_m_f) * jacobian(i,j,ispec_selected_m_f)
!              zgamma = + xix(i,j,ispec_selected_m_f) * jacobian(i,j,ispec_selected_m_f)
!              jacobian1D = sqrt(xgamma**2 + zgamma**2)
!              weight = jacobian1D * wzgll(j)
!              m_f(k,1) = m_f(k,1) + weight*G13(j)*hprime_xx(i,ispec_bd_elmt_elastic_i(k))
!              m_f(k,3) = m_f(k,3) + weight*G33(j)*hprime_xx(i,ispec_bd_elmt_elastic_i(k))
!           
!              exit loop1 
!            endif !'right'
!          
!            if( side_type_elastic(k) == 'B' ) then
!              m_f(k,1) = 0.0
!              m_f(k,3) = 0.0
!              j = 1
!              do i = 1,NGLLX
!                 xixd    = xix(i,j,ispec_selected_m_f)
!                 xizd    = xiz(i,j,ispec_selected_m_f)
!                 gammaxd = gammax(i,j,ispec_selected_m_f)
!                 gammazd = gammaz(i,j,ispec_selected_m_f)
!  
!                 !they are G11(i,1) indeed
!                 G11(i) = m_xx(ispec_bd_edge_elastic,i)*xixd+m_xz(ispec_bd_edge_elastic,i)*xizd 
!                 G13(i) = m_xx(ispec_bd_edge_elastic,i)*gammaxd+m_xz(ispec_bd_edge_elastic,i)*gammazd
!                 G31(i) = m_xz(ispec_bd_edge_elastic,i)*xixd+m_zz(ispec_bd_edge_elastic,i)*xizd
!                 G33(i) = m_xz(ispec_bd_edge_elastic,i)*gammaxd+m_zz(ispec_bd_edge_elastic,i)*gammazd
!              enddo
!
!              do i = 1,NGLLX
!                 xxi = + gammaz(i,j,ispec_selected_m_f) * jacobian(i,j,ispec_selected_m_f)
!                 zxi = - gammax(i,j,ispec_selected_m_f) * jacobian(i,j,ispec_selected_m_f)
!                 jacobian1D = sqrt(xxi**2 + zxi**2)
!                 weight = jacobian1D * wxgll(i)
!                 
!                 m_f(k,1) = m_f(k,1) + weight*G11(i)*hprime_xx(i,ispec_bd_elmt_elastic_i(k))
!                 m_f(k,3) = m_f(k,3) + weight*G31(i)*hprime_xx(i,ispec_bd_elmt_elastic_i(k))
!              enddo
!              i = ispec_bd_elmt_elastic_i(k)
!              xxi = + gammaz(i,j,ispec_selected_m_f) * jacobian(i,j,ispec_selected_m_f)
!              zxi = - gammax(i,j,ispec_selected_m_f) * jacobian(i,j,ispec_selected_m_f)
!              jacobian1D = sqrt(xxi**2 + zxi**2)
!              weight = jacobian1D * wxgll(i)
!              
!              m_f(k,1) = m_f(k,1) + weight*G13(i)*hprime_zz(j,ispec_bd_elmt_elastic_j(k))
!              m_f(k,3) = m_f(k,3) + weight*G33(i)*hprime_zz(j,ispec_bd_elmt_elastic_j(k))
!           
!              exit loop1 
!            endif ! bottom
!
!           if( side_type_elastic(k) == 'T' ) then
!              m_f(k,1) = 0.0
!              m_f(k,3) = 0.0
!              j = NGLLZ
!              do i = 1,NGLLX
!                 xixd    = xix(i,j,ispec_selected_m_f)
!                 xizd    = xiz(i,j,ispec_selected_m_f)
!                 gammaxd = gammax(i,j,ispec_selected_m_f)
!                 gammazd = gammaz(i,j,ispec_selected_m_f)
!  
!                 !they are G11(i,NGLLZ) indeed
!                 G11(i) = m_xx(ispec_bd_edge_elastic,i)*xixd+m_xz(ispec_bd_edge_elastic,i)*xizd 
!                 G13(i) = m_xx(ispec_bd_edge_elastic,i)*gammaxd+m_xz(ispec_bd_edge_elastic,i)*gammazd
!                 G31(i) = m_xz(ispec_bd_edge_elastic,i)*xixd+m_zz(ispec_bd_edge_elastic,i)*xizd
!                 G33(i) = m_xz(ispec_bd_edge_elastic,i)*gammaxd+m_zz(ispec_bd_edge_elastic,i)*gammazd
!              enddo
!
!              do i = 1,NGLLX
!                 xxi = + gammaz(i,j,ispec_selected_m_f) * jacobian(i,j,ispec_selected_m_f)
!                 zxi = - gammax(i,j,ispec_selected_m_f) * jacobian(i,j,ispec_selected_m_f)
!                 jacobian1D = sqrt(xxi**2 + zxi**2)
!                 weight = jacobian1D * wxgll(i)
!                 
!                 m_f(k,1) = m_f(k,1) + weight*G11(i)*hprime_xx(i,ispec_bd_elmt_elastic_i(k))
!                 m_f(k,3) = m_f(k,3) + weight*G31(i)*hprime_xx(i,ispec_bd_elmt_elastic_i(k))
!              enddo
!
!              i = ispec_bd_elmt_elastic_i(k)
!              xxi = + gammaz(i,j,ispec_selected_m_f) * jacobian(i,j,ispec_selected_m_f)
!              zxi = - gammax(i,j,ispec_selected_m_f) * jacobian(i,j,ispec_selected_m_f)
!              jacobian1D = sqrt(xxi**2 + zxi**2)
!              weight = jacobian1D * wxgll(i)
!              
!              m_f(k,1) = m_f(k,1) + weight*G13(i)*hprime_zz(j,ispec_bd_elmt_elastic_j(k))
!              m_f(k,3) = m_f(k,3) + weight*G33(i)*hprime_zz(j,ispec_bd_elmt_elastic_j(k))
!           
!              exit loop1 
!           endif !top
!
!        endif !the ispec_bd_edge_elastic is being located
!
!        enddo loop1 !!!locate ispec_bd_edge_elastic to provide the corresponding moment density tensor
!
!     enddo  !!!end the loop for all recording points
!   endif
! end subroutine calculate_bd_elastic_reconst_m_f
!-------------------------------------------------------------------------
 subroutine write_bd_pnts_reconst()

  use specfem_par, only: it,& !original para
                         fname,f_num,&
                         nspec_bd_pnt_elastic,&
                         trac_f,m_xx,m_xz,m_zz, &
                         record_nt1_reconst,record_nt2_reconst !control time step for recording
 
 
  implicit none
  include "constants.h"

  integer :: k
  integer :: ios
  integer :: length_unf_1
  !integer :: length_unf_2

  if (it < record_nt1_reconst .or. it > record_nt2_reconst ) return

  !inquire (iolength = length_unf_2) grad_pot_bd_pnt_acoustic(:,1),pot_dot_bd_pnt_acoustic(1)

  !for elastic 
  if( nspec_bd_pnt_elastic /= 0 )then
    !calculate moment density tensor point sources
    !call calculate_bd_elastic_reconst_m_f()
    !!!this is the recording length for unformatted recording
    inquire (iolength = length_unf_1) trac_f(:,1),m_xx(1),m_xx(1),m_xx(1)
    f_num=113
    write(fname,"('./OUTPUT_FILES/reconst_record/&
          &elastic_pnts/nt_',i6.6)")it

    !unformatted recording
    open(unit=f_num,file=trim(fname),access='direct',status='new',&
         action='write',iostat=ios,recl=length_unf_1) 
    if( ios /= 0 ) stop 'error saving values at recording points'

    do k = 1, nspec_bd_pnt_elastic
       write(f_num,rec=k) trac_f(:,k),m_xx(k),m_xz(k),m_zz(k)
    enddo

    close(f_num)
  
  endif 

  !for acoustic

 end subroutine write_bd_pnts_reconst
