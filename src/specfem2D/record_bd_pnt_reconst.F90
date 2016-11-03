!!!lcx:this file is following 'record_bd_pnt.F90', recording the information along the local model boundary
!!(i.e., traction T_i and f_i = \partial_k m_{ki}, here m_{kl}=u_i n_j C_{ijkl} is the moment density tensor)


!!this subroutine is to record the stress tensor along the local model boundary. 
!!And It seems we can just keep 'record_bd_elmnt_elastic', but not need a new one.
 subroutine record_bd_elmnt_elastic_reconst(ispec,i,j,&
            sigma_xx,sigma_xy,sigma_xz,sigma_zz,sigma_zy)

   use specfem_par, only: elastic,npnt,ispec_selected_bd_pnt,ispec_selected_bd_pnt_i,ispec_selected_bd_pnt_j,&
                          trac_bd_pnt_elastic_reconst,nx_bd_pnt_elastic,nz_bd_pnt_elastic,&
                          it,record_nt1,record_nt2 !control time step for recording

   implicit none
   include "constants.h"
   real(kind=CUSTOM_REAL), intent(in) :: sigma_xx,sigma_xy,sigma_xz,sigma_zz,sigma_zy
   integer, intent(in) :: ispec,i,j
   integer :: ipnt,ispec_bd_pnt_elastic   

   if (it < record_nt1 .or. it > record_nt2 ) return

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

           else if( side_type_elastic(ispec_bd_pnt_elastic) == 'R' ) !Rigth

             xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
             zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
             jacobian1D = sqrt(xgamma**2 + zgamma**2)
             weight = jacobian1D * wzgll(j)

             trac_f(:,ispec_bd_pnt_elastic) = trac_bd_pnt_elastic_reconst(:,ispec_bd_pnt_elastic)*weight
             
           else if( side_type_elastic(ispec_bd_pnt_elastic) == 'B' ) !Bottom

             xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
             zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
             jacobian1D = sqrt(xxi**2 + zxi**2)
             weight = jacobian1D * wxgll(i)

             trac_f(:,ispec_bd_pnt_elastic) = trac_bd_pnt_elastic_reconst(:,ispec_bd_pnt_elastic)*weight

           else if( side_type_elastic(ispec_bd_pnt_elastic) == 'T' ) !Top

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

 end subroutine record_bd_elmnt_elastic_reconst


 subroutine record_bd_elmnt_elastic_reconst_m(ispec,i,j,&
                                              dux_dxl,dux_dzl,duz_dxl,duz_dzl,duy_dxl,duy_dzl)

   use specfem_par, only: it,& !original para
                          ispec_selected_bd_pnt, ispec_selected_bd_pnt_i, ispec_selected_bd_pnt_j,&
                          record_nt1,record_nt2 !control time step for recording

   implicit none
   include "constants.h"

   integer, intent(in) :: ispec,i,j
   real(kind=CUSTOM_REAL), intent(in) :: dux_dxl,duy_dxl,duz_dxl,dux_dzl,duy_dzl,duz_dzl
   !integer :: k   
   integer :: ispec_bd_pnt_elastic

   if (it < record_nt1 .or. it > record_nt2 ) return
   
   ispec_bd_pnt_elastic = 0

   loop1:do ispec_bd_pnt_elastic = 1, nspec_bd_pnt_elastic
      
      !locate the corresponding recording point
         ispec_bd_pnt_elastic = ispec_bd_pnt_elastic + 1
         !do k = 1,nspec_bd_elmt_elastic_pure
         ! I think it will be an economic way to read the i,j of the recording pointfrom file 'boundary_points'
         !ispec = ispec_selected_bd_pnt(ipnt)
         !i_pnt = ispec_selected_bd_pnt_i(ipnt)
         !j_pnt = ispec_selected_bd_pnt_j(ipnt)
         !iglob = ibool(i_pnt,j_pnt,ispec)
         
      !point is not the corresponding recording point. search the next one 

         if ( ispec_bd_elmt_elastic(ispec_bd_pnt_elastic) == ispec .and. ispec_bd_elmt_elastic_i(ispec_bd_pnt_elastic) == i &
            .and. ispec_bd_elmt_elastic_j(ispec_bd_pnt_elastic) == j ) then

        ! ! get unrelaxed elastic parameters of current spectral element
        ! lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
        ! mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec))
        ! lambdaplus2mu_unrelaxed_elastic = poroelastcoef(3,1,kmato(ispec))
        ! lambdalplusmul_unrelaxed_elastic = lambdal_unrelaxed_elastic + mul_unrelaxed_elastic

        ! if( assign_external_model ) then
        !   cpl = vpext(i_pnt,j_pnt,ispec)
        !   csl = vsext(i_pnt,j_pnt,ispec)
        !   rhol = rhoext(i_pnt,j_pnt,ispec)
        !   mul_unrelaxed_elastic = rhol*csl*csl
        !   lambdal_unrelaxed_elastic = rhol*cpl*cpl - TWO*mul_unrelaxed_elastic
        !   lambdaplus2mu_unrelaxed_elastic = lambdal_unrelaxed_elastic + TWO*mul_unrelaxed_elastic
        !   lambdalplusmul_unrelaxed_elastic = lambdal_unrelaxed_elastic + mul_unrelaxed_elastic
        ! endif

        ! do k = 1,NGLJ
        !    dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j_pnt,ispec))*hprime_xx(i_pnt,k)
        !    duy_dxi = duy_dxi + displ_elastic(2,ibool(k,j_pnt,ispec))*hprime_xx(i_pnt,k)
        !    duz_dxi = duz_dxi + displ_elastic(3,ibool(k,j_pnt,ispec))*hprime_xx(i_pnt,k)
        !    dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i_pnt,k,ispec))*hprime_zz(j_pnt,k)
        !    duy_dgamma = duy_dgamma + displ_elastic(2,ibool(i_pnt,k,ispec))*hprime_zz(j_pnt,k)
        !    duz_dgamma = duz_dgamma + displ_elastic(3,ibool(i_pnt,k,ispec))*hprime_zz(j_pnt,k)
        ! enddo

        ! xixl = xix(i_pnt,j_pnt,ispec)
        ! xizl = xiz(i_pnt,j_pnt,ispec)
        ! gammaxl = gammax(i_pnt,j_pnt,ispec)
        ! gammazl = gammaz(i_pnt,j_pnt,ispec)
 
        ! ! derivatives of displacement
        ! dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
        ! dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl
 
        ! duy_dxl = duy_dxi*xixl + duy_dgamma*gammaxl
        ! duy_dzl = duy_dxi*xizl + duy_dgamma*gammazl
 
        ! duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
        ! duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl
!!!     !mathematic for calculating moment tensor point source
         !1.compute m_{kl} calculation
         !!for isotropic material
         !!m_xx = u_1*n_1*(lambda+2mu) + u_3*n_3*(lamba)
         !!m_zz = u_3*n_3*(lambda+2mu) + u_1*n_1*(lamba)
         !!m_xz = u_1*n_3*mu + u_3*n_1*mu
         !!m_zx = m_xz
         !!m_xy = u_2*n_1*mu !!should we consider m_{y?}
         !!m_zy = u_2*n_3*mu
!!!     tthat: f = - \nabla \cdot m. If we consider that
!!!     )he material properties are continuous (or we don't inculde the ponts at the interface);
!!!     )he geometry is smooth enough (i.e., curvature changes very slowly);
!!!     e the partial derivative only effects on the displacement?
!!!      n simplify the calculation (need to further concern)
         !!f_1 = -{ \partial_1 m_11 + \partial_3 m_31}
         !!    = -{(\partial_1 u_1)*n_1*(lambda+2mu) + (\partial_1 u_3)*n_3*(lamba)
         !!        + (\partial_3 u_1)*n_3*mu + (\partial_3 u_3)*n_1*mu}

         !!f_3 = -{ \partial_1 m_31 + \partial_3 m_33}
         !!    = -{(\partial_1 u_1)*n_3*mu + (\partial_1 u_3)*n_1*mu
         !!        + (\partial_3 u_3)*n_3*(lambda+2mu) + (\partial_3 u_1)*n_1*(lamba)}
         !!f_2 = -{(\partial_1 u_2)*n_1*mu
         !!        (\partial_3 u_2)*n_3*mu}

  
!!i     lentation, code
         !m_xx(ipnt) = displ_elastic(1,iglob)*nx_bd_pnt_elastic(ipnt)*lambdaplus2mu_unrelaxed_elastic &
         !+ displ_elastic(3,iglob)*nz_bd_pnt_elastic(ipnt)*lambdal_unrelaxed_elastic
         !m_xz(ipnt) = displ_elastic(1,iglob)*nz_bd_pnt_elastic(ipnt)*mul_unrelaxed_elastic &
         !+ displ_elastic(3,iglob)*nx_bd_pnt_elastic(ipnt)* mul_unrelaxed_elastic
         !m_zz(ipnt) = displ_elastic(1,iglob)*nx_bd_pnt_elastic(ipnt)*lambdal_unrelaxed_elastic &
         !+ displ_elastic(3,iglob)*nz_bd_pnt_elastic(ipnt)*lambdaplus2mu_unrelaxed_elastic
         !!m_zx(ipnt) = m_xz(ipnt)

         m_f_bd_pnt_elastic(1,ispec_bd_pnt_elastic) = - (dux_dxl*nx_bd_pnt_elastic(ispec_bd_pnt_elastic)*lambdaplus2mu_unrelaxed_elastic &
            + duz_dxl*nz_bd_pnt_elastic(ispec_bd_pnt_elastic)*lambdal_unrelaxed_elastic &
            + dux_dzl*nz_bd_pnt_elastic(ispec_bd_pnt_elastic)*mul_unrelaxed_elastic &
            + duz_dzl*nx_bd_pnt_elastic(ispec_bd_pnt_elastic)*mul_unrelaxed_elastic)
         
         m_f_bd_pnt_elastic(3,ispec_bd_pnt_elastic) = - (duz_dzl*nz_bd_pnt_elastic(ispec_bd_pnt_elastic)*lambdaplus2mu_unrelaxed_elastic & 
            + dux_dzl*nx_bd_pnt_elastic(ispec_bd_pnt_elastic)*lambdal_unrelaxed_elastic &
            + dux_dxl*nz_bd_pnt_elastic(ispec_bd_pnt_elastic)*mul_unrelaxed_elastic &
            + duz_dxl*nx_bd_pnt_elastic(ispec_bd_pnt_elastic)*mul_unrelaxed_elastic)

         m_f_bd_pnt_elastic(2,ispec_bd_pnt_elastic) = - (duy_dxl*nx_bd_pnt_elastic(ispec_bd_pnt_elastic)*mul_unrelaxed_elastic &
            + duy_dzl*nz_bd_pnt_elastic(ispec_bd_pnt_elastic)*mul_unrelaxed_elastic)

        
         exit loop1

       endif !locate the corresponding recording point
           
   enddo loop1  
     

 end subroutine record_bd_elmnt_elastic_reconst_m


 subroutine write_bd_pnts_reconst()

  use specfem_par, only: elastic,it,& !original para
                        npnt,ispec_selected_bd_pnt,fname,f_num,&
                        trac_f,&!m_f_bd_pnt_elastic,&
                        record_nt1,record_nt2 !control time step for recording
 
 
  implicit none
  include "constants.h"

  integer :: k
  integer :: ios
  integer :: length_unf_1
  !integer :: length_unf_2

  if (it < record_nt1 .or. it > record_nt2 ) return

  !inquire (iolength = length_unf_2) grad_pot_bd_pnt_acoustic(:,1),pot_dot_bd_pnt_acoustic(1)

  !for elastic 
  if( nspec_bd_pnt_elastic /= 0 )then
    !!!this is the recording length for unformatted recording
    inquire (iolength = length_unf_1) trac_f(:,1)!,m_f_bd_pnt_elastic(:,1)
    f_num=113
    write(fname,"('./OUTPUT_FILES/reconst_record/&
          &elastic_pnts/nt_',i6.6)")it

    !unformatted recording
    open(unit=f_num,file=trim(fname),access='direct',status='new',&
         action='write',iostat=ios,recl=length_unf_1) 
    if( ios /= 0 ) stop 'error saving values at recording points'

    do k = 1, nspec_bd_pnt_elastic
       write(f_num,rec=k) trac_f(:,k),m_f_bd_pnt_elastic(:,k)
    enddo

    close(f_num)
  
  endif 


 end subroutine write_bd_pnts_reconst
