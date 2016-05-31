!this subroutine will be called to record the boundary (local\global) info

 subroutine record_bd_elemnt_prediction()
 !this subroutine is used to store prediction value of velocity and pot_dot for the
 !selected bd recording element, since the values involved in the acceleration calculation
 !is prediction value but not correction value
    use specfem_par, only: ibool,veloc_elastic,potential_dot_acoustic,& !original para
                           nspec_bd_elmt_elastic_pure,nspec_bd_elmt_acoustic_pure,&
                           vel_bd_elastic, pot_dot_bd_acoustic,&
                           ispec_bd_elmt_elastic_pure,ispec_bd_elmt_acoustic_pure
   
    implicit none
    include "constants.h"

    integer :: i,j,k,iglob,ispec
     
    loop1:do k = 1,nspec_bd_elmt_elastic_pure
             ispec = ispec_bd_elmt_elastic_pure(k)
             do i = 1, NGLLX
                do j = 1, NGLLZ
                iglob = ibool(i,j,ispec)
                vel_bd_elastic(1,i,j,k) = veloc_elastic(1,iglob)
                vel_bd_elastic(2,i,j,k) = veloc_elastic(2,iglob)
                vel_bd_elastic(3,i,j,k) = veloc_elastic(3,iglob)
                enddo
             enddo
     enddo loop1
  
     loop2:do k =1,nspec_bd_elmt_acoustic_pure
              ispec = ispec_bd_elmt_acoustic_pure(k)
              do i = 1, NGLLX 
                 do j = 1, NGLLZ
                 iglob = ibool(i,j,ispec)
                 pot_dot_bd_acoustic(i,j,k) = potential_dot_acoustic(iglob)
                 enddo
              enddo
     enddo loop2
 end subroutine record_bd_elemnt_prediction

 subroutine record_bd_elmnt_elastic(ispec,i,j,&
            sigma_xx,sigma_xy,sigma_xz,sigma_zz,sigma_zy)

   use specfem_par, only: stress_bd_elastic,nspec_bd_elmt_elastic_pure,&
                          ispec_bd_elmt_elastic_pure

   implicit none
   include "constants.h"
   real(kind=CUSTOM_REAL) :: sigma_xx,sigma_xy,sigma_xz,sigma_zz,sigma_zy
   integer :: ispec,i,j,k   
 
   loop1:do k = 1,nspec_bd_elmt_elastic_pure
         if ( ispec_bd_elmt_elastic_pure(k) == ispec ) then
         !record the information for this element
         !  trac_bd_pnt_elastic(1,i,j,k) = tx_store  
         !  trac_bd_pnt_elastic(2,i,j,k) = ty_store
         !  trac_bd_pnt_elastic(3,i,j,k) = tz_store
       
         !is it a better way to store stress tensor, so that we can get rid of the 
         !possible confusion of normal vector(nx_pnt, nz_pnt)?
           stress_bd_elastic(1,i,j,k) = sigma_xx 
           stress_bd_elastic(2,i,j,k) = sigma_xy
           stress_bd_elastic(3,i,j,k) = sigma_xz
           stress_bd_elastic(4,i,j,k) = sigma_zz
           stress_bd_elastic(5,i,j,k) = sigma_zy

           exit loop1
         endif
   enddo loop1 
 end subroutine record_bd_elmnt_elastic


 subroutine record_bd_elmnt_acoustic(ispec,i,j,dux_dxl,dux_dzl)
  
   use specfem_par, only: grad_pot_bd_acoustic,nspec_bd_elmt_acoustic_pure,&
                          ispec_bd_elmt_acoustic_pure
   implicit none
   include "constants.h"
   real(kind=CUSTOM_REAL) :: dux_dxl,dux_dzl
   integer :: ispec,i,j,k

   loop2:do k = 1,nspec_bd_elmt_acoustic_pure
         if ( ispec_bd_elmt_acoustic_pure(k) == ispec ) then

            grad_pot_bd_acoustic(1,i,j,k) = dux_dxl
            grad_pot_bd_acoustic(2,i,j,k) = dux_dzl

            exit loop2
         endif
   enddo loop2
 end subroutine record_bd_elmnt_acoustic


 subroutine write_bd_pnts()

  use specfem_par, only: elastic,acoustic,it,& !original para
                         npnt,ispec_selected_bd_pnt,fname,f_num,&
                         ispec_bd_elmt_elastic_pure, ispec_bd_elmt_acoustic_pure,&
                         hxi_bd_store, hgammar_bd_store,&
                         stress_bd_elastic,vel_bd_elastic,grad_pot_bd_acoustic,pot_dot_bd_acoustic,&
                         stress_bd_pnt_elastic,vel_bd_pnt_elastic,trac_bd_pnt_elastic,&
                         nspec_bd_elmt_elastic_pure,nspec_bd_elmt_acoustic_pure,&
                         nspec_bd_pnt_elastic,nspec_bd_pnt_acoustic,&
                         grad_pot_bd_pnt_acoustic,pot_dot_bd_pnt_acoustic,&
                         nx_pnt,nz_pnt

  implicit none
  include "constants.h"

  integer :: ispec_bd_pnt_elastic = 1, ispec_bd_pnt_acoustic = 1
  integer :: ipnt,ispec,k,kk,i,j
  integer :: ios
  double precision :: hlagrange

  do ipnt = 1, npnt

     ispec =  ispec_selected_bd_pnt(ipnt)

     if ( elastic(ispec) ) then


       loop1: do k = 1,nspec_bd_elmt_elastic_pure
           !here we locate in which element the recording point locates
           if ( ispec_bd_elmt_elastic_pure(k) == ispec ) then

              stress_bd_pnt_elastic(:,ispec_bd_pnt_elastic) = 0.0
              vel_bd_pnt_elastic(:,ispec_bd_pnt_elastic) = 0.0
              !we use the lagrange interpolation to calculate
              !the value at the point based on all the GLL points
              !in that element
              do i = 1,NGLLX
                 do j = 1,NGLLZ
                    hlagrange = hxi_bd_store(ipnt,i)*hgammar_bd_store(ipnt,j)

                    stress_bd_pnt_elastic(:,ispec_bd_pnt_elastic) = &
                    stress_bd_pnt_elastic(:,ispec_bd_pnt_elastic) + &
                    stress_bd_elastic(:,i,j,k)*hlagrange 

                    vel_bd_pnt_elastic(:,ispec_bd_pnt_elastic) = &
                    vel_bd_pnt_elastic(:,ispec_bd_pnt_elastic) + &
                    vel_bd_elastic(:,i,j,k)*hlagrange
                 enddo
              enddo

              exit loop1
           endif
        enddo loop1

        ispec_bd_pnt_elastic = ispec_bd_pnt_elastic + 1

     endif!elastic pnt

     if ( acoustic(ispec) ) then

        grad_pot_bd_pnt_acoustic(:,ispec_bd_pnt_acoustic) = 0.0
        pot_dot_bd_pnt_acoustic(ispec_bd_pnt_acoustic) = 0.0
 
       loop2: do kk = 1,nspec_bd_elmt_acoustic_pure
            if ( ispec_bd_elmt_acoustic_pure(kk) == ispec ) then
               do i = 1,NGLLX
                  do j = 1,NGLLZ
                     hlagrange = hxi_bd_store(ipnt,i)*hgammar_bd_store(ipnt,j)

                     grad_pot_bd_pnt_acoustic(:,ispec_bd_pnt_acoustic) = &
                     grad_pot_bd_pnt_acoustic(:,ispec_bd_pnt_acoustic) + &
                     grad_pot_bd_acoustic(:,i,j,kk)*hlagrange

                     pot_dot_bd_pnt_acoustic(ispec_bd_pnt_acoustic) = &
                     pot_dot_bd_pnt_acoustic(ispec_bd_pnt_acoustic) + &
                     pot_dot_bd_acoustic(i,j,kk)*hlagrange
                  enddo
               enddo

               exit loop2
            endif
        enddo loop2

        ispec_bd_pnt_acoustic = ispec_bd_pnt_acoustic + 1

     endif!acoustic pnt

  enddo

  !store valuse at recording points at every time step
  !we firstly compute the traction for those recording points
  !trac_x
  trac_bd_pnt_elastic(1,:) = nx_pnt(:)*stress_bd_pnt_elastic(:,1) + &
                             nz_pnt(:)*stress_bd_pnt_elastic(:,3)
  !trac_z
  trac_bd_pnt_elastic(2,:) = nx_pnt(:)*stress_bd_pnt_elastic(:,3) + &
                             nz_pnt(:)*stress_bd_pnt_elastic(:,4)
  !trac_y
  trac_bd_pnt_elastic(3,:) = nx_pnt(:)*stress_bd_pnt_elastic(:,2) + &
                             nz_pnt(:)*stress_bd_pnt_elastic(:,5)

  f_num=113
  !for elastic 
  write(fname,"('./OUTPUT_FILES/bg_record/&
        &elastic_pnts/nt_',i6.6)")it
  open(unit=f_num,file=trim(fname),status='new',&
       action='write',iostat=ios) 
  if( ios /= 0 ) stop 'error saving values at recording points'

  do k = 1, nspec_bd_pnt_elastic
     write(f_num,111) trac_bd_pnt_elastic(:,k),vel_bd_pnt_elastic(:,k)
  enddo

  close(f_num)
  !for acoustic
  write(fname,"('./OUTPUT_FILES/bg_record/&
        &acoustic_pnts/nt_',i6.6)")it
  open(unit=f_num,file=trim(fname),status='new',&
       action='write',iostat=ios) 
  if( ios /= 0 ) stop 'error saving values at recording points'
  
  do kk = 1, nspec_bd_pnt_acoustic
     write(f_num,112) grad_pot_bd_pnt_acoustic(:,kk),pot_dot_bd_pnt_acoustic(kk)
  enddo

  close(f_num)

  111 format(6(es12.4,2x)) !112 column
  112 format(3(es12.4,2x)) !36 column
 end subroutine write_bd_pnts
