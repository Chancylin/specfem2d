!!!by lcx: this subroutine will be called by setup_recording_bd() in prepare_timerun_body.F90
!!!this subroutine will determine which elements the recording points
!belong to. And whether it is in fluid or solid region.
      subroutine locate_recording_point(ibool,coord,nspec,nglob,xigll,zigll, &
                          npnt, bd_pnt_xval,bd_pnt_zval,ispec_selected_bd_pnt, &
                          xi_bd_pnt,gamma_bd_pnt,&
                          coorg,knods,ngnod,npgeo, &
                          x_final_bd_pnt,z_final_bd_pnt)

      implicit none

      include "constants.h"

      integer npnt,nspec,nglob,ngnod,npgeo
    
      integer knods(ngnod,nspec)
      double precision coorg(NDIM,npgeo)
    
      integer, dimension(NGLLX,NGLLZ,nspec) :: ibool
    
    ! array containing coordinates of the points
      double precision coord(NDIM,nglob)
    
      integer ipnt,i,j,ispec,iglob,iter_loop,ix_initial_guess,iz_initial_guess
      double precision xi,gamma,dx,dz,dxi,dgamma
    
    ! Gauss-Lobatto-Legendre points of integration
      double precision xigll(NGLLX)
      double precision zigll(NGLLZ)
    
      double precision x,z,xix,xiz,gammax,gammaz,jacobian
    
    ! use dynamic allocation
      double precision distmin_squared
    
    ! recording_bd point information
      integer, dimension(npnt) :: ispec_selected_bd_pnt 
      double precision, dimension(npnt) :: xi_bd_pnt,gamma_bd_pnt
      double precision, dimension(npnt) :: bd_pnt_xval,bd_pnt_zval 
    
    ! tangential detection
      double precision, dimension(npnt)  :: x_final_bd_pnt, z_final_bd_pnt
    
      integer  :: ierror

      open(unit=1,file='DATA/boundary_points',status='old',action='read')
      !loop over all points
      ! loop only on points inside the element
      ! exclude edges to ensure this point is not shared with other elements
      do ipnt=1,npnt

        distmin_squared = HUGEVAL

        do ispec=1,nspec
          do j=2,NGLLZ-1
            do i=2,NGLLX-1
              iglob = ibool(i,j,ispec)
              dist_squared = (bd_pnt_xval(ipnt)-dble(coord(1,iglob)))**2 + &
              (bd_pnt_zval(ipnt)-dble(coord(2,iglob)))**2

              if(dist_squared < distmin_squared) then
                distmin_squared = dist_squared
                ispec_selected_bd_pnt(ipnt) = ispec
                ix_initial_guess = i
                iz_initial_guess = j
              endif
            enddo
          enddo

     
        enddo !end the loop for all elements

        !find the best (xi, gamma) for each recording point
        !start using intial guess in xi and gamma
          xi = xigll(ix_initial_guess)
          gamma = zigll(iz_initial_guess)

         do iter_loop = 1,NUM_ITER
           ! compute coordinates of the new point and derivatives dxi/dx, dxi/dz
           call recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian, &
                       coorg,knods,ispec_selected_bd_pnt(ipnt),ngnod,nspec,npgeo, &
                       .true.)
     
           ! compute distance to target location
           dx = - (x - bd_pnt_xval(ipnt))
           dz = - (z - bd_pnt_zval(ipnt))
     
           ! compute increments
           dxi  = xix*dx + xiz*dz
           dgamma = gammax*dx + gammaz*dz
     
           ! update values
           xi = xi + dxi
           gamma = gamma + dgamma
     
           ! impose that we stay in that element
           ! (useful if user gives a receiver outside the mesh for instance)
           ! we can go slightly outside the [1,1] segment since with finite elements
           ! the polynomial solution is defined everywhere
           ! this can be useful for convergence of itertive scheme with distorted elements
           if (xi > 1.10d0) xi = 1.10d0
           if (xi < -1.10d0) xi = -1.10d0
           if (gamma > 1.10d0) gamma = 1.10d0
           if (gamma < -1.10d0) gamma = -1.10d0
     
         ! end of non linear iterations
         enddo
          ! compute final coordinates of point found

         !by lcx: I think the basic idea here is we want to find the best 
         !(xi,gammar) for the point (x_bd_pnt,z_bd_pnt). Hence we need to 
         ! recompute the jacobian for several iterations. Once a point
         ! which is closed enough is found in (xi,gammar) (nature) coordinate
         ! we will reproject it back to the real coordinate

         call recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian, &
                     coorg,knods,ispec_selected_bd_pnt(ipnt),ngnod,nspec,npgeo, &
                     .true.)
     
         ! store xi,gamma of point found
         xi_bd_pnt(ipnt) = xi
         gamma_bd_pnt(ipnt) = gamma
     
     
         x_final_bd_pnt(ipnt) = x
         z_final_bd_pnt(ipnt) = z
     
      enddo 
      
      close(1)

      end subroutine locate_recording_point
