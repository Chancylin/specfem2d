!this subroutine is used to supply the store traction at boundary and
!apply the absorbing boudnary condition to the scatter waivefield
! (i.e., scatter = totoal - background)
subroutine absorb_scatter_field_solid(ispec,ispecabs,cpl,csl,rhol,accel_elastic,veloc_elastic)
  use specfem_par, only: nglob,ibool,wzgll,wxgll,&
                         xiz,xix,gammaz,gammax,jacobian,&
                         codeabs,codeabs_corner,assign_external_model,&
                         rhoext,vpext,vsext 
  implicit none
  include "constants.h"

  integer :: i,j,iglob  
  integer,intent(in) :: ispec,ispecabs
  real(kind=CUSTOM_REAL), dimension(3,nglob) :: accel_elastic,veloc_elastic
  real(kind=CUSTOM_REAL) :: nx,nz,vx,vy,vz,vn,rho_vp,rho_vs,tx,ty,tz,weight,xxi,zxi,xgamma,zgamma,jacobian1D
  real(kind=CUSTOM_REAL) :: cpl,csl,rhol
  real(kind=CUSTOM_REAL) :: veloc_x_store,veloc_y_store,veloc_z_store,tx_store,ty_store,tz_store
  
  !---left absorbing boundary
  if( codeabs(IEDGE4,ispecabs) ) then
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

       accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tx + tx_store)*weight 
       accel_elastic(2,iglob) = accel_elastic(2,iglob) - (ty + ty_store)*weight
       accel_elastic(3,iglob) = accel_elastic(3,iglob) - (tz + tz_store)*weight
    enddo
  endif
  !---right absorbing boundary
  if( codeabs(IEDGE2,ispecabs) ) then 
    i = NGLLX
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
       nx = + zgamma / jacobian1D
       nz = - xgamma / jacobian1D

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

       accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tx + tx_store)*weight
       accel_elastic(2,iglob) = accel_elastic(2,iglob) - (ty + ty_store)*weight
       accel_elastic(3,iglob) = accel_elastic(3,iglob) - (tz + tz_store)*weight
       
    enddo
  endif
  !---bottom absorbing boundary
  if( codeabs(IEDGE1,ispecabs) ) then
    j = 1
    do i = 1,NGLLX
       iglob = ibool(i,j,ispec)
       ! external velocity model
       if( assign_external_model ) then
         cpl = vpext(i,j,ispec)
         csl = vsext(i,j,ispec)
         rhol = rhoext(i,j,ispec)
       endif
       
       rho_vp = rhol*cpl
       rho_vs = rhol*csl

       xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
       zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
       jacobian1D = sqrt(xxi**2 + zxi**2)
       nx = + zxi / jacobian1D
       nz = - xxi / jacobian1D

       weight = jacobian1D * wxgll(i)

       call pnt_info_interpl_solid(iglob,veloc_x_store,veloc_y_store,veloc_z_store,&
                             tx_store,ty_store,tz_store)

       vx = veloc_elastic(1,iglob) - veloc_x_store
       vy = veloc_elastic(2,iglob) - veloc_y_store
       vz = veloc_elastic(3,iglob) - veloc_z_store

       vn = nx*vx+nz*vz

       tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
       ty = rho_vs*vy
       tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

       if( (codeabs_corner(1,ispecabs) .and. i == 1) .or. (codeabs_corner(2,ispecabs) .and. i == NGLLX) ) then
         tx = 0._CUSTOM_REAL; ty = 0._CUSTOM_REAL; tz = 0._CUSTOM_REAL
         tx_store = 0._CUSTOM_REAL; ty_store = 0._CUSTOM_REAL; tz_store = 0._CUSTOM_REAL
       endif
       
       accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tx + tx_store)*weight
       accel_elastic(2,iglob) = accel_elastic(2,iglob) - (ty + ty_store)*weight
       accel_elastic(3,iglob) = accel_elastic(3,iglob) - (tz + tz_store)*weight
       
    enddo
  endif
  !---top absorbing boundary
  if( codeabs(IEDGE3,ispecabs) ) then
    j = NGLLZ
    do i = 1,NGLLX
       iglob = ibool(i,j,ispec)
       ! external velocity model
       if( assign_external_model ) then
         cpl = vpext(i,j,ispec)
         csl = vsext(i,j,ispec)
         rhol = rhoext(i,j,ispec)
       endif

       
       rho_vp = rhol*cpl
       rho_vs = rhol*csl

       xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
       zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
       jacobian1D = sqrt(xxi**2 + zxi**2)
       nx = - zxi / jacobian1D
       nz = + xxi / jacobian1D

       weight = jacobian1D * wxgll(i)

       call pnt_info_interpl_solid(iglob,veloc_x_store,veloc_y_store,veloc_z_store,&
                             tx_store,ty_store,tz_store)

       vx = veloc_elastic(1,iglob) - veloc_x_store
       vy = veloc_elastic(2,iglob) - veloc_y_store
       vz = veloc_elastic(3,iglob) - veloc_z_store

       vn = nx*vx+nz*vz

       tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
       ty = rho_vs*vy
       tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

       
       if( (codeabs_corner(3,ispecabs) .and. i == 1) .or. (codeabs_corner(4,ispecabs) .and. i == NGLLX) ) then
         tx = 0._CUSTOM_REAL; ty = 0._CUSTOM_REAL; tz = 0._CUSTOM_REAL
         tx_store = 0._CUSTOM_REAL; ty_store = 0._CUSTOM_REAL; tz_store = 0._CUSTOM_REAL
       endif

       accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tx + tx_store)*weight
       accel_elastic(2,iglob) = accel_elastic(2,iglob) - (ty + ty_store)*weight
       accel_elastic(3,iglob) = accel_elastic(3,iglob) - (tz + tz_store)*weight

    enddo
  endif
  
end subroutine absorb_scatter_field_solid


subroutine absorb_scatter_field_fluid(ispec,ispecabs,cpl,rhol,potential_dot_dot_acoustic,potential_dot_acoustic) 
  use specfem_par, only: ibool,nglob,wzgll,wxgll,&
                         xiz,xix,gammaz,gammax,jacobian,&
                         codeabs,codeabs_corner,assign_external_model,&
                         rhoext,vpext,&
                         ibegin_edge1,iend_edge1,ibegin_edge3,iend_edge3, &
                         ibegin_edge4,iend_edge4,ibegin_edge2,iend_edge2
  implicit none
  include "constants.h"

  integer,intent(in) :: ispec,ispecabs
  real(kind=CUSTOM_REAL), dimension(nglob) :: potential_dot_dot_acoustic,potential_dot_acoustic
  real(kind=CUSTOM_REAL) :: nx,nz,weight,xxi,zxi,xgamma,zgamma,jacobian1D
  real(kind=CUSTOM_REAL) :: cpl,rhol
  real(kind=CUSTOM_REAL) :: grad_pot_x_store,grad_pot_z_store,potential_dot_acoustic_store
  integer :: i,j,ibegin,iend,jbegin,jend,iglob

  !test by lcx
  !print *,'test to see the element number: ', ispec
  !if(codeabs(IEDGE4,ispecabs)) then
  !         print *,'element number: ',ispec,' and element type is L'
  !  else if( codeabs(IEDGE2,ispecabs) )then
  !         print *,'element number: ',ispec,' and element type is R'
  !  else if( codeabs(IEDGE1,ispecabs) )then
  !         print *,'element number: ',ispec,' and element type is B'
  !  else if( codeabs(IEDGE3,ispecabs) )then
  !         print *,'element number: ',ispec,' and element type is T'
  !endif
  !--- left absorbing boundary
  if( codeabs(IEDGE4,ispecabs) ) then
    i = 1
    jbegin = ibegin_edge4(ispecabs)
    jend = iend_edge4(ispecabs)
    do j = jbegin,jend
      iglob = ibool(i,j,ispec)
      ! external velocity model
      if( assign_external_model ) then
        cpl = vpext(i,j,ispec)
        rhol = rhoext(i,j,ispec)
      endif
      xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
      zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
      jacobian1D = sqrt(xgamma ** 2 + zgamma ** 2)
      nx = - zgamma / jacobian1D
      nz = + xgamma / jacobian1D
      weight = jacobian1D * wzgll(j)

      call pnt_info_interpl_fluid(iglob,grad_pot_x_store,grad_pot_z_store,potential_dot_acoustic_store)

      !sampling as test
     ! if ( ispec == 2 )then
        ! print *,'j = ',j
        ! print *,'grad_pot_x, grad_pot_z,potential_dot',&
        !         grad_pot_x_store, grad_pot_z_store, potential_dot_acoustic_store
        ! print *,'nx,nz = ',nx,nz
        ! print *,'potential_dot_dot_acoustic(iglob) = ',potential_dot_dot_acoustic(iglob)
        ! print *,'potential_dot_acoustic(iglob) = ',potential_dot_acoustic(iglob)
        ! print *,'weight = ',weight
        ! print *,'rhol = ',rhol
        ! print *,'cpl = ',cpl
     ! endif
      ! adds absorbing boundary contribution
      potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - &
      ( nx*grad_pot_x_store + nz*grad_pot_z_store ) * weight/rhol - & !backgorund field
      ( potential_dot_acoustic(iglob) - potential_dot_acoustic_store ) * weight/cpl/rhol !scattered field

    enddo
  endif

  !--- right absorbing boundary
  if( codeabs(IEDGE2,ispecabs) ) then
    i = NGLLX
    jbegin = ibegin_edge2(ispecabs)
    jend = iend_edge2(ispecabs)
    do j = jbegin,jend
      iglob = ibool(i,j,ispec)
      ! external velocity model
      if( assign_external_model ) then
        cpl = vpext(i,j,ispec)
        rhol = rhoext(i,j,ispec)
      endif
      xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
      zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
      jacobian1D = sqrt(xgamma ** 2 + zgamma ** 2)
      nx = + zgamma / jacobian1D
      nz = - xgamma / jacobian1D
      weight = jacobian1D * wzgll(j)

      call pnt_info_interpl_fluid(iglob,grad_pot_x_store,grad_pot_z_store,potential_dot_acoustic_store)

      ! adds absorbing boundary contribution
      potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - &
      ( nx*grad_pot_x_store + nz*grad_pot_z_store ) * weight/rhol - & !backgorund field
      ( potential_dot_acoustic(iglob) - potential_dot_acoustic_store ) * weight/cpl/rhol !scattered field
    enddo
  endif

  !--- bottom absorbing boundary
  if( codeabs(IEDGE1,ispecabs) ) then
    j = 1
    ibegin = ibegin_edge1(ispecabs)
    iend = iend_edge1(ispecabs)
    ! exclude corners to make sure there is no contradiction on the normal
    if( codeabs_corner(1,ispecabs)) ibegin = 2
    if( codeabs_corner(2,ispecabs)) iend = NGLLX-1
    do i = ibegin,iend
      iglob = ibool(i,j,ispec)
      ! external velocity model
      if( assign_external_model ) then
        cpl = vpext(i,j,ispec)
        rhol = rhoext(i,j,ispec)
      endif
      xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
      zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
      jacobian1D = sqrt(xxi ** 2 + zxi ** 2)
      nx = + zxi / jacobian1D
      nz = - xxi / jacobian1D
      weight = jacobian1D * wxgll(i)

      call pnt_info_interpl_fluid(iglob,grad_pot_x_store,grad_pot_z_store,potential_dot_acoustic_store)

      ! adds absorbing boundary contribution
      potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - &
      ( nx*grad_pot_x_store + nz*grad_pot_z_store ) * weight/rhol - & !backgorund field
      ( potential_dot_acoustic(iglob) - potential_dot_acoustic_store ) * weight/cpl/rhol !scattered field
    
    enddo
  endif 

  !--- top absorbing boundary
  if( codeabs(IEDGE3,ispecabs) ) then
    j = NGLLZ
    ibegin = ibegin_edge3(ispecabs)
    iend = iend_edge3(ispecabs)
    ! exclude corners to make sure there is no contradiction on the normal
    if( codeabs_corner(3,ispecabs)) ibegin = 2
    if( codeabs_corner(4,ispecabs)) iend = NGLLX-1
    do i = ibegin,iend
      iglob = ibool(i,j,ispec)
      ! external velocity model
      if( assign_external_model ) then
        cpl = vpext(i,j,ispec)
        rhol = rhoext(i,j,ispec)
      endif
      xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
      zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
      jacobian1D = sqrt(xxi ** 2 + zxi ** 2)
      nx = - zxi / jacobian1D
      nz = + xxi / jacobian1D
      weight = jacobian1D * wxgll(i)

      call pnt_info_interpl_fluid(iglob,grad_pot_x_store,grad_pot_z_store,potential_dot_acoustic_store)

      ! adds absorbing boundary contribution
      potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - &
      ( nx*grad_pot_x_store + nz*grad_pot_z_store ) * weight/rhol - & !backgorund field
      ( potential_dot_acoustic(iglob) - potential_dot_acoustic_store ) * weight/cpl/rhol !scattered field

    enddo
  endif
end subroutine absorb_scatter_field_fluid

!the subroutines we use to calculate the values for GLL points along the absorbing boundary
!could be improved. Now we just find the nearest point and assign its value to the
!GLL point. If the located point is not that close to the target GLL point,
!the accuracy may not be enough
subroutine pnt_info_interpl_solid(iglob,veloc_x_store,veloc_y_store,veloc_z_store,&
                            tx_store,ty_store,tz_store)

!this subroutine is used to assign the information to the boundary points of local model. 
!Since the profile points we use are those exported from the local model,
!which is already the boundary points, here we just search the closest coordinateand assign the values directly.
  use specfem_par, only: nspec_bd_pnt_elastic,coord,x_final_bd_pnt_elastic,z_final_bd_pnt_elastic,&
                         vel_bd_pnt_elastic,trac_bd_pnt_elastic
                         
  
  implicit none
  include "constants.h"

  integer,intent(in) :: iglob
  real(kind=CUSTOM_REAL) :: veloc_x_store,veloc_y_store,veloc_z_store,tx_store,ty_store,tz_store
  double precision :: x,z
  double precision :: d_min,d_search
  integer :: point_locate, i 
  
  x = dble(coord(1,iglob))
  z = dble(coord(2,iglob))

  !initial distance
  d_min = (x-x_final_bd_pnt_elastic(1))**2 + (z-z_final_bd_pnt_elastic(1))**2
  point_locate = 1

  do i=2,nspec_bd_pnt_elastic
     d_search = (x-x_final_bd_pnt_elastic(i))**2 + (z-z_final_bd_pnt_elastic(i))**2 
     if(d_search < d_min )then
       d_min = d_search
       point_locate = i
     endif

  enddo

  veloc_x_store = vel_bd_pnt_elastic(1,point_locate)
  veloc_y_store = vel_bd_pnt_elastic(2,point_locate)
  veloc_z_store = vel_bd_pnt_elastic(3,point_locate)

  tx_store = trac_bd_pnt_elastic(1,point_locate)
  ty_store = trac_bd_pnt_elastic(2,point_locate)
  tz_store = trac_bd_pnt_elastic(3,point_locate)

end subroutine pnt_info_interpl_solid


subroutine pnt_info_interpl_fluid(iglob,grad_pot_x_store,grad_pot_z_store,potential_dot_acoustic_store)
  use specfem_par, only: nspec_bd_pnt_acoustic,coord,x_final_bd_pnt_acoustic,z_final_bd_pnt_acoustic,&
                         grad_pot_bd_pnt_acoustic,pot_dot_bd_pnt_acoustic

  implicit none
  include "constants.h"

  integer,intent(in) :: iglob
  real(kind=CUSTOM_REAL) :: grad_pot_x_store,grad_pot_z_store,potential_dot_acoustic_store
  double precision :: x,z
  double precision :: d_min,d_search
  integer :: point_locate, i 

  x = dble(coord(1,iglob))
  z = dble(coord(2,iglob))

  !initial distance
  d_min = (x-x_final_bd_pnt_acoustic(1))**2 + (z-z_final_bd_pnt_acoustic(1))**2
  point_locate = 1

  do i=2,nspec_bd_pnt_acoustic
     d_search = (x-x_final_bd_pnt_acoustic(i))**2 + (z-z_final_bd_pnt_acoustic(i))**2 
     if(d_search < d_min )then
       d_min = d_search
       point_locate = i
     endif

  enddo

  grad_pot_x_store = grad_pot_bd_pnt_acoustic(1,point_locate)
  grad_pot_z_store = grad_pot_bd_pnt_acoustic(2,point_locate)
  potential_dot_acoustic_store = pot_dot_bd_pnt_acoustic(point_locate)

end subroutine pnt_info_interpl_fluid
