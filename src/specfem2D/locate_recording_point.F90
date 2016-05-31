!!!by lcx: this subroutine will be called by setup_recording_bd() in prepare_timerun_body.F90
!!!this subroutine will determine which elements the recording points
!belong to. And whether it is in fluid or solid region.
  subroutine locate_recording_point(ibool,coord,nspec,nglob,xigll,zigll, &
                      coorg,knods,ngnod,npgeo&
                      )

  use specfem_par, only: elastic,acoustic,&
                         npnt,ispec_selected_bd_pnt,nspec_bd_pnt_elastic,nspec_bd_pnt_acoustic,&
                        nspec_bd_elmt_elastic_pure,nspec_bd_elmt_acoustic_pure,&
                        ispec_bd_elmt_elastic_pure, ispec_bd_elmt_acoustic_pure,& 
                        hxi_bd_store,hgammar_bd_store,&
                        xi_bd_pnt,gamma_bd_pnt,&
                        bd_pnt_xval,bd_pnt_zval,&
                        nx_pnt,nz_pnt,fname,f_num,&
                        x_final_bd_pnt, z_final_bd_pnt,&
                        x_final_bd_pnt_elastic, z_final_bd_pnt_elastic,&
                        x_final_bd_pnt_acoustic, z_final_bd_pnt_acoustic,&
                        stress_bd_elastic,vel_bd_elastic,grad_pot_bd_acoustic,pot_dot_bd_acoustic,&
                        stress_bd_pnt_elastic,trac_bd_pnt_elastic,&
                        vel_bd_pnt_elastic,grad_pot_bd_pnt_acoustic,pot_dot_bd_pnt_acoustic

  implicit none
  include "constants.h"

  integer nspec,nglob,ngnod,npgeo

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

  double precision hxi_bd_pnt(NGLLX),hpxi_bd_pnt(NGLLX)
  double precision hgamma_bd_pnt(NGLLZ),hpgamma_bd_pnt(NGLLZ)

  double precision x,z,xix,xiz,gammax,gammaz,jacobian

 ! use dynamic allocation
  double precision distmin_squared, dist_squared

  integer  :: ios
  character(len=150) dummystring

  integer  :: k, kk
  !double precision, dimension(:), allocatable :: bd_pnt_xval,bd_pnt_zval 
  integer, dimension(:), allocatable :: ispec_bd_elmt_elastic,ispec_bd_elmt_acoustic
  integer, dimension(:), allocatable :: temp_bd_elmt_elastic,temp_bd_elmt_acoustic

  !figure out the number of recording point firstly
  npnt = 0
  open(unit=1,file='DATA/boundary_points',iostat=ios,status='old',action='read')
  do while(ios == 0)
     read(IIN,"(a)",iostat=ios) dummystring
     if(ios == 0) npnt=npnt+1
  enddo
  close(1)

  allocate(ispec_selected_bd_pnt(npnt))
  allocate(xi_bd_pnt(npnt),gamma_bd_pnt(npnt))
  allocate(x_final_bd_pnt(npnt),z_final_bd_pnt(npnt))

  allocate(bd_pnt_xval(npnt),bd_pnt_zval(npnt))
  allocate(nx_pnt(npnt),nz_pnt(npnt))
  
  nspec_bd_pnt_elastic = 0
  nspec_bd_pnt_acoustic = 0

  open(unit=1,file='DATA/boundary_points',status='old',action='read')
  !loop over all points
  ! loop only on points inside the element
  ! exclude edges to ensure this point is not shared with other elements
  !comments: this could be improved by just comparing with the central point
  do ipnt=1,npnt

    distmin_squared = HUGEVAL

    read(1,*) bd_pnt_xval(ipnt), bd_pnt_zval(ipnt), nx_pnt(ipnt), nz_pnt(ipnt)

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
 
     if(elastic(ispec_selected_bd_pnt(ipnt))) nspec_bd_pnt_elastic=nspec_bd_pnt_elastic +1
     if(acoustic(ispec_selected_bd_pnt(ipnt))) nspec_bd_pnt_acoustic=nspec_bd_pnt_acoustic+1

 
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
  
  !define and store lagrange interpolators at all the boundary recording points
  do ipnt=1,npnt
     call lagrange_any(xi_bd_pnt(ipnt),NGLLX,xigll,hxi_bd_pnt,hpxi_bd_pnt)
     call lagrange_any(gamma_bd_pnt(ipnt),NGLLZ,zigll,hgamma_bd_pnt,hpgamma_bd_pnt)
  
     hxi_bd_store(ipnt,:) = hxi_bd_pnt
     hgammar_bd_store(ipnt,:) = hgamma_bd_pnt 
     
  enddo


  allocate(ispec_bd_elmt_elastic(nspec_bd_pnt_elastic))
  allocate(ispec_bd_elmt_acoustic(nspec_bd_pnt_acoustic))

  allocate(x_final_bd_pnt_elastic(nspec_bd_pnt_elastic),z_final_bd_pnt_elastic(nspec_bd_pnt_elastic))
  allocate(x_final_bd_pnt_acoustic(nspec_bd_pnt_acoustic),z_final_bd_pnt_acoustic(nspec_bd_pnt_acoustic))

  i = 1
  j = 1
  do ipnt=1,npnt
    if(elastic(ispec_selected_bd_pnt(ipnt)))then
      ispec_bd_elmt_elastic(i) = ispec_selected_bd_pnt(ipnt)
    !here we also record the coordinate for the point
      x_final_bd_pnt_elastic(i) = x_final_bd_pnt(ipnt)
      z_final_bd_pnt_elastic(i) = z_final_bd_pnt(ipnt)
      i = i+1
    endif
  
    if(acoustic(ispec_selected_bd_pnt(ipnt)))then
      ispec_bd_elmt_acoustic(j) = ispec_selected_bd_pnt(ipnt)

      x_final_bd_pnt_acoustic(j) = x_final_bd_pnt(ipnt)
      z_final_bd_pnt_acoustic(j) = z_final_bd_pnt(ipnt)

      j = j+1
    endif
  enddo

  !now ispec_bd_elmt_elastic() and ispec_bd_elmt_acoustic() could map
  ! to the real ispec, but we want to remove the duplicate
  !note: need to improve in future: we should consider the case
  !in which the elastic/fluid elment could be zero here !!!
  allocate(temp_bd_elmt_elastic(nspec_bd_pnt_elastic))
  allocate(temp_bd_elmt_acoustic(nspec_bd_pnt_acoustic))
 
  temp_bd_elmt_elastic(1) = ispec_bd_elmt_elastic(1)
  temp_bd_elmt_acoustic(1) = ispec_bd_elmt_acoustic(1)
  
  !delete duplicate
  k=1
  loop1: do i=2,nspec_bd_pnt_elastic
    do j=1,k
       if ( ispec_bd_elmt_elastic(i) == temp_bd_elmt_elastic(j) ) then
          cycle loop1
       endif
    enddo
      k = k + 1
      temp_bd_elmt_elastic(k) = ispec_bd_elmt_elastic(i)
  enddo loop1
  nspec_bd_elmt_elastic_pure = k

  kk=1
  loop2: do i=2,nspec_bd_pnt_acoustic
    do j=1,kk
       if ( ispec_bd_elmt_acoustic(i) == temp_bd_elmt_acoustic(j) ) then
          cycle loop2
       endif
    enddo
      kk = kk + 1
      temp_bd_elmt_acoustic(kk) = ispec_bd_elmt_acoustic(i)
  enddo loop2
  nspec_bd_elmt_acoustic_pure = kk  

  !here we just export the coordinate of the recording point by separating them 
  !into elastic/acoustic
  fname = './OUTPUT_FILES/bg_record/elastic_pnts_profile'
  f_num = 111
  open(unit=f_num,file=trim(fname),status='new',&
       action='write',iostat=ios)
  if( ios /= 0 ) stop 'error saving elastic point profile'

  do i=1,nspec_bd_pnt_elastic
     write(f_num,110) x_final_bd_pnt_elastic(i),  z_final_bd_pnt_elastic(i)
  enddo
  close(f_num)
 
  fname = './OUTPUT_FILES/bg_record/acoustic_pnts_profile' 
  open(unit=f_num,file=trim(fname),status='new',&
       action='write',iostat=ios)
  if( ios /= 0 ) stop 'error saving acoustic point profile' 

  do i=1,nspec_bd_pnt_acoustic
     write(f_num,110) x_final_bd_pnt_acoustic(i), z_final_bd_pnt_acoustic(i)
  enddo
  close(f_num)

  110 format(3(es12.4,2x))!you may need to adjust the format depending on the precision
  deallocate(x_final_bd_pnt_elastic,z_final_bd_pnt_elastic)
  deallocate(x_final_bd_pnt_acoustic,z_final_bd_pnt_acoustic)

  !here k = nspec_bd_elmt_elastic_pure, kk = nspec_bd_elmt_acoustic_pure
  allocate(ispec_bd_elmt_elastic_pure(k))
  ispec_bd_elmt_elastic_pure(1:k) = temp_bd_elmt_elastic(1:k)
  
  allocate(ispec_bd_elmt_acoustic_pure(kk))
  ispec_bd_elmt_acoustic_pure(1:kk) = temp_bd_elmt_acoustic(1:kk)

  deallocate(temp_bd_elmt_elastic,temp_bd_elmt_acoustic)
  deallocate(ispec_bd_elmt_elastic,ispec_bd_elmt_acoustic)

  !allocate array to store info for those elements in which the recording_bd_pnt locates
  !elastic
  !allocate(trac_bd_elastic(3,NGLLX,NGLLZ,k))
  allocate(stress_bd_elastic(5,NGLLX,NGLLZ,k))
  allocate(vel_bd_elastic(3,NGLLX,NGLLZ,k))
  !acoustic
  allocate(grad_pot_bd_acoustic(2,NGLLX,NGLLZ,kk))
  allocate(pot_dot_bd_acoustic(NGLLX,NGLLZ,kk))

  !allocate array to store info for these bd points
  !elastic
  allocate(stress_bd_pnt_elastic(5,nspec_bd_pnt_elastic))
  allocate(trac_bd_pnt_elastic(3,nspec_bd_pnt_elastic))
  allocate(vel_bd_pnt_elastic(3,nspec_bd_pnt_elastic))
  !acoustic
  allocate(grad_pot_bd_pnt_acoustic(2,nspec_bd_pnt_acoustic))
  allocate(pot_dot_bd_pnt_acoustic(nspec_bd_pnt_acoustic))

  end subroutine locate_recording_point
