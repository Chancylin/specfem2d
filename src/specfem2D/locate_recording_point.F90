!!!by lcx: this subroutine will be called by setup_recording_bd() in prepare_timerun_body.F90
!!!this subroutine will determine which elements the recording points
!belong to. And whether it is in fluid or solid region.
  subroutine locate_recording_point(ibool,coord,nspec,nglob,xigll,zigll, &
                      coorg,knods,ngnod,npgeo)

  use specfem_par, only: elastic,acoustic,p_sv,&
                         record_local_bkgd_boundary,npnt,&
                         ispec_selected_bd_pnt,nspec_bd_pnt_elastic,nspec_bd_pnt_acoustic,&
                         nspec_bd_elmt_elastic_pure,nspec_bd_elmt_acoustic_pure,&
                         ispec_bd_elmt_elastic_pure, ispec_bd_elmt_acoustic_pure,& 
                         hxi_bd_store,hgammar_bd_store,&
                         xi_bd_pnt,gamma_bd_pnt,&
                         bd_pnt_xval,bd_pnt_zval,&
                         nx_bd_pnt_elastic,nz_bd_pnt_elastic,&
                         nx_bd_pnt_acoustic,nz_bd_pnt_acoustic,&
                         nx_pnt,nz_pnt,fname,f_num,&
                         x_final_bd_pnt, z_final_bd_pnt,&
                         x_final_bd_pnt_elastic, z_final_bd_pnt_elastic,&
                         x_final_bd_pnt_acoustic, z_final_bd_pnt_acoustic,&
                         stress_bd_elastic,vel_bd_elastic,grad_pot_bd_acoustic,pot_dot_bd_acoustic,&
                         stress_bd_pnt_elastic,trac_bd_pnt_elastic,&
                         vel_bd_pnt_elastic,grad_pot_bd_pnt_acoustic,pot_dot_bd_pnt_acoustic,&
                         !para for reconst
                         record_local_boundary_reconst,bd_pnt_elmnt_num,side_type,side_type_elastic,&
                         ispec_bd_elmt_elastic,ispec_bd_elmt_acoustic,&
                         side_type_acoustic,ispec_bd_elmt_elastic_i,ispec_bd_elmt_elastic_j, &
                         ispec_bd_elmt_acoustic_i,ispec_bd_elmt_acoustic_j, &
                         !ispec_bd_elmt_elastic_pure_edge, ispec_bd_elmt_acoustic_pure_edge,&
                         !ispec_bd_elmt_elastic_pure_side, ispec_bd_elmt_acoustic_pure_side,&
                         !nspec_bd_elmt_elastic_pure_edge,nspec_bd_elmt_acoustic_pure_edge,&
                         trac_bd_pnt_elastic_reconst,trac_f,&
                         m_xx,m_xz,m_zz,m_zx,m_xx_reconst,m_xz_reconst,m_zz_reconst,m_zx_reconst,&
                         m_yx,m_yz,m_yx_reconst,m_yz_reconst,&
                         Grad_pot,grad_pot_x_reconst,grad_pot_z_reconst,Pot_x,Pot_z 

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
  integer, dimension(:), allocatable :: temp_bd_elmt_elastic,temp_bd_elmt_acoustic
  !character, dimension(:), allocatable :: temp_bd_elmt_elastic_side, temp_bd_elmt_acoustic_side
  !temperary variables for reading
  character(len=1) temp_side
  logical :: elastic_flag,acoustic_flag,corner_flag
  double precision :: box_t,box_b,box_l,box_r
 
 !para for recording the information to reconstruct the wavefield
  integer, dimension(:), allocatable :: bd_pnt_i,bd_pnt_j

  !!!geometry bound. Here we play the trick to locate the point in the exact 
  !!element we want. Otherwise, it could happen that the points at the edges
  !!are located in the neighbound element
  box_t = 10000.0
  box_b = -10000.0
  box_l = -63989.8
  box_r = -43989.8
  !!!

  !figure out the number of recording point firstly
  npnt = 0
  open(unit=1,file='DATA/boundary_points',iostat=ios,status='old',action='read')
  do while(ios == 0)
     read(1,"(a)",iostat=ios) dummystring
     if(ios == 0) npnt=npnt+1
  enddo
  close(1)
 
  print *,'total recording points: ',npnt
  
  allocate(ispec_selected_bd_pnt(npnt)) !this is the element index in global model
  allocate(xi_bd_pnt(npnt),gamma_bd_pnt(npnt))
  allocate(hxi_bd_store(npnt,NGLLX),hgammar_bd_store(npnt,NGLLZ))
  allocate(x_final_bd_pnt(npnt),z_final_bd_pnt(npnt))
  
  allocate(bd_pnt_elmnt_num(npnt))
  allocate(side_type(npnt))
  allocate(bd_pnt_i(npnt),bd_pnt_j(npnt))
  allocate(bd_pnt_xval(npnt),bd_pnt_zval(npnt))
  allocate(nx_pnt(npnt),nz_pnt(npnt))
  
  nspec_bd_pnt_elastic = 0
  nspec_bd_pnt_acoustic = 0

  113 format(i5.5,2x,2(i1.1,2x),L1,2x,L1,2x,L1,2x,A1,2x,A1,2x,4(es12.4,2x))
  open(unit=1,file='DATA/boundary_points',status='old',action='read')
  !loop over all points
  ! loop only on points inside the element
  ! exclude edges to ensure this point is not shared with other elements
  !comments: this could be improved by just comparing with the central point
  do ipnt=1,npnt

    distmin_squared = HUGEVAL

    !read(1,*) bd_pnt_xval(ipnt), bd_pnt_zval(ipnt), nx_pnt(ipnt), nz_pnt(ipnt)
    !note now when dealing with the reconstruting problem, you need the type of side
    !(i.e., L, R, T, B), which should be pre-known from the 'boundary_points' profile
    read(1,113)bd_pnt_elmnt_num(ipnt),bd_pnt_i(ipnt),bd_pnt_j(ipnt), elastic_flag,acoustic_flag, corner_flag,&
               side_type(ipnt), temp_side, bd_pnt_xval(ipnt), &
               bd_pnt_zval(ipnt), nx_pnt(ipnt), nz_pnt(ipnt)

    !the following loop is used to locate which element the recording point is in
    !And the algrithm is not perfect, since it could locate the corner recording point
    !in the neighbor element. I guess this is not a serious issue
    do ispec=1,nspec
       
       !!the if sentence is used to distinguish the edge point which could be shared
       !!by the neighbour elastic/fluid elements
       if( (elastic(ispec) .eqv. elastic_flag) .and. (acoustic(ispec) .eqv. acoustic_flag) ) then

      !!add the geometry bound in case it will pick the neighbor element which is
      !!out of the loca region.
      !!but somehow this is unnecessary, considering the traction and velocity should be
      !!continuous in the media which does not have sudden change of property
       !iglob = ibool(2,2,ispec)
       !if((dble(coord(1,iglob)) > box_l) .and. (dble(coord(1,iglob)) < box_r) &
       ! .and. dble(coord(2,iglob)) > box_b .and. dble(coord(2,iglob)) < box_t )then
      !!!!!!!!!!!!!!
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

               !else if (dist_squared == distmin_squared)then

                    ! if((dble(coord(1,iglob)) > box_l) .and. (dble(coord(1,iglob)) < box_r) .and. &
                    !     dble(coord(2,iglob)) > box_b .and. dble(coord(2,iglob)) < box_t )then

                       
             endif
           enddo
         enddo

       ! endif!geometry check 

      endif !check if the elastic/acoustic flag is correct

   enddo !end the loop for all elements
 
     if(elastic(ispec_selected_bd_pnt(ipnt))) nspec_bd_pnt_elastic=nspec_bd_pnt_elastic +1
     if(acoustic(ispec_selected_bd_pnt(ipnt))) nspec_bd_pnt_acoustic=nspec_bd_pnt_acoustic+1

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
 
  enddo !end loop for all recording points 

  close(1)
 
  !!this is a test by lcx: the check could be useful
  !!to check the the final coordinates located by global mesher are
  !!consistent with the import one.
  !fname = './OUTPUT_FILES/bg_record/final_pnts_profile'
  !open(unit=111,file=trim(fname),status='new',&
  !     action='write',iostat=ios)
  !if( ios /= 0 ) stop 'error saving acoustic point profile' 
  !
  !do ipnt=1,npnt
  !   write(111,110) ispec_selected_bd_pnt(ipnt),x_final_bd_pnt(ipnt),z_final_bd_pnt(ipnt)
  !enddo
  !close(111)
  !stop
  
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

  allocate(nx_bd_pnt_elastic(nspec_bd_pnt_elastic),nz_bd_pnt_elastic(nspec_bd_pnt_elastic))
  allocate(nx_bd_pnt_acoustic(nspec_bd_pnt_acoustic),nz_bd_pnt_acoustic(nspec_bd_pnt_acoustic))

  if ( record_local_bkgd_boundary ) then
     
     if ( nspec_bd_pnt_elastic /= 0 ) then
        i = 1
        do ipnt=1,npnt
          if(elastic(ispec_selected_bd_pnt(ipnt)))then
            ispec_bd_elmt_elastic(i) = ispec_selected_bd_pnt(ipnt)  !this is the element index in the global model
            
          !here we also record the coordinate for the point
            nx_bd_pnt_elastic(i) = nx_pnt(ipnt)
            nz_bd_pnt_elastic(i) = nz_pnt(ipnt)
            x_final_bd_pnt_elastic(i) = x_final_bd_pnt(ipnt)
            z_final_bd_pnt_elastic(i) = z_final_bd_pnt(ipnt)
            i = i+1
          endif
        
        enddo
     endif

     if ( nspec_bd_pnt_acoustic /= 0 ) then
        j = 1
        do ipnt=1,npnt
        
          if(acoustic(ispec_selected_bd_pnt(ipnt)))then
            ispec_bd_elmt_acoustic(j) = ispec_selected_bd_pnt(ipnt)

            nx_bd_pnt_acoustic(j) = nx_pnt(ipnt)
            nz_bd_pnt_acoustic(j) = nz_pnt(ipnt)
            x_final_bd_pnt_acoustic(j) = x_final_bd_pnt(ipnt)
            z_final_bd_pnt_acoustic(j) = z_final_bd_pnt(ipnt)

            j = j+1
          endif
        enddo
     endif

  endif

  !now ispec_bd_elmt_elastic() and ispec_bd_elmt_acoustic() could map
  ! to the real ispec, but we want to remove the duplicate
  !note: need to improve in future: we should consider the case
  !in which the elastic/fluid elment could be zero here !!!
  if( record_local_bkgd_boundary ) then  

     if ( nspec_bd_pnt_elastic /= 0 ) then

        allocate(temp_bd_elmt_elastic(nspec_bd_pnt_elastic))
        temp_bd_elmt_elastic(1) = ispec_bd_elmt_elastic(1)

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
            !print *,k,' th pure elastic,','global index as ',temp_bd_elmt_elastic(k)
        enddo loop1
        nspec_bd_elmt_elastic_pure = k
        !test
        print *,'total pure elastic elments are ', nspec_bd_elmt_elastic_pure

        !here k = nspec_bd_elmt_elastic_pure, kk = nspec_bd_elmt_acoustic_pure
        allocate(ispec_bd_elmt_elastic_pure(k))
        ispec_bd_elmt_elastic_pure(1:k) = temp_bd_elmt_elastic(1:k)
        deallocate(temp_bd_elmt_elastic)

     endif

     if ( nspec_bd_pnt_acoustic /= 0 ) then

        allocate(temp_bd_elmt_acoustic(nspec_bd_pnt_acoustic))
        temp_bd_elmt_acoustic(1) = ispec_bd_elmt_acoustic(1)

        kk=1
        loop2: do i=2,nspec_bd_pnt_acoustic
          do j=1,kk
             if ( ispec_bd_elmt_acoustic(i) == temp_bd_elmt_acoustic(j) ) then
                cycle loop2
             endif
          enddo
            kk = kk + 1
            temp_bd_elmt_acoustic(kk) = ispec_bd_elmt_acoustic(i)
            !print *,kk,' th pure elastic,','global index as ',temp_bd_elmt_acoustic(kk)
        enddo loop2
        nspec_bd_elmt_acoustic_pure = kk  
        !test
        print *,'total pure acoustic elments are ', nspec_bd_elmt_acoustic_pure

        allocate(ispec_bd_elmt_acoustic_pure(kk))
        ispec_bd_elmt_acoustic_pure(1:kk) = temp_bd_elmt_acoustic(1:kk)
        deallocate(temp_bd_elmt_acoustic)

     endif
     
      !deallocate(ispec_bd_elmt_elastic,ispec_bd_elmt_acoustic)!this will be used in 'record_local_boundary_reconst'
  endif

  if ( record_local_bkgd_boundary ) then
     !allocate array to store info for those elements in which the recording_bd_pnt locates
     if ( nspec_bd_pnt_elastic /= 0 ) then

        !elastic
        !allocate(trac_bd_elastic(3,NGLLX,NGLLZ,k))
        allocate(stress_bd_elastic(5,NGLLX,NGLLZ,k))
        allocate(vel_bd_elastic(3,NGLLX,NGLLZ,k))
        stress_bd_elastic = 0.0
        vel_bd_elastic = 0.0

        !allocate array to store info for these bd points
        !elastic
        allocate(stress_bd_pnt_elastic(5,nspec_bd_pnt_elastic))
        allocate(trac_bd_pnt_elastic(3,nspec_bd_pnt_elastic))
        allocate(vel_bd_pnt_elastic(3,nspec_bd_pnt_elastic))
        stress_bd_pnt_elastic = 0.0
        trac_bd_pnt_elastic = 0.0
        vel_bd_pnt_elastic = 0.0
     endif
   
     if ( nspec_bd_pnt_acoustic /= 0 ) then
        !acoustic
        allocate(grad_pot_bd_acoustic(2,NGLLX,NGLLZ,kk))
        allocate(pot_dot_bd_acoustic(NGLLX,NGLLZ,kk))
        grad_pot_bd_acoustic = 0.0
        pot_dot_bd_acoustic = 0.0

        !allocate array to store info for these bd points
        !acoustic
        allocate(grad_pot_bd_pnt_acoustic(2,nspec_bd_pnt_acoustic))
        allocate(pot_dot_bd_pnt_acoustic(nspec_bd_pnt_acoustic))
        grad_pot_bd_pnt_acoustic = 0.0
        pot_dot_bd_pnt_acoustic = 0.0
     endif

  endif

  if( record_local_boundary_reconst )then
  !when recording the information for wavefield reconstruction, we just need
  !the exact recording points provided from 'boundary_points' based on the local model
  !Thus, unnecessary to locate the points again 

     allocate(ispec_bd_elmt_elastic_i(nspec_bd_pnt_elastic))
     allocate(ispec_bd_elmt_elastic_j(nspec_bd_pnt_elastic))
     allocate(side_type_elastic(nspec_bd_pnt_elastic))

     allocate(ispec_bd_elmt_acoustic_i(nspec_bd_pnt_acoustic))
     allocate(ispec_bd_elmt_acoustic_j(nspec_bd_pnt_acoustic))
     allocate(side_type_acoustic(nspec_bd_pnt_acoustic))


     i = 1
     j = 1
     do ipnt=1,npnt
       if(elastic(ispec_selected_bd_pnt(ipnt)))then
         ispec_bd_elmt_elastic(i) = bd_pnt_elmnt_num(ipnt)  !!this is the element index in the local model
         
         ispec_bd_elmt_elastic_i(i) = bd_pnt_i(ipnt)
         ispec_bd_elmt_elastic_j(i) = bd_pnt_j(ipnt)
         side_type_elastic(i) = side_type(ipnt)
       !here we also record the coordinate for the point
         nx_bd_pnt_elastic(i) = nx_pnt(ipnt)
         nz_bd_pnt_elastic(i) = nz_pnt(ipnt)
         x_final_bd_pnt_elastic(i) = bd_pnt_xval(ipnt)
         z_final_bd_pnt_elastic(i) = bd_pnt_zval(ipnt)
         i = i+1
       endif
     
       if(acoustic(ispec_selected_bd_pnt(ipnt)))then
         ispec_bd_elmt_acoustic(j) = bd_pnt_elmnt_num(ipnt)
         ispec_bd_elmt_acoustic_i(j) = bd_pnt_i(ipnt)
         ispec_bd_elmt_acoustic_j(j) = bd_pnt_j(ipnt)
         side_type_acoustic(j) = side_type(ipnt)

         nx_bd_pnt_acoustic(j) = nx_pnt(ipnt)
         nz_bd_pnt_acoustic(j) = nz_pnt(ipnt)
         x_final_bd_pnt_acoustic(j) = bd_pnt_xval(ipnt)
         z_final_bd_pnt_acoustic(j) = bd_pnt_zval(ipnt)

         j = j+1
       endif
     enddo

!
!  !!we need to sort out what element edges we will operate. 
!  !!m11, m13, m31, m33. All of them will be array(nelement_edge_recording,5)
!  !!e.g., m11
!  !!for top edge: m11^i5
!  !!for bottom edge: m11^i1
!  !!for left edge: m11^1j
!  !!for right edge: m11^5j
!  !!the reason of this step is particularly due to the corner elements which hold two possible edges
! 
!     if ( nspec_bd_pnt_elastic /= 0 ) then
!
!        allocate(temp_bd_elmt_elastic(nspec_bd_pnt_elastic))
!        allocate(temp_bd_elmt_elastic_side(nspec_bd_pnt_elastic))
!        temp_bd_elmt_elastic(1) = ispec_bd_elmt_elastic(1)
!        temp_bd_elmt_elastic_side(1) = side_type_elastic(1) 
!        !delete duplicate
!        k=1
!        loop3: do i=2,nspec_bd_pnt_elastic
!          do j=1,k
!             if ( ispec_bd_elmt_elastic(i) == temp_bd_elmt_elastic(j) &
!                  .and. side_type_elastic(i) == temp_bd_elmt_elastic_side(j) ) then
!                cycle loop3
!             endif
!          enddo
!            k = k + 1
!            temp_bd_elmt_elastic(k) = ispec_bd_elmt_elastic(i)
!            temp_bd_elmt_elastic_side(k) = side_type_elastic(i)
!            !print *,k,' th pure elastic,','global index as ',temp_bd_elmt_elastic(k)
!        enddo loop3
!        nspec_bd_elmt_elastic_pure_edge = k
!        !test
!        print *,'total pure elastic elments (edge) are ', nspec_bd_elmt_elastic_pure_edge
!
!        allocate(ispec_bd_elmt_elastic_pure_edge(k))
!        allocate(ispec_bd_elmt_elastic_pure_side(k))
!        ispec_bd_elmt_elastic_pure_edge(1:k) = temp_bd_elmt_elastic(1:k)
!        ispec_bd_elmt_elastic_pure_side(1:k) = temp_bd_elmt_elastic_side(1:k)
!        deallocate(temp_bd_elmt_elastic)
!        deallocate(temp_bd_elmt_elastic_side)
!
!       if( .TRUE. ) then
!          open(117,file='./OUTPUT_FILES/reconst_record/pure_elastic_elments_edge',status='unknown',&
!               action='write',iostat=ios) 
!          if( ios /= 0 ) stop 'error saving elastic point profile'
!          write(117,*) nspec_bd_pnt_elastic 
!          write(117,*) nspec_bd_elmt_elastic_pure_edge
!          do kk = 1, nspec_bd_elmt_elastic_pure_edge
!             write(117,117) ispec_bd_elmt_elastic_pure_edge(kk), ispec_bd_elmt_elastic_pure_side(kk) 
!          enddo 
!
!          close(117)
!       endif
!  117  format(i5,2x,A1)
!     endif
!  
!
!     if ( nspec_bd_pnt_acoustic /= 0 ) then
!
!        allocate(temp_bd_elmt_acoustic(nspec_bd_pnt_acoustic))
!        allocate(temp_bd_elmt_acoustic_side(nspec_bd_pnt_acoustic))
!        temp_bd_elmt_acoustic(1) = ispec_bd_elmt_acoustic(1)
!        temp_bd_elmt_acoustic_side(1) = side_type_acoustic(1) 
!        !delete duplicate
!        k=1
!        loop4: do i=2,nspec_bd_pnt_acoustic
!          do j=1,k
!             if ( ispec_bd_elmt_acoustic(i) == temp_bd_elmt_acoustic(j) &
!                  .and. side_type_acoustic(i) == temp_bd_elmt_acoustic_side(j) ) then
!                cycle loop4
!             endif
!          enddo
!            k = k + 1
!            temp_bd_elmt_acoustic(k) = ispec_bd_elmt_acoustic(i)
!            temp_bd_elmt_acoustic_side(k) = side_type_acoustic(i)
!        enddo loop4
!        nspec_bd_elmt_acoustic_pure_edge = k
!        !test
!        print *,'total pure acoustic elments (edge) are ', nspec_bd_elmt_acoustic_pure_edge
!
!        allocate(ispec_bd_elmt_acoustic_pure_edge(k))
!        allocate(ispec_bd_elmt_acoustic_pure_side(k))
!        ispec_bd_elmt_acoustic_pure_edge(1:k) = temp_bd_elmt_acoustic(1:k)
!        ispec_bd_elmt_acoustic_pure_side(1:k) = temp_bd_elmt_acoustic_side(1:k)
!        deallocate(temp_bd_elmt_acoustic)
!        deallocate(temp_bd_elmt_acoustic_side)
!
!     endif
!
  endif


  !here we just export the coordinate of the recording point by separating them 
  !into elastic/acoustic
  if( record_local_bkgd_boundary ) then

    if ( nspec_bd_pnt_elastic /= 0 ) then
        fname = './OUTPUT_FILES/bg_record/elastic_pnts_profile'
        f_num = 111
        open(unit=f_num,file=trim(fname),status='new',&
             action='write',iostat=ios)
        if( ios /= 0 ) stop 'error saving elastic point profile'

        do i=1,nspec_bd_pnt_elastic
           !write(f_num,110) ispec_bd_elmt_elastic(i), x_final_bd_pnt_elastic(i), z_final_bd_pnt_elastic(i)
           write(f_num,110) ispec_bd_elmt_elastic(i), x_final_bd_pnt_elastic(i), z_final_bd_pnt_elastic(i)
                            ! nx_bd_pnt_elastic(i),nz_bd_pnt_elastic(i)
        enddo
        close(f_num)
        !stop 'here we see the normal vector. DO NOT foget to modify the writing format 100'
    endif

    if ( nspec_bd_pnt_acoustic /= 0 ) then
       fname = './OUTPUT_FILES/bg_record/acoustic_pnts_profile' 
       open(unit=f_num,file=trim(fname),status='new',&
            action='write',iostat=ios)
       if( ios /= 0 ) stop 'error saving acoustic point profile' 

       do i=1,nspec_bd_pnt_acoustic
          !print *,'i = ',i
          !print *,ispec_bd_elmt_acoustic(i)
          !print *,x_final_bd_pnt_acoustic(i)
          !print *,z_final_bd_pnt_acoustic(i)
          write(f_num,110) ispec_bd_elmt_acoustic(i), x_final_bd_pnt_acoustic(i), z_final_bd_pnt_acoustic(i)
       enddo
       close(f_num)
    endif
    
  endif

!!export the elastic/acoustic points profiles, which will be used to reconstruct the wavefield
  if( record_local_boundary_reconst ) then
    
     fname = './OUTPUT_FILES/reconst_record/elastic_pnts_profile'
     f_num = 111
     open(unit=f_num,file=trim(fname),status='new',&
          action='write',iostat=ios)
     if( ios /= 0 ) stop 'error saving elastic point profile'

     if ( nspec_bd_pnt_elastic /= 0 ) then
        do i=1,nspec_bd_pnt_elastic
           !write(f_num,110) ispec_bd_elmt_elastic(i), x_final_bd_pnt_elastic(i), z_final_bd_pnt_elastic(i)
           write(f_num,111) ispec_bd_elmt_elastic(i), ispec_bd_elmt_elastic_i(i), ispec_bd_elmt_elastic_j(i),&
                x_final_bd_pnt_elastic(i), z_final_bd_pnt_elastic(i)
           ! nx_bd_pnt_elastic(i),nz_bd_pnt_elastic(i)
        enddo
     endif
     close(f_num)

     fname = './OUTPUT_FILES/reconst_record/acoustic_pnts_profile' 
     open(unit=f_num,file=trim(fname),status='new',&
          action='write',iostat=ios)
     if( ios /= 0 ) stop 'error saving acoustic point profile' 

     if ( nspec_bd_pnt_acoustic /= 0 ) then
        do i=1,nspec_bd_pnt_acoustic
           write(f_num,111) ispec_bd_elmt_acoustic(i), ispec_bd_elmt_acoustic_i(i), ispec_bd_elmt_acoustic_j(i),& 
                x_final_bd_pnt_acoustic(i), z_final_bd_pnt_acoustic(i)
        enddo
     endif
     close(f_num)

  endif

  110 format(i5,2(es12.4,2x))!you may need to adjust the format depending on the precision
  111 format(i5,2x,i1,2x,i1,2x,2(es12.4,2x))!you may need to adjust the format depending on the precision

  deallocate(x_final_bd_pnt_elastic,z_final_bd_pnt_elastic)
  deallocate(x_final_bd_pnt_acoustic,z_final_bd_pnt_acoustic)

  



  if ( record_local_boundary_reconst ) then
     if ( nspec_bd_pnt_elastic /= 0 ) then
        !elastic
        allocate(trac_bd_pnt_elastic_reconst(3,nspec_bd_pnt_elastic))
        allocate(trac_f(3,nspec_bd_pnt_elastic))

        trac_bd_pnt_elastic_reconst = 0.0
        trac_f = 0.0
        
        !P-SV and SH cases will record different boundary information
        if( p_sv ) then
           allocate(m_xx(nspec_bd_pnt_elastic))
           allocate(m_xz(nspec_bd_pnt_elastic))
           allocate(m_zz(nspec_bd_pnt_elastic))
           allocate(m_zx(nspec_bd_pnt_elastic))
           allocate(m_xx_reconst(nspec_bd_pnt_elastic))
           allocate(m_xz_reconst(nspec_bd_pnt_elastic))
           allocate(m_zz_reconst(nspec_bd_pnt_elastic))
           allocate(m_zx_reconst(nspec_bd_pnt_elastic))
        else
           !should we also separate the allocation of memory for P-SV and SH cases
           allocate(m_yx(nspec_bd_pnt_elastic))
           allocate(m_yz(nspec_bd_pnt_elastic))
           allocate(m_yx_reconst(nspec_bd_pnt_elastic))
           allocate(m_yz_reconst(nspec_bd_pnt_elastic))
        endif

     endif

     if( nspec_bd_pnt_acoustic /= 0 ) then
        !acoustic
        if( p_sv ) then
           allocate(Grad_pot(nspec_bd_pnt_acoustic))
           allocate(grad_pot_x_reconst(nspec_bd_pnt_acoustic),grad_pot_z_reconst(nspec_bd_pnt_acoustic))
           allocate(Pot_x(nspec_bd_pnt_acoustic),Pot_z(nspec_bd_pnt_acoustic))
        endif
     endif

  endif

  deallocate(bd_pnt_elmnt_num)
  deallocate(bd_pnt_i,bd_pnt_j)

  deallocate(bd_pnt_xval,bd_pnt_zval)
  end subroutine locate_recording_point
