!!!by lcx: this subroutine will be called by setup_recording_bd() in prepare_timerun_body.F90
!!!this subroutine will determine which elements the recording points
!belong to. And whether it is in fluid or solid region.
  subroutine locate_recording_point(ibool,coord,nspec,nglob,xigll,zigll, &
                      coorg,knods,ngnod,npgeo)


    use mpi
    
    use specfem_par, only: myrank,nproc,elastic,acoustic,p_sv,&
                           record_local_bkgd_boundary,npnt,npnt_local,num_pnt_elastic,num_pnt_acoustic,&
                           bg_record_elastic, bg_record_acoustic,&
                           nspec_bd_pnt_elastic_clt, nspec_bd_pnt_acoustic_clt,&
                           ispec_selected_bd_pnt,nspec_bd_pnt_elastic,nspec_bd_pnt_acoustic,&
                           nspec_bd_elmt_elastic_pure,nspec_bd_elmt_acoustic_pure,&
                           ispec_bd_elmt_elastic_pure, ispec_bd_elmt_acoustic_pure,& 
                           hxi_bd_store,hgammar_bd_store,&
                           xi_bd_pnt,gamma_bd_pnt,&
                           ! bd_pnt_xval,bd_pnt_zval,&
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
                           Grad_pot,grad_pot_x_reconst,grad_pot_z_reconst,Pot_x,Pot_z,& 
                           temp_record_elastic,temp_record_acoustic

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
  double precision distmin_squared, dist_squared, dist_glob_squared

  integer  :: ios
  character(len=150) dummystring

  integer  :: k, kk
  !double precision, dimension(:), allocatable :: bd_pnt_xval,bd_pnt_zval 
  integer, dimension(:), allocatable :: temp_bd_elmt_elastic,temp_bd_elmt_acoustic
  !character, dimension(:), allocatable :: temp_bd_elmt_elastic_side, temp_bd_elmt_acoustic_side
  !temperary variables for reading
  character(len=1) temp_side
  logical, dimension(:), allocatable :: in_element
  double precision :: box_t,box_b,box_l,box_r

  !dynamic array for wavefield record in global model
  integer, dimension(:), allocatable :: bd_pnt_elmnt_num_total
  integer, dimension(:), allocatable :: bd_pnt_i_total, bd_pnt_j_total
  character, dimension(:), allocatable :: side_type_total
  integer, dimension(:), allocatable :: ispec_selected_bd_pnt_total
  logical, dimension(:), allocatable :: elastic_flag_total, acoustic_flag_total,corner_flag_total
  logical,dimension(:), allocatable :: elastic_flag,acoustic_flag!,corner_flag
  double precision, dimension(:), allocatable :: bd_pnt_xval_total, bd_pnt_zval_total
  double precision, dimension(:), allocatable :: nx_pnt_total, nz_pnt_total
  double precision, dimension(:), allocatable ::xi_bd_pnt_total, gamma_bd_pnt_total
  double precision, dimension(:), allocatable :: x_final_bd_pnt_total, z_final_bd_pnt_total

 
 !para for recording the information to reconstruct the wavefield
  integer, dimension(:), allocatable :: bd_pnt_i, bd_pnt_j
  integer, dimension(:), allocatable :: bd_pnt_i_elastic, bd_pnt_j_elastic
  integer, dimension(:), allocatable :: bd_pnt_i_acoustic, bd_pnt_j_acoustic

  !para for mpi
  integer :: ierror
  integer :: is_proc_recording_pnt
  integer :: nb_proc_recording_pnt
  integer, dimension(1:nproc) :: allgather_is_proc_recording_pnt
  integer, dimension(1) :: locate_is_proc_recording_pnt

  print *,myrank,'this has ran' 

  
  !!!geometry bound. Here we play the trick to locate the point in the exact 
  !!element we want. Otherwise, it could happen that the points at the edges
  !!are located in the neighbound element
  box_t = 10000.0
  box_b = -10000.0
  box_l = -63989.8
  box_r = -43989.8
  !!!

  ! print *, nspec, ' in rank ', myrank
  !figure out the number of recording point firstly
  npnt = 0
  if( record_local_bkgd_boundary ) then

     open(unit=1,file='DATA/boundary_points_total',iostat=ios,status='old',action='read')
     do while(ios == 0)
        read(1,"(a)",iostat=ios) dummystring
        if(ios == 0) npnt=npnt+1
     enddo
     close(1)

     print *,'total recording points from profile: ',npnt

  else if( record_local_boundary_reconst ) then

     write(fname,"('./DATA/boundary_points',i5.5)") myrank
     open(unit=1,file=trim(fname),iostat=ios,status='old',action='read')
     do while(ios == 0)
        read(1,"(a)",iostat=ios) dummystring
        if(ios == 0) npnt=npnt+1
     enddo
     close(1)

     print *,'total recording points from profile: ',npnt, ' for rank', myrank

  else
     stop 'flag is incorrect for the usage of locate_recording_points()'
     
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! #ifdef USE_MPI
  !   allocate(found_element(npnt))
  !   found_element = .false.
  !   call MPI_BARRIER(MPI_COMM_WORLD, ierror)
  ! #endif
  allocate(in_element(npnt))
  allocate(elastic_flag_total(npnt),acoustic_flag_total(npnt),corner_flag_total(npnt))
  allocate(ispec_selected_bd_pnt_total(npnt)) !this is the element index in global model
  allocate(xi_bd_pnt_total(npnt),gamma_bd_pnt_total(npnt))
  allocate(x_final_bd_pnt_total(npnt),z_final_bd_pnt_total(npnt))
  
  allocate(bd_pnt_elmnt_num_total(npnt))
  allocate(side_type_total(npnt))
  allocate(bd_pnt_i_total(npnt),bd_pnt_j_total(npnt))
  allocate(bd_pnt_xval_total(npnt),bd_pnt_zval_total(npnt))
  allocate(nx_pnt_total(npnt),nz_pnt_total(npnt))

  ! !build the booking
  ! allocate(booking_record = )
  ! booking_record = 
  
  in_element = .true.
  ispec_selected_bd_pnt_total = 0 

  nspec_bd_pnt_elastic = 0
  nspec_bd_pnt_acoustic = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
113 format(i5.5,2x,2(i1.1,2x),L1,2x,L1,2x,L1,2x,A1,2x,A1,2x,4(es12.4,2x))
  
  if( record_local_bkgd_boundary ) then

     open(unit=1,file='DATA/boundary_points_total',status='old',action='read')
     ! loop over all points
     ! loop only on points inside the element
     ! exclude edges to ensure this point is not shared with other elements
     !comments: this could be improved by just comparing with the central point
     ipnt_locate: do ipnt=1,npnt

        distmin_squared = HUGEVAL

        is_proc_recording_pnt = 0
        !read(1,*) bd_pnt_xval(ipnt), bd_pnt_zval(ipnt), nx_pnt(ipnt), nz_pnt(ipnt)
        !note now when dealing with the reconstruting problem, you need the type of side
        !(i.e., L, R, T, B), which should be pre-known from the 'boundary_points' profile

        !lcx: you may want to use derived datatype or structure to process the data.
        !and then use map/union to quickly acess to the dataset
        read(1,113)bd_pnt_elmnt_num_total(ipnt),bd_pnt_i_total(ipnt),bd_pnt_j_total(ipnt),&
             elastic_flag_total(ipnt),acoustic_flag_total(ipnt),&
             corner_flag_total(ipnt),side_type_total(ipnt), temp_side,&
             bd_pnt_xval_total(ipnt), bd_pnt_zval_total(ipnt),&
             nx_pnt_total(ipnt), nz_pnt_total(ipnt)


        ! !recording point already found in other processors
        ! if( found_element(ipnt) ) cycle ipnt_locate

        !the following loop is used to locate which element the recording point is in
        !And the algrithm is not perfect, since it could locate the corner recording point
        !in the neighbor element. I guess this is not a serious issue
        do ispec=1,nspec

           !!the if sentence is used to distinguish the edge point which could be shared
           !!by the neighbour elastic/fluid elements
           if( (elastic(ispec) .eqv. elastic_flag_total(ipnt) ) .and. (acoustic(ispec) .eqv. acoustic_flag_total(ipnt) ) ) then

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
                    dist_squared = (bd_pnt_xval_total(ipnt)-dble(coord(1,iglob)))**2 + &
                         (bd_pnt_zval_total(ipnt)-dble(coord(2,iglob)))**2

                    if(dist_squared < distmin_squared) then
                       distmin_squared = dist_squared
                       ispec_selected_bd_pnt_total(ipnt) = ispec
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! we add the global communication to solve the allocating problem when
        ! the recording point is at the interface of different partitions
#ifdef USE_MPI
        ! global minimum distance computed over all processes
        call MPI_ALLREDUCE (distmin_squared, dist_glob_squared, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)
#else
        dist_glob_squared = distmin_squared
#endif

        ! check if this process contains the source
        if ( abs(sqrt(dist_glob_squared) - sqrt(distmin_squared)) < TINYVAL ) is_proc_recording_pnt = 1
        
#ifdef USE_MPI
        ! determining the number of processes that contain the recording point
        ! (useful when the source is located on an interface)
        call MPI_ALLREDUCE (is_proc_recording_pnt, nb_proc_recording_pnt, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
#else
        nb_proc_recording_pnt = is_proc_recording_pnt
#endif

#ifdef USE_MPI
        ! when several processes contain the source, we elect one of them (minimum rank).
        if ( nb_proc_recording_pnt > 1 ) then

           call MPI_ALLGATHER(is_proc_recording_pnt, 1, MPI_INTEGER, allgather_is_proc_recording_pnt(1), &
                1, MPI_INTEGER, MPI_COMM_WORLD, ierror)
           
           locate_is_proc_recording_pnt = maxloc(allgather_is_proc_recording_pnt) - 1

           if ( myrank /= locate_is_proc_recording_pnt(1) ) then
              is_proc_recording_pnt = 0
           endif
           nb_proc_recording_pnt = 1

        endif

#endif

        if( is_proc_recording_pnt == 0 ) then
           in_element(ipnt) = .false.
           cycle ipnt_locate
        endif
        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        !if ispec_selected_bd_pnt_total(ipnt) is zero, which means the flag condition
        !cannot be met when the previous loop try to locate it. That means the recording
        !point is not in the current partition. Then just jump to next recording point
        if( ispec_selected_bd_pnt_total(ipnt) == 0 )then
           in_element(ipnt) = .false.
           cycle ipnt_locate
        endif

        ! print *, ipnt, ' and ', ispec_selected_bd_pnt_total(ipnt), ' from rank ', myrank
        if(elastic(ispec_selected_bd_pnt_total(ipnt))) nspec_bd_pnt_elastic=nspec_bd_pnt_elastic +1
        if(acoustic(ispec_selected_bd_pnt_total(ipnt))) nspec_bd_pnt_acoustic=nspec_bd_pnt_acoustic+1

        !find the best (xi, gamma) for each recording point
        !start using intial guess in xi and gamma
        xi = xigll(ix_initial_guess)
        gamma = zigll(iz_initial_guess)

        do iter_loop = 1,NUM_ITER
           ! compute coordinates of the new point and derivatives dxi/dx, dxi/dz
           call recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian, &
                coorg,knods,ispec_selected_bd_pnt_total(ipnt),ngnod,nspec,npgeo, &
                .true.)

           ! compute distance to target location
           dx = - (x - bd_pnt_xval_total(ipnt))
           dz = - (z - bd_pnt_zval_total(ipnt))

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
           !lcx: we can use this as criteria to check whether the recording points is in the element

           ! if (xi > 1.10d0) xi = 1.10d0
           ! if (xi < -1.10d0) xi = -1.10d0
           ! if (gamma > 1.10d0) gamma = 1.10d0
           ! if (gamma < -1.10d0) gamma = -1.10d0


           ! end of non linear iterations
        enddo

! !since now we locate the point according to the global_min_dist, the if-sentence here is unnecessary any more.
!         !check whether the recording points is in the located element
!         if( xi > 1.10d0 .or. xi < -1.10d0 .or. gamma > 1.10d0 .or. gamma < -1.10d0)then
!            in_element(ipnt) = .false.

!            !since we already use the elastic/acoustic_flag, the calculation of
!            !elastic/acoustic elements here is unnecessary
!            if( elastic(ispec_selected_bd_pnt_total(ipnt)) )then
!               nspec_bd_pnt_elastic = nspec_bd_pnt_elastic - 1
!            else if( acoustic(ispec_selected_bd_pnt_total(ipnt)) )then
!               nspec_bd_pnt_acoustic = nspec_bd_pnt_acoustic - 1
!            endif

!         endif

        ! compute final coordinates of point found


        !by lcx: I think the basic idea here is we want to find the best 
        !(xi,gammar) for the point (x_bd_pnt,z_bd_pnt). Hence we need to 
        ! recompute the jacobian for several iterations. Once a point
        ! which is closed enough is found in (xi,gammar) (nature) coordinate
        ! we will reproject it back to the real coordinate

        call recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian, &
             coorg,knods,ispec_selected_bd_pnt_total(ipnt),ngnod,nspec,npgeo, &
             .true.)

        ! store xi,gamma of point found
        xi_bd_pnt_total(ipnt) = xi
        gamma_bd_pnt_total(ipnt) = gamma


        x_final_bd_pnt_total(ipnt) = x
        z_final_bd_pnt_total(ipnt) = z


     enddo ipnt_locate !end loop for all recording points 

     close(1)

     
  else if( record_local_boundary_reconst ) then

     write(fname,"('./DATA/boundary_points',i5.5)") myrank
     open(unit=1,file=trim(fname),iostat=ios,status='old',action='read')
     
     do ipnt=1,npnt
        
        read(1,113)bd_pnt_elmnt_num_total(ipnt),bd_pnt_i_total(ipnt),bd_pnt_j_total(ipnt),&
             elastic_flag_total(ipnt),acoustic_flag_total(ipnt),&
             corner_flag_total(ipnt),side_type_total(ipnt), temp_side,&
             bd_pnt_xval_total(ipnt), bd_pnt_zval_total(ipnt),&
             nx_pnt_total(ipnt), nz_pnt_total(ipnt)

     enddo

     close(1)

     x_final_bd_pnt_total = bd_pnt_xval_total
     z_final_bd_pnt_total = bd_pnt_zval_total

  endif

 
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

  !figure out whether the points are in the current partition
  ! do ipnt=1,npnt
  !    if( xi > 1.10d0 .or. xi < -1.10d0 .or. gamma > 1.10d0 .or. gamma < -1.10d0)then
  !       in_element(ipnt) = .false.

  !       if( elastic(ispec_selected_bd_pnt(ipnt)) )then
  !          nspec_bd_pnt_elastic = nspec_bd_pnt_elastic - 1
  !       else if( acoustic(ispec_selected_bd_pnt(ipnt)) )then
  !          nspec_bd_pnt_acoustic = nspec_bd_pnt_acoustic - 1
  !       endif

  !    endif
  ! enddo

  !calculate the number of recording points located in this partition
  !actually, for recording in local model, the npnt_local is just the npnt read from the profile, i.e., npnt_local = npnt
  npnt_local = count(in_element)

  !for test purpose
  if( record_local_boundary_reconst ) then
     if(npnt_local /= npnt ) stop 'something wrong'
  endif

  print *, 'recording points in partition', myrank, ' is ', npnt_local

  
  !create the arrays for recording points in this partition allocate(xi_bd_pnt_total(npnt),gamma_bd_pnt_total(npnt))

  !for both in global and local model
  allocate(ispec_selected_bd_pnt(npnt_local))
  allocate(elastic_flag(npnt_local),acoustic_flag(npnt_local))
  allocate(bd_pnt_elmnt_num(npnt_local))
  allocate(nx_pnt(npnt_local),nz_pnt(npnt_local))

  !for recording in global model
  allocate(xi_bd_pnt(npnt_local),gamma_bd_pnt(npnt_local))
  allocate(x_final_bd_pnt(npnt_local),z_final_bd_pnt(npnt_local))

  !for recording in local model
  allocate(side_type(npnt_local))
  allocate(bd_pnt_i(npnt_local),bd_pnt_j(npnt_local))

  !unnecessary
  ! allocate(bd_pnt_xval(npnt_local),bd_pnt_zval(npnt_local))

  ispec_selected_bd_pnt = pack(ispec_selected_bd_pnt_total,in_element)
  elastic_flag          = pack(elastic_flag_total,in_element)
  acoustic_flag         = pack(acoustic_flag_total,in_element)
  side_type             = pack(side_type_total,in_element) 
  bd_pnt_elmnt_num      = pack(bd_pnt_elmnt_num_total,in_element) 
  bd_pnt_i              = pack(bd_pnt_i_total,in_element) 
  bd_pnt_j              = pack(bd_pnt_j_total,in_element) 
  ! bd_pnt_xval           = pack(bd_pnt_xval_total,in_element) 
  ! bd_pnt_zval           = pack(bd_pnt_zval_total,in_element) 
  nx_pnt                = pack(nx_pnt_total,in_element) 
  nz_pnt                = pack(nz_pnt_total,in_element) 
  x_final_bd_pnt        = pack(x_final_bd_pnt_total,in_element)
  z_final_bd_pnt        = pack(z_final_bd_pnt_total,in_element)
  
  !define and store lagrange interpolators at all the boundary recording points
  if( record_local_bkgd_boundary ) then

     xi_bd_pnt             = pack(xi_bd_pnt_total,in_element)
     gamma_bd_pnt          = pack(gamma_bd_pnt_total,in_element)

     allocate(hxi_bd_store(npnt_local,NGLLX),hgammar_bd_store(npnt_local,NGLLZ))

     do ipnt=1,npnt_local
        call lagrange_any(xi_bd_pnt(ipnt),NGLLX,xigll,hxi_bd_pnt,hpxi_bd_pnt)
        call lagrange_any(gamma_bd_pnt(ipnt),NGLLZ,zigll,hgamma_bd_pnt,hpgamma_bd_pnt)

        hxi_bd_store(ipnt,:) = hxi_bd_pnt
        hgammar_bd_store(ipnt,:) = hgamma_bd_pnt 

     enddo
  endif

  nspec_bd_pnt_elastic  = count(elastic_flag)
  nspec_bd_pnt_acoustic = count(acoustic_flag)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !this should be written as a subroutine in future
#ifdef USE_MPI
  
  !this calculates the offset and prepare for the mpi writing
  !call the subroutine to figure out the accumulative points number in partition
  call calculate_accumulative_pnts(myrank,nspec_bd_pnt_elastic,num_pnt_elastic)
  call calculate_accumulative_pnts(myrank,nspec_bd_pnt_acoustic,num_pnt_acoustic)

  !need to calculate the total nspec_bd_pnt_elastic and nspec_bd_pnt_acoustic of all the partitions
  call MPI_ALLREDUCE(nspec_bd_pnt_elastic, nspec_bd_pnt_elastic_clt, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
  call MPI_ALLREDUCE(nspec_bd_pnt_acoustic, nspec_bd_pnt_acoustic_clt, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
  
  
  if( myrank == 0 ) then
     print *, 'offset: nspec_bd_pnt_elastic_clt = ', nspec_bd_pnt_elastic_clt, &
          ' from rank ', myrank
     print *, 'offset: nspec_bd_pnt_acoustic_clt = ', nspec_bd_pnt_acoustic_clt, &
          ' from rank ', myrank
  endif
  ! !for test
  ! print *, 'offset: nspec_bd_pnt_elastic = ', nspec_bd_pnt_elastic, &
  !      ' num_pnt_elastic = ', num_pnt_elastic, ' from rank ', myrank
  ! print *, 'offset: nspec_bd_pnt_acoustic = ', nspec_bd_pnt_acoustic, &
  !      ' num_pnt_acoustic = ', num_pnt_acoustic, ' from rank ', myrank
  ! ! stop 'bad testing'
  
  !the reason doing this is: some processors may have not any recording
  !points. so we need to exclude them when we do the writting
  call build_commu_bg_record(nspec_bd_pnt_elastic, nspec_bd_pnt_acoustic,  bg_record_elastic, bg_record_acoustic)

#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !this is the target elements where the recording points locate in
  !we try to seperate them into elastic and acoustic these are the necessary
  !information to calculate the traction by lagrange interpolation in the target
  !elements in global simulation
  allocate(ispec_bd_elmt_elastic(nspec_bd_pnt_elastic)) 
  allocate(ispec_bd_elmt_acoustic(nspec_bd_pnt_acoustic))
  allocate(bd_pnt_i_elastic(nspec_bd_pnt_elastic))
  allocate(bd_pnt_j_elastic(nspec_bd_pnt_elastic))
  allocate(bd_pnt_i_acoustic(nspec_bd_pnt_acoustic))
  allocate(bd_pnt_j_acoustic(nspec_bd_pnt_acoustic))
  
  allocate(x_final_bd_pnt_elastic(nspec_bd_pnt_elastic),z_final_bd_pnt_elastic(nspec_bd_pnt_elastic))
  allocate(x_final_bd_pnt_acoustic(nspec_bd_pnt_acoustic),z_final_bd_pnt_acoustic(nspec_bd_pnt_acoustic))

  allocate(nx_bd_pnt_elastic(nspec_bd_pnt_elastic),nz_bd_pnt_elastic(nspec_bd_pnt_elastic))
  allocate(nx_bd_pnt_acoustic(nspec_bd_pnt_acoustic),nz_bd_pnt_acoustic(nspec_bd_pnt_acoustic))


  if ( record_local_bkgd_boundary ) then
     
     if ( nspec_bd_pnt_elastic /= 0 ) then

        nx_bd_pnt_elastic      = pack(nx_pnt,elastic_flag)
        nz_bd_pnt_elastic      = pack(nz_pnt,elastic_flag)
        x_final_bd_pnt_elastic = pack(x_final_bd_pnt,elastic_flag)
        z_final_bd_pnt_elastic = pack(z_final_bd_pnt,elastic_flag)
        ispec_bd_elmt_elastic  = pack(ispec_selected_bd_pnt,elastic_flag)

        ! old algorithm, no using array operation
        ! i = 1
        ! do ipnt=1,npnt_local
        !   if(elastic(ispec_selected_bd_pnt(ipnt)))then
        !     ispec_bd_elmt_elastic(i) = ispec_selected_bd_pnt(ipnt)  !this is the element index in the global model
            
        !   !here we also record the coordinate for the point
        !     nx_bd_pnt_elastic(i) = nx_pnt(ipnt)
        !     nz_bd_pnt_elastic(i) = nz_pnt(ipnt)
        !     x_final_bd_pnt_elastic(i) = x_final_bd_pnt(ipnt)
        !     z_final_bd_pnt_elastic(i) = z_final_bd_pnt(ipnt)
        !     i = i+1
        !   endif
        
        ! enddo
     endif

     if ( nspec_bd_pnt_acoustic /= 0 ) then
        nx_bd_pnt_acoustic      = pack(nx_pnt,acoustic_flag)
        nz_bd_pnt_acoustic      = pack(nz_pnt,acoustic_flag)
        x_final_bd_pnt_acoustic = pack(x_final_bd_pnt,acoustic_flag)
        z_final_bd_pnt_acoustic = pack(z_final_bd_pnt,acoustic_flag)
        ispec_bd_elmt_acoustic  = pack(ispec_selected_bd_pnt,acoustic_flag)

        ! old algorithm, no using array operation
        ! j = 1
        ! do ipnt=1,npnt
        
        !   if(acoustic(ispec_selected_bd_pnt(ipnt)))then
        !     ispec_bd_elmt_acoustic(j) = ispec_selected_bd_pnt(ipnt)

        !     nx_bd_pnt_acoustic(j) = nx_pnt(ipnt)
        !     nz_bd_pnt_acoustic(j) = nz_pnt(ipnt)
        !     x_final_bd_pnt_acoustic(j) = x_final_bd_pnt(ipnt)
        !     z_final_bd_pnt_acoustic(j) = z_final_bd_pnt(ipnt)

        !     j = j+1
        !   endif
        ! enddo
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

  !allocate array to store info for those elements in which the recording_bd_pnt locates
  if ( record_local_bkgd_boundary ) then

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

  !here we just export the element index in global model partition, the
  !coordinate of the recording point by separating them
  !into elastic/acoustic
  if( record_local_bkgd_boundary ) then

#ifdef USE_MPI
     write(fname,"('./OUTPUT_FILES/bg_record/elastic_pnts_profile',i5.5)") myrank
#else
     fname = './OUTPUT_FILES/bg_record/elastic_pnts_profile' 
#endif
     f_num = 111
     open(unit=f_num,file=trim(fname),status='new',&
          action='write',iostat=ios)
     if( ios /= 0 ) stop 'error saving elastic point profile'

     if ( nspec_bd_pnt_elastic /= 0 ) then
        do i=1,nspec_bd_pnt_elastic
           !write(f_num,110) ispec_bd_elmt_elastic(i), x_final_bd_pnt_elastic(i), z_final_bd_pnt_elastic(i)
           write(f_num,110) ispec_bd_elmt_elastic(i), x_final_bd_pnt_elastic(i), z_final_bd_pnt_elastic(i)
           ! nx_bd_pnt_elastic(i),nz_bd_pnt_elastic(i)
        enddo
     endif

     close(f_num)
     !stop 'here we see the normal vector. DO NOT foget to modify the writing format 100'

#ifdef USE_MPI
     write(fname,"('./OUTPUT_FILES/bg_record/acoustic_pnts_profile',i5.5)") myrank
#else
     fname = './OUTPUT_FILES/bg_record/acoustic_pnts_profile' 
#endif

     open(unit=f_num,file=trim(fname),status='new',&
          action='write',iostat=ios)
     if( ios /= 0 ) stop 'error saving acoustic point profile' 

     if ( nspec_bd_pnt_acoustic /= 0 ) then
        do i=1,nspec_bd_pnt_acoustic
           !print *,'i = ',i
           !print *,ispec_bd_elmt_acoustic(i)
           !print *,x_final_bd_pnt_acoustic(i)
           !print *,z_final_bd_pnt_acoustic(i)
           write(f_num,110) ispec_bd_elmt_acoustic(i), x_final_bd_pnt_acoustic(i), z_final_bd_pnt_acoustic(i)
        enddo
     endif

     close(f_num)

  endif

  
  !When recording the information for wavefield reconstruction, we just need
  !the exact recording points provided from 'boundary_points' based on the local model
  !Thus, unnecessary to locate the points again 
  if( record_local_boundary_reconst )then

     if( nspec_bd_pnt_elastic /= 0 ) then
        
        nx_bd_pnt_elastic      = pack(nx_pnt,elastic_flag)
        nz_bd_pnt_elastic      = pack(nz_pnt,elastic_flag)
        x_final_bd_pnt_elastic = pack(x_final_bd_pnt,elastic_flag)
        z_final_bd_pnt_elastic = pack(z_final_bd_pnt,elastic_flag)
        
        bd_pnt_i_elastic       = pack(bd_pnt_i,elastic_flag)
        bd_pnt_j_elastic       = pack(bd_pnt_j,elastic_flag)

        allocate(ispec_bd_elmt_elastic_i(nspec_bd_pnt_elastic))
        allocate(ispec_bd_elmt_elastic_j(nspec_bd_pnt_elastic))
        allocate(side_type_elastic(nspec_bd_pnt_elastic))

        ispec_bd_elmt_elastic   = pack(bd_pnt_elmnt_num,elastic_flag)
        ispec_bd_elmt_elastic_i = bd_pnt_i_elastic
        ispec_bd_elmt_elastic_j = bd_pnt_j_elastic
        side_type_elastic       = pack(side_type,elastic_flag)

     endif

     if( nspec_bd_pnt_acoustic /= 0 ) then

        nx_bd_pnt_acoustic      = pack(nx_pnt,acoustic_flag)
        nz_bd_pnt_acoustic      = pack(nz_pnt,acoustic_flag)
        x_final_bd_pnt_acoustic = pack(x_final_bd_pnt,acoustic_flag)
        z_final_bd_pnt_acoustic = pack(z_final_bd_pnt,acoustic_flag)
        
        bd_pnt_i_acoustic       = pack(bd_pnt_i,acoustic_flag)
        bd_pnt_j_acoustic       = pack(bd_pnt_j,acoustic_flag)

        allocate(ispec_bd_elmt_acoustic_i(nspec_bd_pnt_acoustic))
        allocate(ispec_bd_elmt_acoustic_j(nspec_bd_pnt_acoustic))
        allocate(side_type_acoustic(nspec_bd_pnt_acoustic))

        ispec_bd_elmt_acoustic   = pack(bd_pnt_elmnt_num,acoustic_flag)
        ispec_bd_elmt_acoustic_i = bd_pnt_i_acoustic
        ispec_bd_elmt_acoustic_j = bd_pnt_j_acoustic
        side_type_acoustic       = pack(side_type,acoustic_flag)

     endif
     
     

     
     ! i = 1
     ! j = 1
     ! do ipnt=1,npnt
     !   if(elastic(ispec_selected_bd_pnt(ipnt)))then
     !     ispec_bd_elmt_elastic(i) = bd_pnt_elmnt_num(ipnt)  !!this is the element index in the local model
         
     !     ispec_bd_elmt_elastic_i(i) = bd_pnt_i(ipnt)
     !     ispec_bd_elmt_elastic_j(i) = bd_pnt_j(ipnt)
     !     side_type_elastic(i) = side_type(ipnt)
     !   !here we also record the coordinate for the point
     !     nx_bd_pnt_elastic(i) = nx_pnt(ipnt)
     !     nz_bd_pnt_elastic(i) = nz_pnt(ipnt)
     !     x_final_bd_pnt_elastic(i) = bd_pnt_xval(ipnt)
     !     z_final_bd_pnt_elastic(i) = bd_pnt_zval(ipnt)
     !     i = i+1
     !   endif
     
     !   if(acoustic(ispec_selected_bd_pnt(ipnt)))then
     !     ispec_bd_elmt_acoustic(j) = bd_pnt_elmnt_num(ipnt)
     !     ispec_bd_elmt_acoustic_i(j) = bd_pnt_i(ipnt)
     !     ispec_bd_elmt_acoustic_j(j) = bd_pnt_j(ipnt)
     !     side_type_acoustic(j) = side_type(ipnt)

     !     nx_bd_pnt_acoustic(j) = nx_pnt(ipnt)
     !     nz_bd_pnt_acoustic(j) = nz_pnt(ipnt)
     !     x_final_bd_pnt_acoustic(j) = bd_pnt_xval(ipnt)
     !     z_final_bd_pnt_acoustic(j) = bd_pnt_zval(ipnt)

     !     j = j+1
     !   endif
     ! enddo

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

  if ( record_local_boundary_reconst ) then

     if ( nspec_bd_pnt_elastic /= 0 ) then
        !elastic
        allocate(trac_bd_pnt_elastic_reconst(3,nspec_bd_pnt_elastic))
        allocate(trac_f(3,nspec_bd_pnt_elastic))

        trac_bd_pnt_elastic_reconst = 0.0
        trac_f = 0.0
        

        !P-SV and SH cases will record different boundary information
        if( p_sv ) then
           !temporary array for MPI writing   
           allocate(temp_record_elastic(3,nspec_bd_pnt_elastic))

           allocate(m_xx(nspec_bd_pnt_elastic))
           allocate(m_xz(nspec_bd_pnt_elastic))
           allocate(m_zz(nspec_bd_pnt_elastic))
           allocate(m_zx(nspec_bd_pnt_elastic))
           allocate(m_xx_reconst(nspec_bd_pnt_elastic))
           allocate(m_xz_reconst(nspec_bd_pnt_elastic))
           allocate(m_zz_reconst(nspec_bd_pnt_elastic))
           allocate(m_zx_reconst(nspec_bd_pnt_elastic))
        else

           allocate(temp_record_elastic(2,nspec_bd_pnt_elastic))
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
           allocate(temp_record_acoustic(2,nspec_bd_pnt_acoustic))
        endif
     endif

  endif

  !!export the elastic/acoustic points profiles, which will be used to reconstruct the wavefield
  if( record_local_boundary_reconst ) then
    
#ifdef USE_MPI
     write(fname,"('./OUTPUT_FILES/reconst_record/elastic_pnts_profile',i5.5)") myrank
#else
     fname = './OUTPUT_FILES/reconst_record/elastic_pnts_profile'
#endif
     
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

#ifdef USE_MPI
     write(fname,"('./OUTPUT_FILES/reconst_record/acoustic_pnts_profile',i5.5)") myrank
#else
     fname = './OUTPUT_FILES/reconst_record/acoustic_pnts_profile' 
#endif

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




  !deallocate the intermediate arrays for wavefield record in global model
  !total
  deallocate(side_type_total)
  deallocate(ispec_selected_bd_pnt_total)
  deallocate(elastic_flag_total, acoustic_flag_total)
  deallocate(corner_flag_total)

  deallocate(bd_pnt_elmnt_num_total)
  deallocate(bd_pnt_i_total, bd_pnt_j_total)
  deallocate(bd_pnt_xval_total, bd_pnt_zval_total)

  deallocate(xi_bd_pnt_total, gamma_bd_pnt_total)

  deallocate(nx_pnt_total, nz_pnt_total)
  deallocate(x_final_bd_pnt_total, z_final_bd_pnt_total)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !local
  deallocate(side_type)
  !deallocate(ispec_selected_bd_pnt)
  deallocate(elastic_flag,acoustic_flag)

  deallocate(bd_pnt_elmnt_num)
  deallocate(bd_pnt_i,bd_pnt_j)
  !deallocate(bd_pnt_xval,bd_pnt_zval)

  deallocate(xi_bd_pnt, gamma_bd_pnt)

  deallocate(nx_pnt, nz_pnt)
  deallocate(x_final_bd_pnt, z_final_bd_pnt)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !elastic/acoustic
  deallocate(bd_pnt_i_elastic,bd_pnt_j_elastic)
  deallocate(bd_pnt_i_acoustic,bd_pnt_j_acoustic)
  
  !the coordinate of recording points are not needed since we already know the
  !target element of location and the lagrange coefficients
  deallocate(x_final_bd_pnt_elastic,z_final_bd_pnt_elastic)
  deallocate(x_final_bd_pnt_acoustic,z_final_bd_pnt_acoustic)

  end subroutine locate_recording_point
