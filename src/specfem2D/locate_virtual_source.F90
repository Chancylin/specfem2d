
!========================================================================
! This subroutine is following the original "locate_source_force.F90"
! to locate the positions of the virtual mutiple sources to reconstruct
! the wavefield.
! Note that the original subroutine 'locate_source_force' is not perfect
! since it could locate the point in the elastic element, which we expect
! to be in the fluid region
!========================================================================

!----
!----

  subroutine locate_virtual_source(elastic_flag,acoustic_flag,ibool,coord,nspec,nglob,xigll,zigll,&
               x_source,z_source,ispec_selected_source,myrank,nproc, &
               xi_source,gamma_source,coorg,knods,ngnod,npgeo,iglob_source)

  use specfem_par, only : elastic,acoustic,&
                          AXISYM,is_on_the_axis,xiglj

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  include "constants.h"

  logical,intent(in) :: elastic_flag,acoustic_flag

  integer nspec,nglob,ngnod,npgeo

  integer knods(ngnod,nspec)
  double precision coorg(NDIM,npgeo)

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool

! array containing coordinates of the points
  double precision coord(NDIM,nglob)

  integer i,j,ispec,iglob,iter_loop,ix_initial_guess,iz_initial_guess,number_of_iterations

  double precision x_source,z_source,dist_squared
  double precision xi,gamma,dx,dz,dxi,dgamma

! Gauss-Lobatto-Legendre points of integration
  double precision xigll(NGLLX)
  double precision zigll(NGLLZ)

  double precision x,z,xix,xiz,gammax,gammaz,jacobian
  double precision distmin_squared,final_distance

! source information
  integer ispec_selected_source,iglob_source
  integer, intent(in)  :: nproc, myrank
  double precision xi_source,gamma_source

! #ifdef USE_MPI
!   integer, dimension(1:nproc)  :: allgather_is_proc_source
!   integer, dimension(1)  :: locate_is_proc_source
!   integer  :: ierror
! #endif

! !geometry bounda. Here we play a trick to locate the virtual sources in the local mesh side
! !instead of the global mesh side, regarding the special case that the local mesh share the 
! !same boundary (or GLL points) with the global mesh. Hopefully we won't use this in future
  double precision :: box_t,box_b,box_l,box_r
  ! logical, intent(out) :: element_locate

  box_t = 5000.0
  box_b = -5000.0
  box_l = -5000.0
  box_r = 5000.0
  
  ! element_locate = .FALSE.

! **************
  if (myrank == 0 .or. nproc == 1) then
    write(IOUT,*)
    write(IOUT,*) '*******************************'
    write(IOUT,*) ' locating virtual multiple sourcs'
    write(IOUT,*) '*******************************'
    write(IOUT,*)
  endif

! set distance to huge initial value
  distmin_squared = HUGEVAL

  do ispec = 1,nspec
     !add this condition to make sure the expected elastic/acoustic points will
     !locate in the elastic/acoustic region
     if( (elastic(ispec) .eqv. elastic_flag) .and. (acoustic(ispec) .eqv. acoustic_flag) ) then

     !geometry bounda
     iglob = ibool(2,2,ispec)

     if((dble(coord(1,iglob)) > box_l) .and. (dble(coord(1,iglob)) < box_r) &
         .and. dble(coord(2,iglob)) > box_b .and. dble(coord(2,iglob)) < box_t )then

! loop only on points inside the element
! exclude edges to ensure this point is not shared with other elements
        do j = 2,NGLLZ-1
           do i = 2,NGLLX-1
   
              iglob = ibool(i,j,ispec)
   
              !  we compare squared distances instead of distances themselves to significantly speed up calculations
              dist_squared = (x_source-dble(coord(1,iglob)))**2 + (z_source-dble(coord(2,iglob)))**2
   
              ! keep this point if it is closer to the source
              if(dist_squared < distmin_squared) then
                 iglob_source = iglob
                 distmin_squared = dist_squared
                 ispec_selected_source = ispec
                 ix_initial_guess = i
                 iz_initial_guess = j
              endif
   
           enddo
        enddo

       endif !geometry bound
     endif
! end of loop on all the spectral elements
  enddo
  

!   ! this step comparing the local minumum distance in order to find the global
!   ! minimum is unnecessary, because we have allocated the recording points (here
!   ! as virtual sources) at step 2 of our workflow
! #ifdef USE_MPI
!   ! global minimum distance computed over all processes
!   call MPI_ALLREDUCE (distmin_squared, dist_glob_squared, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)
! #else
!   dist_glob_squared = distmin_squared
! #endif

! ! check if this process contains the source
!   if ( abs(sqrt(dist_glob_squared) - sqrt(distmin_squared)) < TINYVAL ) is_proc_source = 1

! #ifdef USE_MPI
!   ! determining the number of processes that contain the source
!   ! (useful when the source is located on an interface)
!   call MPI_ALLREDUCE (is_proc_source, nb_proc_source, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
! #else
!   nb_proc_source = is_proc_source
! #endif


! #ifdef USE_MPI
!   ! when several processes contain the source, we elect one of them (minimum rank).
!   if ( nb_proc_source > 1 ) then

!      call MPI_ALLGATHER(is_proc_source, 1, MPI_INTEGER, allgather_is_proc_source(1), &
!                         1, MPI_INTEGER, MPI_COMM_WORLD, ierror)
!      locate_is_proc_source = maxloc(allgather_is_proc_source) - 1

!      if ( myrank /= locate_is_proc_source(1) ) then
!         is_proc_source = 0
!      endif
!      nb_proc_source = 1

!   endif

! #endif


! ****************************************
! find the best (xi,gamma) for each source
! ****************************************

! use initial guess in xi and gamma
  if (AXISYM) then
    if (is_on_the_axis(ispec_selected_source)) then
      xi = xiglj(ix_initial_guess)
    else
      xi = xigll(ix_initial_guess)
    endif
  else
    xi = xigll(ix_initial_guess)
  endif
  gamma = zigll(iz_initial_guess)

! iterate to solve the non linear system
  if(USE_BEST_LOCATION_FOR_SOURCE) then
    number_of_iterations = NUM_ITER
  else
    number_of_iterations = 0 ! this means that the loop below will not be executed, i.e. we will not iterate
  endif

  do iter_loop = 1,NUM_ITER

! recompute jacobian for the new point
    call recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian, &
                  coorg,knods,ispec_selected_source,ngnod,nspec,npgeo, &
                  .true.)

! compute distance to target location
    dx = - (x - x_source)
    dz = - (z - z_source)

! compute increments
    dxi  = xix*dx + xiz*dz
    dgamma = gammax*dx + gammaz*dz

! update values
    xi = xi + dxi
    gamma = gamma + dgamma

! impose that we stay in that element
! (useful if user gives a source outside the mesh for instance)
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
  call recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian, &
                    coorg,knods,ispec_selected_source,ngnod,nspec,npgeo, &
                    .true.)

! store xi,gamma of point found
  xi_source = xi
  gamma_source = gamma

! compute final distance between asked and found
  final_distance = sqrt((x_source-x)**2 + (z_source-z)**2)

  write(IOUT,*)
  write(IOUT,*) 'virtual source:'

  if(final_distance == HUGEVAL) call exit_MPI('error locating force source')

  write(IOUT,*) '            original x: ',sngl(x_source)
  write(IOUT,*) '            original z: ',sngl(z_source)
  write(IOUT,*) 'closest estimate found: ',sngl(final_distance),' m away'
#ifdef USE_MPI
  write(IOUT,*) ' in rank ',myrank
#endif
  write(IOUT,*) ' in element ',ispec_selected_source
  write(IOUT,*) ' at xi,gamma coordinates = ',xi_source,gamma_source
  write(IOUT,*)

  write(IOUT,*)
  write(IOUT,*) 'end of virtual source detection'
  write(IOUT,*)

  end subroutine locate_virtual_source


  subroutine locate_virtual_source_only_locate(elastic_flag,acoustic_flag,ibool,coord,nspec,nglob,&
               x_source,z_source,ispec_selected_source,myrank,nproc,element_locate)

  use specfem_par, only : elastic,acoustic

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  include "constants.h"

  logical,intent(in) :: elastic_flag,acoustic_flag

  integer nspec,nglob

  double precision coord(NDIM,nglob)
  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool


  integer i,j,ispec,iglob,ix_initial_guess,iz_initial_guess
  double precision x_source,z_source
  double precision dist_squared,distmin_squared,dist_glob_squared

! source information
  integer ispec_selected_source,iglob_source
  integer :: is_proc_source,nb_proc_source
  integer, intent(in)  :: nproc, myrank

#ifdef USE_MPI
  integer, dimension(1:nproc)  :: allgather_is_proc_source
  integer, dimension(1)  :: locate_is_proc_source
  integer  :: ierror
#endif
  logical :: element_locate
  
  element_locate = .FALSE.


! set distance to huge initial value
  distmin_squared = HUGEVAL
  is_proc_source = 0

  do ispec = 1,nspec
     !add this condition to make sure the expected elastic/acoustic points will
     !locate in the elastic/acoustic region
     if( (elastic(ispec) .eqv. elastic_flag) .and. (acoustic(ispec) .eqv. acoustic_flag) ) then


     ! if((dble(coord(1,iglob)) > box_l) .and. (dble(coord(1,iglob)) < box_r) &
     !     .and. dble(coord(2,iglob)) > box_b .and. dble(coord(2,iglob)) < box_t )then

! loop only on points inside the element
! exclude edges to ensure this point is not shared with other elements
        do j = 2,NGLLZ-1
           do i = 2,NGLLX-1
   
              iglob = ibool(i,j,ispec)
   
              !  we compare squared distances instead of distances themselves to significantly speed up calculations
              dist_squared = (x_source-dble(coord(1,iglob)))**2 + (z_source-dble(coord(2,iglob)))**2
   
              ! keep this point if it is closer to the source
              if(dist_squared < distmin_squared) then
                 iglob_source = iglob
                 distmin_squared = dist_squared
                 ispec_selected_source = ispec
                 ix_initial_guess = i
                 iz_initial_guess = j
              endif
   
           enddo
        enddo

       ! endif !geometry bound
     endif
! end of loop on all the spectral elements
  enddo
  

  ! this step comparing the local minumum distance in order to find the global
  ! minimum is unnecessary, because we have allocated the recording points (here
  ! as virtual sources) at step 2 of our workflow
#ifdef USE_MPI
  ! global minimum distance computed over all processes
  call MPI_ALLREDUCE (distmin_squared, dist_glob_squared, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)
#else
  dist_glob_squared = distmin_squared
#endif

! check if this process contains the source
  if ( abs(sqrt(dist_glob_squared) - sqrt(distmin_squared)) < TINYVAL ) is_proc_source = 1

#ifdef USE_MPI
  ! determining the number of processes that contain the source
  ! (useful when the source is located on an interface)
  call MPI_ALLREDUCE (is_proc_source, nb_proc_source, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
#else
  nb_proc_source = is_proc_source
#endif

#ifdef USE_MPI
  ! when several processes contain the source, we elect one of them (minimum rank).
  if ( nb_proc_source > 1 ) then

     call MPI_ALLGATHER(is_proc_source, 1, MPI_INTEGER, allgather_is_proc_source(1), &
          1, MPI_INTEGER, MPI_COMM_WORLD, ierror)
     locate_is_proc_source = maxloc(allgather_is_proc_source) - 1

     if ( myrank /= locate_is_proc_source(1) ) then
        is_proc_source = 0
     endif
     nb_proc_source = 1

  endif

#endif

  if( is_proc_source == 1 ) then
     element_locate = .true.
  else
     element_locate = .false.
  endif

     ! print *,'locate', element_locate, ' in rank ', myrank

  endsubroutine locate_virtual_source_only_locate
