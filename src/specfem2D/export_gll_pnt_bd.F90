!this subroutine is to export the coordinate of GLL points at boudaries
!note that it will need the database generated from the local model
subroutine export_gll_pnt_bd()

#ifdef USE_MPI
  use mpi
#endif
  
  use specfem_par, only: nelemabs,ibool,xiz,xix,gammax,gammaz,jacobian,&
                         coord,numabs,codeabs,codeabs_corner,f_num,typeabs, &
                         elastic,acoustic,myrank,ier
  implicit none
  include "constants.h"

  real(kind=CUSTOM_REAL) :: nx,nz,nx_alt,nz_alt,xgamma,zgamma,xxi,zxi,jacobian1D
  integer :: i,j,ispec,ispecabs,iglob
  character(len=256) :: bd_name
  character(len=1) :: geom_side
  character(len=1) :: elastic_flag, acoustic_flag, corner_flag

#ifdef USE_MPI
  write(bd_name,"('./DATA/boundary_points',i5.5)") myrank
#else
 bd_name = './DATA/boundary_points' 
#endif
  
  f_num = 111
  open(unit=f_num,file=trim(bd_name),status='unknown',action='write')
  !recored info:
  !elment_index_in_local_model/ /boundary_type_defined_by_code/ /
  !boundary_type_in_real_model/ /x/ /z/ /nx/ /nz
  !Otc, 17 2016. Update: add gll_point_index_i gll_point_index_j, which will help when
  !recording information for wavefield reconstruction
  !elment_index_in_local_model/ /gll_point_index_i/ /gll_point_index_j/ /boundary_type_defined_by_code/ /
  !boundary_type_in_real_model/ /x/ /z/ /nx/ /nz
  do ispecabs = 1,nelemabs

     ispec = numabs(ispecabs)
    
     if(elastic(ispec))then
         elastic_flag = 'T'
       else 
         elastic_flag = 'F'
     endif

     if(acoustic(ispec))then
          acoustic_flag = 'T'
       else
          acoustic_flag = 'F'
     endif
     
!the tricky problem we need to deal with is that if the external mesh has the 
!rotated elements, which means the left/right/top/bottom may not match
!the true geometry. And in this case, how can we make sure the corresponding
!normal vector is correct, thus the implementation of absorbing boundary condition
     !--left boundary
     if( codeabs(IEDGE4,ispecabs) ) then
       i = 1
       do j = 1,NGLLZ
          iglob = ibool(i,j,ispec)

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)
          nx = - zgamma / jacobian1D
          nz = + xgamma / jacobian1D
         
         ! write(f_num,113) coord(1,iglob),coord(2,iglob),nx,nz
          if( typeabs(ispecabs) == 1 ) then
                 geom_side='B'
            else if ( typeabs(ispecabs) == 2 ) then
                 geom_side='R'
            else if ( typeabs(ispecabs) == 3 ) then
                 geom_side='T'
            else
                 geom_side='L' 
          endif
 
         !judge whether this GLL point is corner
         !if corner, then use average
          !if( (codeabs_corner(1,ispecabs) .and. j == NGLLZ) ) then
          !if(ispec==10) then
          !   print *,ispecabs
          !   print *,'left: codeabs_corner(3,ispecabs) = ',codeabs_corner(3,ispecabs)
          !endif

          if( codeabs_corner(1,ispecabs) .and. j == 1 ) then
            !Left-Bottom
                 corner_flag = 'T'
                 
                 xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                 zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
                 jacobian1D = sqrt(xxi**2 + zxi**2)
                 nx_alt = + zxi / jacobian1D
                 nz_alt = - xxi / jacobian1D

                 nx = 0.5*nx + 0.5*nx_alt
                 nz = 0.5*nz + 0.5*nz_alt

            else if ( codeabs_corner(3,ispecabs) .and. j == NGLLZ ) then
              !Left-Top
                 corner_flag = 'T'

                 xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                 zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
                 jacobian1D = sqrt(xxi**2 + zxi**2)
                 nx_alt = - zxi / jacobian1D
                 nz_alt = + xxi / jacobian1D
                 
                 nx = 0.5*nx + 0.5*nx_alt
                 nz = 0.5*nz + 0.5*nz_alt

            else
                 corner_flag = 'F'
          endif

          write(f_num,113) ispec,i,j,elastic_flag,acoustic_flag,corner_flag,&
                           'L',geom_side,coord(1,iglob),coord(2,iglob),nx,nz

       enddo
     endif
     
     !--right boundary
     if( codeabs(IEDGE2,ispecabs) ) then
       i = NGLLX
       do j = 1,NGLLZ
          iglob = ibool(i,j,ispec)
  
          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)
          nx = + zgamma / jacobian1D
          nz = - xgamma / jacobian1D

         !write(f_num,113) coord(1,iglob),coord(2,iglob),nx,nz
          if( typeabs(ispecabs) == 1 ) then
                 geom_side='B'
            else if ( typeabs(ispecabs) == 2 ) then
                 geom_side='R'
            else if ( typeabs(ispecabs) == 3 ) then
                 geom_side='T'
            else
                 geom_side='L' 
          endif
         !judge whether this GLL point is corner
         !if corner, then use average
          if( (codeabs_corner(2,ispecabs) .and. j == 1) ) then
            !Right-Bottom
                 corner_flag = 'T'
                 
                 xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                 zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
                 jacobian1D = sqrt(xxi**2 + zxi**2)
                 nx_alt = + zxi / jacobian1D
                 nz_alt = - xxi / jacobian1D

                 nx = 0.5*nx + 0.5*nx_alt
                 nz = 0.5*nz + 0.5*nz_alt

            else if ( (codeabs_corner(4,ispecabs) .and. j == NGLLZ) ) then
              !Right-Top
                 corner_flag = 'T'

                 xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                 zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
                 jacobian1D = sqrt(xxi**2 + zxi**2)
                 nx_alt = - zxi / jacobian1D
                 nz_alt = + xxi / jacobian1D
                 
                 nx = 0.5*nx + 0.5*nx_alt
                 nz = 0.5*nz + 0.5*nz_alt

            else
                 corner_flag = 'F'
          endif
          write(f_num,113) ispec,i,j,elastic_flag,acoustic_flag,corner_flag,&
                           'R',geom_side,coord(1,iglob),coord(2,iglob),nx,nz
          
       enddo
     endif    

    !--bottom boundary
    if( codeabs(IEDGE1,ispecabs) ) then

      corner_flag = 'F'

      j = 1
      do i = 1,NGLLX
         iglob = ibool(i,j,ispec)

         xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
         zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
         jacobian1D = sqrt(xxi**2 + zxi**2)
         nx = + zxi / jacobian1D
         nz = - xxi / jacobian1D
      
    !exclude corners to make sure there is no contradiction on the normal
         if( (codeabs_corner(1,ispecabs) .and. i == 1) .or. (codeabs_corner(2,ispecabs) .and. i == NGLLX) ) then 
            ! print *,'detect a corner GLL'
      !skip this corner point because it has aleady been recorded
           else
             !write(f_num,113) coord(1,iglob),coord(2,iglob),nx,nz
             if( typeabs(ispecabs) == 1 ) then
                    geom_side='B'
               else if ( typeabs(ispecabs) == 2 ) then
                    geom_side='R'
               else if ( typeabs(ispecabs) == 3 ) then
                    geom_side='T'
               else
                    geom_side='L' 
             endif
             write(f_num,113) ispec,i,j,elastic_flag,acoustic_flag,corner_flag,&
                              'B',geom_side,coord(1,iglob),coord(2,iglob),nx,nz
         endif

      enddo
    endif     
   
    !--top boundary
    if( codeabs(IEDGE3,ispecabs) ) then

      corner_flag = 'F'

      j = NGLLZ
      do i = 1,NGLLX
         iglob = ibool(i,j,ispec)

         xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
         zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
         jacobian1D = sqrt(xxi**2 + zxi**2)
         nx = - zxi / jacobian1D
         nz = + xxi / jacobian1D

    !exclude corners to make sure thers is no contradiction on the normal
         if( (codeabs_corner(3,ispecabs) .and. i == 1) .or. (codeabs_corner(4,ispecabs) .and. i == NGLLX) ) then
             !print *,'detect a corner GLL'
         !skip recording
           else
               !write(f_num,113) coord(1,iglob),coord(2,iglob),nx,nz
               if( typeabs(ispecabs) == 1 ) then
                      geom_side='B'
                 else if ( typeabs(ispecabs) == 2 ) then
                      geom_side='R'
                 else if ( typeabs(ispecabs) == 3 ) then
                      geom_side='T'
                 else
                      geom_side='L' 
               endif
               write(f_num,113) ispec,i,j,elastic_flag,acoustic_flag,corner_flag,&
                                'T',geom_side,coord(1,iglob),coord(2,iglob),nx,nz
         endif

      enddo
    endif
    
   
  enddo

  close(f_num)
  
#ifdef USE_MPI
  print *, myrank, ' have sucessfully obtained the coordinate of GLL points at boundary'
  !MPI finish here
  call MPI_Barrier(MPI_COMM_WORLD,ier)
  if( myrank == 0 ) print *, 'The program will exit here'
  call MPI_FINALIZE(ier)
#else
  print *, 'we have sucessfully obtained the coordinate of GLL points at boundary'
  print *, 'program will stop here'
#endif

  stop


  113 format(i5.5,2x,2(i1.1,2x),A1,2x,A1,2x,A1,2x,A1,2x,A1,2x,4(es12.4,2x)) !make sure the writting format is proper
end subroutine export_gll_pnt_bd
