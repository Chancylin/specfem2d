!this subroutine is to export the coordinate of GLL points at boudaries
!note that it will need the database generated from the local model
subroutine export_gll_pnt_bd()
  use specfem_par, only: nelemabs,ibool,xiz,xix,gammax,gammaz,jacobian,&
                         coord,numabs,codeabs,codeabs_corner,f_num,typeabs, &
                         elastic,acoustic
  implicit none
  include "constants.h"

  real(kind=CUSTOM_REAL) :: nx,nz,xgamma,zgamma,xxi,zxi,jacobian1D
  integer :: i,j,ispec,ispecabs,iglob
  character(len=1) :: geom_side
  character(len=1) :: elastic_flag, acoustic_flag

  f_num = 111
  open(unit=f_num,file='DATA/boundary_points',status='unknown',action='write')
  !recored info:
  !elment_index_in_local_model/ /boundary_type_defined_by_code/ /
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
          write(f_num,113) ispec,elastic_flag,acoustic_flag,'L',geom_side,coord(1,iglob),coord(2,iglob),nx,nz

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
          write(f_num,113) ispec,elastic_flag,acoustic_flag,'R',geom_side,coord(1,iglob),coord(2,iglob),nx,nz
          
       enddo
     endif    

    !--bottom boundary
    if( codeabs(IEDGE1,ispecabs) ) then
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
             write(f_num,113) ispec,elastic_flag,acoustic_flag,'B',geom_side,coord(1,iglob),coord(2,iglob),nx,nz
         endif

      enddo
    endif     
   
    !--top boundary
    if( codeabs(IEDGE3,ispecabs) ) then
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
               write(f_num,113) ispec,elastic_flag,acoustic_flag,'T',geom_side,coord(1,iglob),coord(2,iglob),nx,nz
         endif

      enddo
    endif
    
   
  enddo

  close(f_num)

  print *, 'We have sucessfully obtained the coordinate of GLL points at boundary,'
  print *, 'which will be further used in the hybrid method.'
  print *, 'The program will exit here'
  stop

  113 format(i3.3,2x,A1,2x,A1,2x,A1,2x,A1,2x,4(es12.4,2x)) !48 column, make sure the writting format is proper
end subroutine export_gll_pnt_bd
