!by lcx: this subroutine will be called by prepare_timerun()
  subroutine setup_recording_bd()

  use specfem_par, only: ibool,coord,nspec,nglob,xigll,zigll,&
                         coorg,knods,ngnod,npgeo

  call locate_recording_point(ibool,coord,nspec,nglob,xigll,zigll, &
                              coorg,knods,ngnod,npgeo)
  end subroutine setup_recording_bd
