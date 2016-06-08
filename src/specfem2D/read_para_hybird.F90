!by lcx:the subroutine prepares parameters to control recording or
!reading back the local/global bd info
  subroutine read_para_hybird()
    use specfem_par, only: record_local_bkgd_boundary,export_gll_pnt_local
    character(len=256)datlin

    open(111,file='./DATA/switch_solver',action='read',status='old') 
    read(111,*) datlin
    read(111,*) export_gll_pnt_local
    read(111,*) datlin
    read(111,*) record_local_bkgd_boundary
    close(111)
  end subroutine read_para_hybird
