!by lcx:the subroutine prepares parameters to control recording or
!reading back the local/global bd info
  subroutine read_para_hybrid()

    use specfem_par, only: record_local_bkgd_boundary,export_gll_pnt_local,&
                           supply_local_bkgd_boundary,virtual_ab_bd,&
                           record_nt1,record_nt2,deltat_record,&
                           read_nt1,read_nt2,deltat_read,&
                     !!para for reconstructing
                           record_local_boundary_reconst,supply_reconst,&
                           record_nt1_reconst,record_nt2_reconst,deltat_record_reconst,&
                           read_nt1_reconst,read_nt2_reconst,deltat_read_reconst

    character(len=256)datlin

    open(111,file='./DATA/switch_solver',action='read',status='old') 
    read(111,*) datlin
    read(111,*) export_gll_pnt_local
    read(111,*) datlin
    read(111,*) record_local_bkgd_boundary
    read(111,*) datlin
    read(111,*) deltat_record,record_nt1,record_nt2
    read(111,*) datlin
    read(111,*) supply_local_bkgd_boundary
    read(111,*) datlin
    read(111,*) deltat_read,read_nt1,read_nt2
    read(111,*) datlin
    read(111,*) record_local_boundary_reconst
    read(111,*) datlin
    read(111,*) deltat_record_reconst,record_nt1_reconst,record_nt2_reconst
    read(111,*) datlin
    read(111,*) supply_reconst
    read(111,*) datlin
    read(111,*) deltat_read_reconst,read_nt1_reconst,read_nt2_reconst
    

    close(111)
    if(supply_local_bkgd_boundary) then
       virtual_ab_bd = .true.
      else
       virtual_ab_bd = .false.
    endif

   !!add some info to check the flag
   print *,"warning: make sure what your purpose is at this step"
   if( export_gll_pnt_local ) print *,"the code will export the boundary GLL of local model"
   if( record_local_bkgd_boundary ) print *,"the code will record local/boundary info"
   if( supply_local_bkgd_boundary ) print *,"the code will read local/boundary info"
   if( record_local_boundary_reconst ) print *,"the code will record the excitation for reconstruction"
   if( supply_reconst ) print *,"the code will supply excitations for reconstruction"


  end subroutine read_para_hybrid
