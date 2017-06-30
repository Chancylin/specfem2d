!by lcx:the subroutine prepares parameters to control recording or
!reading back the local/global bd info
  subroutine read_para_hybrid()

    use specfem_par, only: myrank,record_local_bkgd_boundary,export_gll_pnt_local,&
                           supply_local_bkgd_boundary,virtual_ab_bd,&
                           num_step_output,num_step_input,&
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
    read(111,*) num_step_output 
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


    ! here the code check whether the velues of parameters are reasonble
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if( record_nt1 < 1 .or. read_nt1 < 1 .or. record_nt1_reconst < 1 &
         .or. read_nt1_reconst < 1 )then
       stop 'any starting step of reading/recoding must not be smaller than 1'
    endif
    
    if( read_nt1*deltat_read < record_nt1*deltat_record &
         .or. read_nt2*deltat_read > record_nt2*deltat_record )then
       stop 'time window to supply the excitation in local model &
            & should be within recording time window of global model'
    endif

    if( deltat_read > deltat_record )then
       stop 'time step in local simulation could not be larger than global simulation'
    endif
    
    if( record_nt1 /= read_nt1_reconst .or. record_nt2 /= read_nt2_reconst )then
       print *,'we strongly suggest two global simulation (step 2 and 4) should &
            & have same starting and ending step'
    endif
    
     if( read_nt1 /= record_nt1_reconst .or. read_nt2 /= record_nt2_reconst )then
       print *,'we strongly suggest the supplying and recording processes in local &
            & simulation should have same starting and ending step'
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !!add some info to check the flag
    if( myrank == 0 )then
       print *,"warning: make sure what your purpose is at this step"
       if( export_gll_pnt_local ) print *,"the code will export the boundary GLL of local model"
       if( record_local_bkgd_boundary ) print *,"the code will record local/boundary info"
       if( supply_local_bkgd_boundary ) print *,"the code will read local/boundary info"
       if( record_local_boundary_reconst ) print *,"the code will record the excitation for reconstruction"
       if( supply_reconst ) print *,"the code will supply excitations for reconstruction"
    endif

    !!so far, we just let the input interval equal with the output interval
    num_step_input = num_step_output


  end subroutine read_para_hybrid
