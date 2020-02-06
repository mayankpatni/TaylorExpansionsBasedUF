!     Last change:  MP  25 Dec 2019 
!===========================================================================
!      Talyor Expansion based Unified Formulation (for beam-like structures)
!---------------------------------------------------------------------------
!      UF_TE PROGRAMME
!===========================================================================
    PROGRAM UF_TE
    use var_inputdata
    use var_analysis
    use read_input_data
    use integ_1D_beam
    use integ_2D_cs_TE
    use k_mat
    use mod_gauss
    use mod_math
    use load_vector
    use solution
    use solver
    use post_process
    
    implicit none
    integer::i
    
    write(*,90)'Step','System clock time','CPU time'
    CALL SYSTEM_CLOCK(COUNT_RATE=t_rate); CALL SYSTEM_CLOCK(COUNT=ti_sys); call cpu_time(ti)
80  format(A22,11X,F15.2,A1,F15.2,A1)
90  format(A22,11X,A17,A17)
100 format(A15,7X,A15)    
!!!!!!!!!! Open I/O files !!!!!!!
    open(1,file="input/input.dat")
    open(2,file="input/cs_mesh.msh")
    open(4,file="input/load_inp.dat")
    open(5,file="input/post_p.dat")
    
    open(11,file="output/output.txt")
    open(13,file="output/uvect.txt")
    open(14,file="output/post_p_out.dat")

!!!!!!!!! Input processing !!!!!!!!
    CALL SYSTEM_CLOCK(COUNT_RATE=t0_rate); CALL SYSTEM_CLOCK(COUNT=t1_sys); call cpu_time(t1)
    
    call read_input                                 ! call read_input subroutine (defined in read_input.f90)
    call phy_to_norCS_1                             ! call phy_to_norCS_1 subroutine (defined in MOD_math.f90)
    
    CALL SYSTEM_CLOCK(COUNT=t2_sys); call cpu_time(t2)
    write(*,80)'Input processing',real(t2_sys-t1_sys)/t_rate,'s',t2-t1,'s'
    
!!!!! Memory allocation and initialisation !!!!!!
    allocate (sf(nne),dsf(nne),ft(mexp),ftx(mexp),ftz(mexp),trcmat(6,6),fvect(tdof))
    fvect=0.0d0                                     ! force vector
    trcmat=0.0d0                                    ! material stiffness C matrix (6x6)
    
!!!!!!!! Load vector !!!!!!!!
    CALL SYSTEM_CLOCK(COUNT_RATE=t0_rate); CALL SYSTEM_CLOCK(COUNT=t1_sys); call cpu_time(t1)
    
    call loadvector                                 ! call loadvector subroutine (defined in load_vector.f90)
  
    CALL SYSTEM_CLOCK(COUNT=t2_sys); call cpu_time(t2)
    write(*,80)'Load vector assembly',real(t2_sys-t1_sys)/t0_rate,'s',t2-t1,'s'
    
!!!!!! K matrix assembly !!!!!!!!!
    CALL SYSTEM_CLOCK(COUNT_RATE=t0_rate); CALL SYSTEM_CLOCK(COUNT=t1_sys); call cpu_time(t1)
    
    call kmat_assembly                              ! call kmat_assembly subroutine (defined in K_mat.f90)
    
    CALL SYSTEM_CLOCK(COUNT=t2_sys); call cpu_time(t2)
    write(*,80)'K matrix assembly',real(t2_sys-t1_sys)/t_rate,'s',t2-t1,'s'
    
!!!!!! Memory allocation and initialisation !!!!!!!!
    allocate (uvect(tdof),del_uvect(tdof),pvect_res(tdof))
    uvect=0.0d0                                     ! displacement vector       
    del_uvect=0.0d0                                 ! delta u vector (u(n+1)-u(n))
    pvect_res=0.0d0                                 ! residual force vector
    
!!!!!! Solution !!!!!
    CALL SYSTEM_CLOCK(COUNT_RATE=t0_rate); CALL SYSTEM_CLOCK(COUNT=t1_sys); call cpu_time(t1)
     
    call lin_static_TE                              ! call lin_static_TE subroutine (defined in solution.f90)
    
    CALL SYSTEM_CLOCK(COUNT=t2_sys); call cpu_time(t2)
    write(*,80)'Solution',real(t2_sys-t1_sys)/t0_rate,'s',t2-t1,'s'
    
!!!! Output and Memory de-allocation!!!!!!!
    write(13,100)'fvect','uvect'
    do i=1,tdof
        write(13,*)fvect(i),uvect(i)
    enddo
    deallocate (del_uvect,pvect_res,fvect)
   
!!!!!! Post-process !!!!!!
    CALL SYSTEM_CLOCK(COUNT_RATE=t0_rate); CALL SYSTEM_CLOCK(COUNT=t1_sys); call cpu_time(t1)
    
    call post_processing2                           ! call post_processing2 subroutine (defined in post_process.f90)

    CALL SYSTEM_CLOCK(COUNT=t2_sys); call cpu_time(t2)
    write(*,80)'Post-process points',real(t2_sys-t1_sys)/t0_rate,'s',t2-t1,'s'
   
!!!!! Memory de-allocation of arrays/vectors !!!!!!!
    deallocate (sf,dsf,ft,ftx,ftz,trcmat,uvect)

!!!!!! Close I/O files !!!!!
    close(1)
    close(2)
    close(4)
    close(5)
    close(11)
    close(13)
    close(14)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    CALL SYSTEM_CLOCK(COUNT=tf_sys); call cpu_time(tf)
    write(*,80)'Full analysis',real(tf_sys-ti_sys)/t_rate,'s',tf-ti,'s'
    
	stop
    END PROGRAM UF_TE