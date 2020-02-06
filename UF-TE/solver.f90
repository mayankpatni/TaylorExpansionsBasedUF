module solver
contains

!----------------------------------------------------------------
!   FUNCTION solve_disp_sparse FOR calculating displacement vector
!----------------------------------------------------------------
 SUBROUTINE solve_disp_sparse 
    use var_inputdata
    use var_analysis
    use MKL_PARDISO
    
      IMPLICIT NONE 
!-----Internal variables------------------------------------------------
      REAL*8,allocatable            ::RESULTS(:)                  !SOLUTION
      INTEGER                       ::I
      
       !SOLVER PARAMTERS
      TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE:: PT(:)   ! INTERNAL SOLVER MEMORY POINTER
      
      INTEGER                 :: MAXFCT               ! MAXIMUM NUMBER OF FACTORS IN MEMORY (1)
      INTEGER                 :: MNUM                 ! NUMBER OF MATRIX ( 1<=MAXFCT) (1) 
      INTEGER                 :: MTYPE                ! MATRIX TYPE(11- REAL UNSYMMETRIC, 2 -REAL SYMM)
      INTEGER                 :: PHASE                ! CONTROLS THE EXECUTION OF SOLVER 
      INTEGER                 :: NRHS                 ! NUMBER OF RHS
      INTEGER                 :: ERROR                ! ERROR FLAG
      INTEGER                 :: MSGLVL               ! INFORMATION FROM SOLVER (0/1)
      INTEGER                 :: IDUM(1)              !
      INTEGER,ALLOCATABLE     :: IPARM(:)             ! ARRAY WITH VARIOUS PARAMETERS FOR MKL PARADISO (SIZE:64)
      REAL*8                  :: DDUM(1)
      
      !SETTING UP PARADISO CONTROL PARAMETERS
      
      ALLOCATE(IPARM(64))
      IPARM =     0
      IPARM(1) =  1 ! no solver default
      IPARM(2) =  2 ! fill-in reordering from METIS
      IPARM(4) =  0 ! no iterative-direct algorithm
      IPARM(5) =  0 ! no user fill-in reducing permutation
      IPARM(6) =  0 ! =0 solution on the first n components of x
      IPARM(8) =  2 ! numbers of iterative refinement steps
      IPARM(10) = 13 ! perturb the pivot elements with 1E-13
      IPARM(11) = 1 ! use nonsymmetric permutation and scaling MPS
      IPARM(13) = 1 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
      IPARM(14) = 0 ! Output: number of perturbed pivots
      IPARM(18) =-1 ! Output: number of nonzeros in the factor LU
      IPARM(19) =-1 ! Output: Mflops for LU factorization
      IPARM(20) = 0 ! Output: Numbers of CG Iterations 
      !IPARM(59) = 2 ! OOC 
      
      ERROR       = 0 
      MSGLVL      = 0 !PRINT STATISTICAL INFORMATION
      MAXFCT      = 1
      NRHS        = 1
      MNUM        = 1
      MTYPE       = 2
      
      !IF (M_TYPE == 'SP') MTYPE= 2
      !IF (M_TYPE == 'SI') MTYPE=-2  
      !IF (M_TYPE == 'NS') MTYPE=11
      
      ALLOCATE(PT(64))
      DO I =1,64
          PT(I)%DUMMY = 0
      ENDDO
      
      allocate(RESULTS(tdof))
      RESULTS = 0.0D0
      uvect=0.0
     
      !SOLVER - PHASE 1 - REORDERING AND SYMBOLIC FACTORIZATION - MEMEORY ALLOCATION IS DONE
      PHASE = 11
      CALL PARDISO (PT,MAXFCT,MNUM,MTYPE,PHASE,tdof,k_val,row_ptr,col_idx,IDUM,NRHS,IPARM,MSGLVL,DDUM,DDUM,ERROR)
      IF (error.ne.0) THEN
        WRITE(*,*) 'PHASE 1 - The following ERROR was detected: ', error
        goto 1000
      ENDIF
      
      !SOLVER - PHASE 2 - FACTORIZATION
      PHASE = 22
      CALL PARDISO (PT,MAXFCT,MNUM,MTYPE,PHASE,tdof,k_val,row_ptr,col_idx,IDUM,NRHS,IPARM,MSGLVL,DDUM,DDUM,ERROR)
      IF (error.ne.0) THEN
        WRITE(*,*) 'PHASE 2 - The following ERROR was detected: ', error
        goto 1000
      ENDIF
      
      !SOLVER - PHASE 3 - SOLUTION
      IPARM(8) = 2 ! MAXIMUM NUMBER IF INTERATIVE REFINEMENT STEPS
      PHASE = 33
      CALL PARDISO (PT,MAXFCT,MNUM,MTYPE,PHASE,tdof,k_val,row_ptr,col_idx,IDUM,NRHS,IPARM,MSGLVL,pvect_res,RESULTS,ERROR)
      del_uvect=RESULTS
      IF (error.ne.0) THEN
        WRITE(*,*) 'PHASE 3 - The following ERROR was detected: ', error
        goto 1000
      ENDIF
      
1000  CONTINUE    
      !SOLVER - PHASE 4 - TERMINATION AND RELEASE
      PHASE = -1
      CALL PARDISO (PT,MAXFCT,MNUM,MTYPE,PHASE,tdof,k_val,row_ptr,col_idx,IDUM,NRHS,IPARM,MSGLVL,DDUM,DDUM,ERROR)
      
      deallocate(RESULTS,k_val,row_ptr,col_idx)
      RETURN
 END SUBROUTINE
 
end module

