module k_mat
    contains
    
!----------------------------------------------------------------------------------------
!   FUNCTION kmat for Stiffness (K) matrix assembly (Taylor Expansion) - Linear analysis
!----------------------------------------------------------------------------------------
    subroutine kmat_assembly
    use var_inputdata
    use var_analysis
    use integ_1D_beam
    use integ_2D_cs_TE
    
    implicit none
    integer       ::i,ele
     
    ALLOCATE(MTR_R(tdof))
    DO i=1,tdof
        ALLOCATE(MTR_R(i)%VAL(INT(tdof*3.0/100)))
        ALLOCATE(MTR_R(i)%COL(INT(tdof*3.0/100)))
        MTR_R(i)%VAL=0.0D0
        MTR_R(i)%NNZ=1
        MTR_R(i)%COL=0
        MTR_R(i)%COL(1)=i
    ENDDO
    cpos=0;  cpos(11)=1; cpos(22)=2; cpos(33)=3; cpos(23)=4; cpos(32)=4; cpos(13)=5; cpos(31)=5; cpos(12)=6; cpos(21)=6
    allocate(ninj_test(ne*nne*nne*3*3*3*3))
        do ele=1,ane
            allocate(ftfs_test(mexp*mexp*3*3))
            call FEM_1D_int(ele)
            call cuf_2D_int_TE(ele)
            call k_mat_initial(ele)   
            DEALLOCATE (ftfs_test)
        enddo
    DEALLOCATE (ninj_test)
    call kmat_BC
    call kmat_sparse
    DEALLOCATE(MTR_R)
    
return
end subroutine
    
!-----------------------------------------------------------------------------------------------------
!     FUNCTION k_mat_initial for calculation of K matrix - and initial assembly of sparse matrix (local K)
!-----------------------------------------------------------------------------------------------------
   	subroutine k_mat_initial(ele)
    use var_inputdata
    use var_analysis
    
        implicit none
        integer::i,j,t,s,count1,count2,ii,jj,k,ele,count,cntt
        real*8,allocatable::k_nu(:)
        allocate(k_nu(9))
  
        do k=1,ne
            count=3*mexp*(k-1)*(nne-1)
            do i=1,nne
                do j=1,nne    
                    do t=1,mexp
                      do s=1,mexp
                        count1=(3*mexp*(j-1))+(3*s-2)+count
                        count2=(3*mexp*(i-1))+(3*t-2)+count
                        if (count2.ge.count1) then 
                            call kmat0_ll(i,j,k,t,s,k_nu)
                        
                            cntt=1
                            do ii=count1,count1+2
                                do jj=count2,count2+2
                                    if (jj.ge.ii) then
                                        call ASSEMBLY(ii,jj,k_nu(cntt))
                                    endif
                                    cntt=cntt+1
                                enddo
                            enddo
                          endif 
                        enddo
                    enddo  
                enddo
            enddo 
        enddo
        
      deallocate(k_nu)
     return
    end subroutine 

!-----------------------------------------------------------------------------
!     FUNCTION kmat0_ll for calculation of linear part of tangent K matrix 
!------------------------------------------------------------------------------
    subroutine kmat0_ll(i,j,k,t,s,k_nu0)
    use var_inputdata
    use var_analysis
    
    implicit none
    integer::i,j,k,t,s,aa,bb,hh,gg,ii
    real*8,allocatable,dimension(:)::k_nu0
    
    k_nu0=0.0d0
       
    hh=9*(mexp*(t-1)+(s-1))
    do aa=1,3
        do bb=1,3
            gg=81*(nne*nne*(k-1)+nne*(i-1)+(j-1))+27*(aa-1)+9*(bb-1)
            k_nu0(3*(aa-1)+bb)=DOT_PRODUCT(ftfs_test(hh+1:hh+9),ninj_test(gg+1:gg+9))
        enddo
    enddo
    
    return
    end subroutine 
    
!------------------------------------------------------------------------
!   FUNCTION ASSEMBLY FOR initial assembly of sparse matrix
!------------------------------------------------------------------------
    SUBROUTINE ASSEMBLY(IR,IC,VAL)
      use var_inputdata
      use var_analysis

      IMPLICIT NONE
      INTEGER               ::  I,NZ,MAXCOL
      real*8                ::  VAL
      INTEGER               ::  IC,IR
      INTEGER               ::  IA,IB,IAB,A,B,AB,SO

      NZ=MTR_R(IR)%NNZ
      MAXCOL=MTR_R(IR)%COL(NZ)
      
      IF(IC.GT.MAXCOL) THEN
         MTR_R(IR)%VAL(NZ+1)=MTR_R(IR)%VAL(NZ+1)+VAL
         MTR_R(IR)%COL(NZ+1)=IC
         MTR_R(IR)%NNZ=MTR_R(IR)%NNZ+1
      ELSE
        IA=1
        IB=MTR_R(IR)%NNZ
        IAB=INT((IA+IB)/2)
        A=MTR_R(IR)%COL(IA)
        B=MTR_R(IR)%COL(IB)
        AB=MTR_R(IR)%COL(IAB)

10      IF((IB-IA).LT.10) GOTO 20
        
        IF (IC.GE.A.AND.IC.LT.AB) THEN
            IA=IA
            IB=IAB
            IAB=INT((IA+IB)/2)
            A=MTR_R(IR)%COL(IA)
            B=MTR_R(IR)%COL(IB)
            AB=MTR_R(IR)%COL(IAB)
            GOTO 10
        ELSEIF (IC.GE.AB.AND.IC.LE.B) THEN
            IA=IAB
            IB=IB
            IAB=INT((IA+IB)/2)
            A=MTR_R(IR)%COL(IA)
            B=MTR_R(IR)%COL(IB)
            AB=MTR_R(IR)%COL(IAB)
            GOTO 10
        ENDIF

20      DO I=IA,IB
            IF (IC.EQ.MTR_R(IR)%COL(I))THEN
                MTR_R(IR)%VAL(I)=MTR_R(IR)%VAL(I)+VAL
                EXIT
            ELSEIF (IC.LT.MTR_R(IR)%COL(I) )THEN
                MTR_R(IR)%VAL(I+1:NZ+1)=MTR_R(IR)%VAL(I:NZ)
                MTR_R(IR)%COL(I+1:NZ+1)=MTR_R(IR)%COL(I:NZ)
                MTR_R(IR)%VAL(I)=0.0D0
                MTR_R(IR)%VAL(I)=MTR_R(IR)%VAL(I)+VAL
                MTR_R(IR)%COL(I)=IC
                MTR_R(IR)%NNZ=MTR_R(IR)%NNZ+1
                EXIT
            ENDIF
        ENDDO
      ENDIF

      SO=SIZE(MTR_R(IR)%VAL)
      IF(SO-MTR_R(IR)%NNZ<1) THEN
          ALLOCATE (MTR_R(IR)%TEMP(SO))

          MTR_R(IR)%TEMP=0.0d0
          MTR_R(IR)%TEMP=MTR_R(IR)%VAL
          DEALLOCATE (MTR_R(IR)%VAL)
          ALLOCATE  (MTR_R(IR)%VAL(SO*2))
          MTR_R(IR)%VAL=0.0D0
          MTR_R(IR)%VAL(1:SO)=MTR_R(IR)%TEMP(1:SO)

          MTR_R(IR)%TEMP=0.0d0
          MTR_R(IR)%TEMP=MTR_R(IR)%COL
          DEALLOCATE (MTR_R(IR)%COL)
          ALLOCATE  (MTR_R(IR)%COL(SO*2))
          MTR_R(IR)%COL=0
          MTR_R(IR)%COL(1:SO)=INT(MTR_R(IR)%TEMP(1:SO))
          DEALLOCATE (MTR_R(IR)%TEMP)
      ENDIF
      
     return
    end subroutine
    
!----------------------------------------------------------------
!   FUNCTION kmat_sparse FOR FORMATION OF global sparse MATRIX
!----------------------------------------------------------------
    SUBROUTINE  kmat_sparse
      use var_inputdata
      use var_analysis
      
      IMPLICIT NONE
      INTEGER                :: I,J,JJ
      REAL*8                 :: MFACTOR

      tot_nnz=0

      DO I=1,tdof
        tot_nnz=tot_nnz+MTR_R(I)%NNZ
      ENDDO

      ALLOCATE(k_val(tot_nnz))
      ALLOCATE(col_idx(tot_nnz))
      ALLOCATE(row_ptr(tdof+1))
      
      col_idx=0
      k_val=0.0D0
      row_ptr=0

      JJ=1
      row_ptr(1)=1
      DO I=1,tdof
          row_ptr(I+1)=row_ptr(I)+MTR_R(I)%NNZ
          DO J=1,MTR_R(I)%NNZ
            k_val(JJ)=MTR_R(I)%VAL(J)
            col_idx(JJ)=MTR_R(I)%COL(J)
            JJ=JJ+1
          ENDDO
      ENDDO

      MFACTOR=tot_nnz/(tdof**2.0D0)*100.0D0
    
      return
    end subroutine

    
!-----------------------------------------------------------
!      FUNCTION kmat_BC FOR getting Kmat after applying BCs
!-----------------------------------------------------------
	subroutine kmat_BC  
    use var_inputdata
    use var_analysis
    
    implicit none
    integer::row_nz,i,j,beam_nd,row_nnz,ii,ij,j1
    
    do i=1,num_bcinp
        beam_nd=node_bc(i,1)
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Fixed BC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (node_bc(i,2).eq.0.and.node_bc(i,3).eq.0) then           
           if (beam_nd.eq.1) then
               do j=1,3*mexp
                    row_nz=MTR_R(j)%NNZ
                    MTR_R(j)%VAL(2:row_nz)=0.0d0
                    MTR_R(j)%VAL(1)=1.0d0
              enddo 
           elseif (beam_nd.eq.tnodes) then  
               do j=(tdof-3*mexp)+1,tdof
                  row_nz=MTR_R(j)%NNZ
                  MTR_R(j)%VAL(2:row_nz)=0.0d0
                  MTR_R(j)%VAL(1)=1.0d0
                  
                  do ii=1,tdof
                    row_nnz=MTR_R(ii)%NNZ
                    do ij=2,row_nnz
                        if (MTR_R(ii)%COL(ij).eq.j) then
                            MTR_R(ii)%VAL(ij)=0.0d0
                        endif
                    enddo
                  enddo
               enddo
           else
               do j=3*mexp*(beam_nd-1)+1,3*mexp*beam_nd
                    row_nz=MTR_R(j)%NNZ
                    MTR_R(j)%VAL(2:row_nz)=0.0d0
                    MTR_R(j)%VAL(1)=1.0d0
                  
                    do ii=1,tdof
                        row_nnz=MTR_R(ii)%NNZ
                        do ij=2,row_nnz
                            if (MTR_R(ii)%COL(ij).eq.j) then
                                MTR_R(ii)%VAL(ij)=0.0d0
                            endif
                        enddo
                    enddo 
               enddo
           endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Simply-supported BC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        elseif (node_bc(i,2).eq.100.and.node_bc(i,3).eq.100) then
            if (beam_nd.eq.1) then                                                       ! first beam node
               do j=1,3*mexp
                   if (mod(j-1,3).eq.0.or.mod(j,3).eq.0) then                           !!!!!! x and z constrained !!!!!!
                       row_nz=MTR_R(j)%NNZ
                       MTR_R(j)%VAL(2:row_nz)=0.0d0
                       MTR_R(j)%VAL(1)=1.0d0
                       
                       do ii=1,tdof
                        row_nnz=MTR_R(ii)%NNZ
                        do ij=2,row_nnz
                            if (MTR_R(ii)%COL(ij).eq.j) then
                                MTR_R(ii)%VAL(ij)=0.0d0
                            endif
                        enddo
                      enddo
                   endif
               enddo 
           elseif (beam_nd.eq.tnodes) then                                               ! last beam node
               do j=(tdof-3*mexp)+1,tdof
                    if (mod(j-1,3).eq.0.or.mod(j,3).eq.0) then                          !!!!!! x and z constrained !!!!!!
                      row_nz=MTR_R(j)%NNZ
                      MTR_R(j)%VAL(2:row_nz)=0.0d0
                      MTR_R(j)%VAL(1)=1.0d0
                  
                      do ii=1,tdof
                        row_nnz=MTR_R(ii)%NNZ
                        do ij=2,row_nnz
                            if (MTR_R(ii)%COL(ij).eq.j) then
                                MTR_R(ii)%VAL(ij)=0.0d0
                            endif
                        enddo
                      enddo
                    endif                    
                enddo
           else
               do j=3*(beam_nd-1)*mexp+1,3*beam_nd*mexp
                   if (mod(j+1,3).eq.0) then                                           !!!!!! y constrained !!!!!! 
                        row_nz=MTR_R(j)%NNZ
                        MTR_R(j)%VAL(2:row_nz)=0.0d0
                        MTR_R(j)%VAL(1)=1.0d0
                  
                        do ii=1,tdof
                            row_nnz=MTR_R(ii)%NNZ
                            do ij=2,row_nnz
                                if (MTR_R(ii)%COL(ij).eq.j) then
                                    MTR_R(ii)%VAL(ij)=0.0d0
                                endif
                            enddo
                        enddo 
                   endif
               enddo
           endif
        endif
    enddo
              
    return
    end subroutine
    
end module

