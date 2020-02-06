module mod_math
contains

!---------------------------------------------------------------------------------------
!      FUNCTION phy_to_norCS_1 FOR transforming physical coordinate to normal coordinate
!---------------------------------------------------------------------------------------
	subroutine phy_to_norCS_1
    use var_inputdata
    use var_analysis
    
        implicit none
        integer::normcnt,phycnt,ii1,ii2,ii3,ii4,nod
        real*8::norm1,norm2,dnorm,xphynode,zphynode,alpha1,beta1,ft1(4)
        
        tpts=(phynormesh+1)*(phynormesh+1)
        allocate(normcs(tpts,2),phycs(tpts*ane,3))
        phycs=0.0
        
        norm1=-1.0d0
        norm2=1.0d0
        dnorm=(norm2-norm1)/phynormesh
        normcnt=0
        phycnt=0
          
        do ii1=1,phynormesh+1
            do ii2=1,phynormesh+1
                normcnt=normcnt+1
                normcs(normcnt,1)=norm1+dnorm*(ii2-1)
                normcs(normcnt,2)=norm1+dnorm*(ii1-1)
            enddo
        enddo

        do ii3=1,ane
            do ii1=1,tpts
                alpha1=normcs(ii1,1)
                beta1=normcs(ii1,2)
                ft1(1)=0.25d0*(1.0-alpha1)*(1.0-beta1)   
                ft1(2)=0.25d0*(1.0+alpha1)*(1.0-beta1)
                ft1(3)=0.25d0*(1.0+alpha1)*(1.0+beta1)
                ft1(4)=0.25d0*(1.0-alpha1)*(1.0+beta1) 
                phycnt=phycnt+1
                phycs(phycnt,1)=ii3
                
                do ii2=1,4
                    nod=elemcon(ii3,ii2+1)
                    xphynode=cscord(nod,2)
                    zphynode=cscord(nod,3)
                    phycs(phycnt,2)=phycs(phycnt,2)+xphynode*ft1(ii2)
                    phycs(phycnt,3)=phycs(phycnt,3)+zphynode*ft1(ii2)
                enddo
            enddo
        enddo
        
      return 
    end subroutine


!-----------------------------------------------------------------------------------------------------------
!      FUNCTION phy_to_norCS_2_post FOR transforming physical coordinate to normal coordinate - post processing
!-----------------------------------------------------------------------------------------------------------
	subroutine phy_to_norCS_2_post(xp,zp,alpha_p,beta_p,ele_cs_p,alpbet_cnt)
    use var_inputdata
        implicit none
        integer::phycnt,ii1,minlocat,ii2,alpbet_cnt
        real*8::mini,xp,zp
        real*8,allocatable,dimension(:)::phycs1,alpha_p,beta_p
        integer,allocatable,dimension(:)::ele_cs_p
        allocate(phycs1(tpts*ane))
        
        phycnt=tpts*ane
        phycs1=0.0
        ele_cs_p=0
        alpha_p=0.0
        beta_p=0.0
        
        do ii1=1,phycnt
            phycs1(ii1)=((phycs(ii1,2)-xp)**2)+((phycs(ii1,3)-zp)**2)
        enddo
        
        mini=minval(phycs1)
        
        alpbet_cnt=0
        do ii2=1,phycnt
            if (phycs1(ii2).eq.mini) then
                alpbet_cnt=alpbet_cnt+1
                ele_cs_p(alpbet_cnt)=phycs(ii2,1)
                minlocat=ii2-(ele_cs_p(alpbet_cnt)-1)*tpts
                alpha_p(alpbet_cnt)=normcs(minlocat,1)
                beta_p(alpbet_cnt)=normcs(minlocat,2)
            endif
        enddo
        
        
        deallocate(phycs1)
      return 
    end subroutine
    
    
!-----------------------------------------------------------------------------------------
!      FUNCTION phy_to_norBeam FOR transforming physical coordinate to normal coordinate
!-----------------------------------------------------------------------------------------
	subroutine phy_to_norBeam(yp,eta,ele_beam)
    use var_inputdata
        implicit none
        integer::ii1,ele_beam,ynod1,ynod2
        real*8::eta,yp1,yp2,yp
        logical::FLAG
        ele_beam=0
        
        FLAG=.TRUE.
        do ii1=1,ne
            ynod1=((nne-1)*(ii1-1))+1
            ynod2=((nne-1)*ii1)+1
        
            if (yp.le.y(ynod2).and.FLAG .eq. .TRUE.) then
                ele_beam=ii1
                yp1=y(ynod1)
                yp2=y(ynod2)
                eta=-1.0d0+(2.0d0*(yp-yp1)/(yp2-yp1))
                FLAG=.FALSE.
            endif
        enddo
        
      return 
    end subroutine
 
!--------------------------------------
!       matrix print
!--------------------------------------
		SUBROUTINE PRINTF(A,N,M) 
		IMPLICIT NONE
		INTEGER:: N,I,J, M
		real*8, DIMENSION (N,M):: A
       
	   DO I = 1, N
	  
        WRITE(11,1003)  A(I,:)
	1003 format(//,24(en18.7,1x))
       END DO
	   
       WRITE(11,*)       
        END SUBROUTINE 
        
!--------------------------------------
!       integer matrix print
!--------------------------------------
		SUBROUTINE PRINTI(A,N,M) 
		IMPLICIT NONE
		INTEGER:: N,I,J, M
		integer, DIMENSION (N,M):: A
       
	   DO I = 1, N
	  
        WRITE(11,1001)  A(I,:)
	1001 format(//,24(I4,1x))
       END DO
	   
       WRITE(11,*)       
        END SUBROUTINE 
         

! --------------------------------------------------------------------
! INTEGER FUNCTION  FindMinimum():
!    This function returns the location of the minimum in the section
! between Start and End.
! --------------------------------------------------------------------
   INTEGER FUNCTION  FindMinimum(x, Start, endd)
      IMPLICIT  NONE
      real*8, DIMENSION(1:), INTENT(IN)  :: x
      INTEGER, INTENT(IN)                :: Start, endd
      real*8                             :: Minimum
      INTEGER                            :: Location
      INTEGER                            :: i

      Minimum  = x(Start)		                ! assume the first is the min
      Location = Start			                ! record its position
      DO i = Start+1, endd		                ! start with next elements
         IF (x(i) < Minimum) THEN	            !   if x(i) less than the min?
            Minimum  = x(i)		                !      Yes, a new minimum found
            Location = i                        !      record its position
         END IF
      END DO
      FindMinimum = Location        	        ! return the position
   END FUNCTION  FindMinimum

! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

   SUBROUTINE  Swap(a, b)
      IMPLICIT  NONE
      real*8, INTENT(INOUT)     :: a, b
      real*8                    :: Temp

      Temp = a
      a    = b
      b    = Temp
   END SUBROUTINE  Swap

! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

   SUBROUTINE  Sort(x, Size)
      IMPLICIT  NONE
      real*8, DIMENSION(1:), INTENT(INOUT)  :: x
      INTEGER, INTENT(IN)                   :: Size
      INTEGER                               :: i
      INTEGER                               :: Location

      DO i = 1, Size-1			                            ! except for the last
         Location = FindMinimum(x, i, Size)	                ! find min from this to last
         CALL  Swap(x(i), x(Location))	                    ! swap this and the minimum
      END DO
   END SUBROUTINE  Sort
   
   
    
!-----------------------------------------------------------------------------------------
!      FUNCTION inverse FOR matrix inversion
!-----------------------------------------------------------------------------------------  
  subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
real*8::a(n,n), c(n,n)
real*8::L(n,n), U(n,n), b(n), d(n), x(n)
real*8::coeff
integer::i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
  end subroutine inverse
  
!-----------------------------------------------------------------------------------------
!      FUNCTION matinv3 FOR matrix inversion of 3x3 matrix
!-----------------------------------------------------------------------------------------    
    subroutine matinv3(A,B)
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    real*8 :: A(3,3), B(3,3)  !! Matrix
    real*8 :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
  end subroutine matinv3
        
end module
