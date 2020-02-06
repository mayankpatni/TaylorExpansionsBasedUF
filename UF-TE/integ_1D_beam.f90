module integ_1D_beam
    
contains

!-------------------------------------------------------------------
!    FUNCTION FEM_1D_int for evaluating 1D integration along beam
!-------------------------------------------------------------------
   subroutine FEM_1D_int(ele)
    use var_inputdata
    use var_analysis
    use mod_gauss
    use read_input_data
    
    implicit none
    integer::j,kk,i1,i2,ele
    real*8::wt,eta,eta_phy

    allocate(gpt(ngp),wpt(ngp))   
    ninj_test=0.0d0
    
    call gauss(ngp)
    
    do kk=1,ne
        do j=1,ngp
    	    eta=gpt(j)
     	    wt=wpt(j)
            eta_phy=jacob(kk)*(eta+1.0d0)+y((nne-1)*(kk-1)+1)
            
            call mat_stiff(ele,eta_phy)
 	        call shape1(eta)
  	        call lineint(wt,kk)
        enddo
    enddo
    
    deallocate(gpt,wpt)
    return
   end subroutine
   
!----------------------------------------------------------------
!   FUNCTION shape FOR SHAPE FUNCTIONS and their derivatives
!----------------------------------------------------------------
	subroutine shape1(eta)
    use var_inputdata
    use var_analysis
    
    implicit none
    real*8::eta
    
    if (nne.eq.2) then
        sf(1)=0.5d0*(1.0-eta)
        sf(2)=0.5d0*(1.0+eta)
        dsf(1)=(-0.5d0)
        dsf(2)=0.5d0
 
    elseif (nne.eq.3) then
        sf(1)=(0.5d0*eta*(eta-1.0))
        sf(2)=(1.0d0+eta)*(1.0d0-eta)
        sf(3)=(0.5d0*eta*(eta+1.0))
        dsf(1)=(eta-0.5d0)
        dsf(2)=-2.0d0*eta
        dsf(3)=(eta+0.5d0)
        
    elseif (nne.eq.4) then
        sf(1)=(-9.0d0/16.0d0)*(eta-1.0d0)*(eta+(1.0d0/3.0d0))*(eta-(1.0d0/3.0d0))
        sf(2)=(27.0d0/16.0d0)*(eta+1.0d0)*(eta-(1.0d0/3.0d0))*(eta-1.0d0)
        sf(3)=(-27.0d0/16.0d0)*(eta+1.0d0)*(eta+(1.0d0/3.0d0))*(eta-1.0d0)
        sf(4)=(9.0d0/16.0d0)*(eta+1.0d0)*(eta+(1.0d0/3.0d0))*(eta-(1.0d0/3.0d0))
        
        dsf(1)=(-9.0d0/16.0d0)*(3.0d0*(eta**2)-(2.0d0*eta)-(1.0d0/9.0d0))
        dsf(2)=(27.0d0/16.0d0)*(3.0d0*(eta**2)-((2.0d0/3.0d0)*eta)-1.0d0)
        dsf(3)=(-27.0d0/16.0d0)*(3.0d0*(eta**2)+((2.0d0/3.0d0)*eta)-1.0d0)
        dsf(4)=(9.0d0/16.0d0)*(3.0d0*(eta**2)+(2.0d0*eta)-(1.0d0/9.0d0))
        
    elseif (nne.eq.5) then
        sf(1)=(2.0d0/3.0d0)*(eta-1)*(eta+0.5d0)*(eta-0.5d0)*eta
        sf(2)=(-8.0d0/3.0d0)*(eta+1)*(eta-1)*(eta-0.5d0)*eta
        sf(3)=4.0*(eta+1)*(eta-1)*(eta+0.5d0)*(eta-0.5d0)
        sf(4)=(-8.0d0/3.0d0)*(eta+1)*(eta-1)*(eta+0.5d0)*eta
        sf(5)=(2.0d0/3.0d0)*(eta+1)*(eta+0.5d0)*(eta-0.5d0)*eta
        
        dsf(1)=(8.0d0/3.0d0)*((eta**3)-0.75d0*(eta**2)-0.125d0*eta+0.0625d0)
        dsf(2)=(-32.0d0/3.0d0)*((eta**3)-0.375d0*(eta**2)-0.5d0*eta+0.125d0)
        dsf(3)=16.0d0*((eta**3)-0.625d0*eta)
        dsf(4)=(-32.0d0/3.0d0)*((eta**3)+0.375d0*(eta**2)-0.5*eta-0.125d0)
        dsf(5)=(8.0d0/3.0d0)*((eta**3)+0.75d0*(eta**2)-0.125*eta-0.0625d0)
    endif
        
	return
    end subroutine
    
    
!-------------------------------------------------------------
!     FUNCTION lineint for the line integrals
!-------------------------------------------------------------
   	subroutine lineint(wt,kk)
    use var_inputdata
    use var_analysis
    
    implicit none
    integer::ii,jj,kk,iii,jjj,aa,bb,cc,dd,ee1,ff1,gg
    real*8::wt
    
    allocate(sf_test(nne,3))
    
    sf_test(:,1)=sf(:)
    sf_test(:,2)=dsf(:)*jacobinv(kk)
    sf_test(:,3)=sf(:) 
    
    do aa=1,3
        do bb=1,3
            do cc=1,3
                do dd=1,3
                    ee1=min(cpos(10*aa+dd),cpos(10*bb+cc)); ff1=max(cpos(10*aa+dd),cpos(10*bb+cc))
                    do ii=1,nne
                        do jj=1,nne
                            gg=81*(nne*nne*(kk-1)+nne*(ii-1)+(jj-1))+27*(aa-1)+9*(bb-1)+3*(cc-1)+dd 
                            ninj_test(gg)=ninj_test(gg)+(trcmat(ee1,ff1)*wt*jacob(kk)*sf_test(ii,cc)*sf_test(jj,dd))
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo

    deallocate(sf_test)
   
return
end subroutine
    
    
end module


