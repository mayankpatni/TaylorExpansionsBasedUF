module load_vector
contains

!---------------------------------------------------------
!      FUNCTION loadvector FOR getting the load vector
!---------------------------------------------------------
	subroutine loadvector  
    use var_inputdata
    use var_analysis
    use integ_2D_cs_TE
    use integ_1D_beam
    use mod_gauss
    use mod_math
    
    implicit none
    integer::i9,j9,k9,i11,i10,ele_beam,ele_cs,count,loadtyp,i
    real*8::xp,zp,yp,eta,alpha,beta,wth,wtk,xmm,zmm,tjacdet
    real*8::loadx,loady,loadz
        
    do i11=1,np
        xp=loadinp(1,i11)                                   ! xp - x coordinate in physical csys
        yp=loadinp(2,i11)                                   ! yp - y coordinate in physical csys
        zp=loadinp(3,i11)                                   ! zp - z coordinate in physical csys
        loadx=loadinp(4,i11) 
        loady=loadinp(5,i11)  
        loadz=loadinp(6,i11)                     
        loadtyp=loadinp(7,i11)                              ! loadtype: point or cross-section surface
        
        call phy_to_norBeam(yp,eta,ele_beam)                ! phy_to_norBeam subroutine for mapping a coordinate from physical to normal csys (defined in MOD_math.f90)
                                                            ! yp - input, eta and ele_beam are outputs: corresponding normal coordinate and beam element number 
        call shape1(eta)                                    ! beam shape functions (defined in integ_1D_beam.f90)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Point Load !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (loadtyp.eq.1) then
            call expfun(xp,zp)                              ! cross-section expansion functions (defined in integ_2D_TE.f90)
                
            do i9=1,nne 
                do i10=1,mexp
                    count=(3*mexp*(ele_beam-1)*(nne-1))+(3*mexp*(i9-1))+(3*i10)                     ! position in the load vector
                    fvect(count-2)=fvect(count-2)+loadx*ft(i10)*sf(i9)
                    fvect(count-1)=fvect(count-1)+loady*ft(i10)*sf(i9)
                    fvect(count)=fvect(count)+loadz*ft(i10)*sf(i9)
                enddo
            enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! cross-section Surface Load !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        elseif (loadtyp.eq.2) then
            
            allocate(ft_int(mexp,ane),gpt(agp),wpt(agp))
            ft_int=0.0d0
            call gauss(agp)
            
            do ele_cs=1,ane
                do j9=1,agp
                    do k9=1,agp
    	                alpha=gpt(j9)                                                               ! gauss points, alpha and beta
                        beta=gpt(k9)
                        wth=wpt(j9)                                                                 ! gauss weights, wth and wtk
                        wtk=wpt(k9)
                        call expfun1(ele_cs,alpha,beta,xmm,zmm,tjacdet)                             ! calculates cross-section element's jacobian (tjacdet) (defined in integ_2D_TE.f90) 
                        call expfun(xmm,zmm)                                                        ! cross-section shape functions (defined in integ_2D_TE.f90)
                        do i9=1,mexp
                            ft_int(i9,ele_cs)=ft_int(i9,ele_cs)+(tjacdet*wth*wtk*ft(i9))            ! cross-section integration
                        enddo 
                    enddo
                enddo
            enddo
            
            do ele_cs=1,ane
                do i9=1,nne
                    do i10=1,mexp                        
                        count=(3*mexp*(ele_beam-1)*(nne-1))+(3*mexp*(i9-1))+(3*i10)                 ! position in the load vector
                        fvect(count-2)=fvect(count-2)+loadx*ft_int(i10,ele_cs)*sf(i9)
                        fvect(count-1)=fvect(count-1)+loady*ft_int(i10,ele_cs)*sf(i9)
                        fvect(count)=fvect(count)+loadz*ft_int(i10,ele_cs)*sf(i9)
                    enddo
                enddo
            enddo
            
            deallocate(gpt,wpt,ft_int) 
        endif
    enddo
    
    return
    end subroutine
    
end module