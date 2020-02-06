module integ_2D_cs_TE
contains

!--------------------------------------------------------------------------------------------------
!     FUNCTION expansion for expansion function and their derivatives and surface integral for TE
!--------------------------------------------------------------------------------------------------
   	subroutine cuf_2D_int_TE(ele)
    use var_inputdata
    use var_analysis
    use mod_gauss
    
    implicit none
    integer::i5,j5,i3,j3,j4,ele,cs,cc,dd,hh
    real*8::xmm,zmm,alpha,beta,wth,wtk,tjacdet
        
    ALLOCATE(ft_test(mexp,3))
    ftfs_test=0.0d0
        
    allocate(gpt(agp),wpt(agp))   
    call gauss(agp)

	do i3=1,agp
        do j3=1,agp
    	    alpha=gpt(i3)
            beta=gpt(j3)
            wth=wpt(i3)
            wtk=wpt(j3)
                    
            call expfun1(ele,alpha,beta,xmm,zmm,tjacdet)
            call expfun(xmm,zmm)
            ft_test(:,1)=ftx(:)
            ft_test(:,2)=ft(:)
            ft_test(:,3)=ftz(:)
                            
            do i5=1,mexp
                 do j5=1,mexp
                     do cc=1,3
                         do dd=1,3 
                              hh=9*(mexp*(i5-1)+(j5-1))+3*(cc-1)+dd                                   
                              ftfs_test(hh)=ftfs_test(hh)+(tjacdet*wth*wtk*ft_test(i5,cc)*ft_test(j5,dd))
                         enddo
                      enddo  
                 enddo
            enddo
        enddo
    enddo
    
     
    DEALLOCATE (ft_test)
    deallocate(gpt,wpt)
    return
    end subroutine

!-----------------------------------------------------------------
!     FUNCTION expfun1 for cross section integration for TE 
!-----------------------------------------------------------------
   	subroutine expfun1(ele,alpha,beta,xmm,zmm,tjacdet)
    use var_inputdata
        implicit none
        integer::j4,nodenum,ele
        real*8::alpha,beta,xmm,zmm,xmm1,zmm1,xmm2,zmm2,tjacdet
        real*8,ALLOCATABLE, DIMENSION (:)::x1,z1,tfta,tftb,tft
        ALLOCATE(x1(4),z1(4),tfta(4),tftb(4),tft(4))
         
         do j4=1,4
            nodenum=elemcon(ele,j4+1)
            x1(j4)=cscord(nodenum,2)
            z1(j4)=cscord(nodenum,3)
         enddo
         
        tft(1)=0.25d0*(1.0d0-alpha)*(1.0d0-beta)
        tft(2)=0.25d0*(1.0d0+alpha)*(1.0d0-beta)
        tft(3)=0.25d0*(1.0d0+alpha)*(1.0d0+beta)
        tft(4)=0.25d0*(1.0d0-alpha)*(1.0d0+beta)
          
        tfta(1)=0.25d0*(-1.0d0)*(1.0d0-beta)
        tfta(2)=0.25d0*(1.0d0)*(1.0d0-beta)
        tfta(3)=0.25d0*(1.0d0)*(1.0d0+beta)
        tfta(4)=0.25d0*(-1.0d0)*(1.0d0+beta)
          
        tftb(1)=0.25d0*(1.0d0-alpha)*(-1.0d0)
        tftb(2)=0.25d0*(1.0d0+alpha)*(-1.0d0)
        tftb(3)=0.25d0*(1.0d0+alpha)*(1.0d0)
        tftb(4)=0.25d0*(1.0d0-alpha)*(1.0d0)
        
        xmm=tft(1)*x1(1)+tft(2)*x1(2)+tft(3)*x1(3)+tft(4)*x1(4)
        zmm=tft(1)*z1(1)+tft(2)*z1(2)+tft(3)*z1(3)+tft(4)*z1(4)
        xmm1=tfta(1)*x1(1)+tfta(2)*x1(2)+tfta(3)*x1(3)+tfta(4)*x1(4)
        xmm2=tftb(1)*x1(1)+tftb(2)*x1(2)+tftb(3)*x1(3)+tftb(4)*x1(4)
        zmm1=tfta(1)*z1(1)+tfta(2)*z1(2)+tfta(3)*z1(3)+tfta(4)*z1(4)
        zmm2=tftb(1)*z1(1)+tftb(2)*z1(2)+tftb(3)*z1(3)+tftb(4)*z1(4)
        tjacdet=xmm1*zmm2-xmm2*zmm1
        
      DEALLOCATE (tfta,tftb,tft,x1,z1)          
     return
    end subroutine 
    
!------------------------------------------------------------------------------------
!     FUNCTION expfun for evaluating expansion function at given coordinates for TE
!------------------------------------------------------------------------------------
   	subroutine expfun(xmm,zmm)
    use var_inputdata
    use var_analysis
    
        implicit none
        integer::firstterm,lastterm,counter,i1,i2
        real*8::xmm,zmm
        
         do i1=0,nexp
                firstterm=((i1**2)+i1+2)/2
                lastterm=((i1+1)*(i1+2))/2
                counter=i1
                do i2=firstterm,lastterm
                    ft(i2)=(xmm**counter)*(zmm**(i1-counter))
                    ftx(i2)=(counter*(xmm**(abs(counter-1))))*(zmm**(i1-counter))
                    ftz(i2)=(xmm**counter)*((i1-counter)*(zmm**(abs(i1-counter-1))))
                    counter=counter-1
                enddo
            enddo
     return
    end subroutine
    
end module
