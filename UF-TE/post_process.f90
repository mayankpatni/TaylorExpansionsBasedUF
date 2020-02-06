module post_process
contains

!------------------------------------------------------------------------------
!      FUNCTION post_processing2 FOR post processing at desired points
!------------------------------------------------------------------------------
	subroutine post_processing2
    use var_inputdata
    use var_analysis
    use mod_math
    
    implicit none
    integer::k3,k4,ele_cs,alpbet_cnt,k2,i,npost_p
    character*1::mark
    real*8::xm1,zm1,alpha,beta
    real*8,ALLOCATABLE, DIMENSION (:)::nodeval,alpha_p,beta_p
    real*8,allocatable,dimension(:,:)::post_p,u_G,eps_G,sig_G,sig_Le,sig_Ls,xyz_post
    integer, ALLOCATABLE, DIMENSION (:)::ele_cs_p
    allocate (alpha_p(4),beta_p(4),ele_cs_p(4))
    
    read(5,*)mark
    read(5,*)npost_p
    allocate (post_p(npost_p,2))
    read(5,*)mark
    do i=1,npost_p
        read(5,*)post_p(i,1),post_p(i,2)
    enddo   
    
    allocate (strn_p_G(tnodes,9,npost_p),strs_p_Le(tnodes,9,npost_p),strs_p_Ls(tnodes,9,npost_p),strs_p_G(tnodes,9,npost_p))
    allocate (disp_p(tnodes,6,npost_p))
    
    disp_p=0.0d0
    strn_p_G=0.0d0
    strs_p_Le=0.0d0
    strs_p_Ls=0.0d0
    strs_p_G=0.0d0
    
    allocate (eps_G(tnodes,6),u_G(tnodes,3),xyz_post(tnodes,3))
    allocate (sig_G(tnodes,6),sig_Le(tnodes,6),sig_Ls(tnodes,6),nodeval(tnodes))
        
        
    do k3=1,npost_p
        xm1=post_p(k3,1)
        zm1=post_p(k3,2)
        
        call phy_to_norCS_2_post(xm1,zm1,alpha_p,beta_p,ele_cs_p,alpbet_cnt)   
        
        do k4=1,alpbet_cnt
            alpha=alpha_p(k4)
            beta=beta_p(k4)
            ele_cs=ele_cs_p(k4)
            
            call post_disp_strs_strn(xm1,zm1,alpha,beta,ele_cs,xyz_post,u_G,eps_G,sig_G,sig_Le,sig_Ls,nodeval)
         
            do k2=1,tnodes
               disp_p(k2,1,k3)=xyz_post(k2,1)
               disp_p(k2,2,k3)=xyz_post(k2,2) 
               disp_p(k2,3,k3)=xyz_post(k2,3)
                   
               do i=1,3
                   disp_p(k2,i+3,k3)=(disp_p(k2,i+3,k3)*(k4-1)+u_G(k2,i))/(k4*1.0d0)
               enddo
                   
               strn_p_G(k2,1,k3)=xyz_post(k2,1)
               strn_p_G(k2,2,k3)=xyz_post(k2,2)
               strn_p_G(k2,3,k3)=xyz_post(k2,3)
                   
               strs_p_G(k2,1,k3)=xyz_post(k2,1)
               strs_p_G(k2,2,k3)=xyz_post(k2,2)
               strs_p_G(k2,3,k3)=xyz_post(k2,3)
                   
               strs_p_Le(k2,1,k3)=xyz_post(k2,1)
               strs_p_Le(k2,2,k3)=xyz_post(k2,2)
               strs_p_Le(k2,3,k3)=xyz_post(k2,3)
                   
               strs_p_Ls(k2,1,k3)=xyz_post(k2,1)
               strs_p_Ls(k2,2,k3)=xyz_post(k2,2) 
               strs_p_Ls(k2,3,k3)=xyz_post(k2,3)
                   
               do i=1,6
                   strn_p_G(k2,i+3,k3)=(strn_p_G(k2,i+3,k3)*(k4-1)+eps_G(k2,i))/(k4*1.0d0)                      
                   strs_p_G(k2,i+3,k3)=(strs_p_G(k2,i+3,k3)*(k4-1)+sig_G(k2,i))/(k4*1.0d0)
                   strs_p_Le(k2,i+3,k3)=(strs_p_Le(k2,i+3,k3)*(k4-1)+sig_Le(k2,i))/(k4*1.0d0)
                   strs_p_Ls(k2,i+3,k3)=(strs_p_Ls(k2,i+3,k3)*(k4-1)+sig_Ls(k2,i))/(k4*1.0d0)
               enddo
            enddo 
        enddo
    enddo
    
    deallocate (u_G,eps_G,sig_G,sig_Le,sig_Ls,nodeval,xyz_post)
    
    write(14,*)'uz at along y at bottom, mid and top'
    do k3=1,tnodes
        write(14,*)disp_p(k3,6,1),disp_p(k3,6,18),disp_p(k3,6,35)
    enddo
    
    write(14,*)'through-thickness sigyy and tauyz (at 50% of the beam length from left end)'
    do k3=1,npost_p
        write(14,*)strs_p_G(1+((tnodes-1)*50/100),5,k3),strs_p_G(1+((tnodes-1)*50/100),7,k3)
    enddo
    
    DEALLOCATE (disp_p,strs_p_G,strn_p_G,strs_p_Le,strs_p_Ls,post_p)
       
    return
    end subroutine

!----------------------------------------------------------------------------
!      FUNCTION post_disp_strs_strn FOR calculating disp, strains and stress
!----------------------------------------------------------------------------
	subroutine post_disp_strs_strn(xm1,zm1,alpha,beta,ele_cs,xyz_post,u_G,eps_G,sig_G,sig_Le,sig_Ls,nodeval)        
    use var_inputdata
    use var_analysis
    use integ_1D_beam
    use integ_2D_cs_TE
    use read_input_data
    use mod_math
    
    implicit none
    integer::i15,j15,o15,cn,ele_cs,i,ele_beam,count
    real*8::alpha,beta,eta,unode(3),uvw_xyz(9),epsnode(6),signode(6),signode_Le(6),signode_Ls(6),etainc,epsnode_mat(3,3),signode_mat(3,3),signode_Le_mat(3,3),signode_Ls_mat(3,3)
    real*8::eta_phy,xm1,zm1
    real*8,ALLOCATABLE, DIMENSION (:)::nodeval
    real*8,allocatable,dimension(:,:)::eps_G,sig_G,sig_Le,sig_Ls,u_G,xyz_post
        
    etainc=2.0d0/(nne-1)
    nodeval=0.0d0
    u_G=0.0d0
    eps_G=0.0d0
    sig_G=0.0d0
    sig_Le=0.0d0
    sig_Ls=0.0d0
    xyz_post=0.0d0
    
    call expfun(xm1,zm1)
    do ele_beam=1,ne
        do o15=1,nne
            cn=o15+(ele_beam-1)*(nne-1)
            eta=-1.0d0+((o15-1)*etainc) 
            eta_phy=((y((nne-1)*ele_beam+1)-y((nne-1)*(ele_beam-1)+1))/2.0d0)*(eta+1.0d0)+y((nne-1)*(ele_beam-1)+1)
            call shape1(eta)
            call mat_stiff(ele_cs,eta_phy)
            unode=0.0d0
            epsnode=0.0d0   
            uvw_xyz=0.0d0
        
            do i15=1,nne
                do j15=1,mexp
                    count=(3*mexp*(ele_beam-1)*(nne-1))+(3*mexp*(i15-1))+(3*j15)
                    unode(1)=unode(1)+uvect(count-2)*ft(j15)*sf(i15)
                    unode(2)=unode(2)+uvect(count-1)*ft(j15)*sf(i15)
                    unode(3)=unode(3)+uvect(count)*ft(j15)*sf(i15)
                    
                    epsnode(1)=epsnode(1)+ftx(j15)*sf(i15)*uvect(count-2)
                    epsnode(2)=epsnode(2)+ft(j15)*dsf(i15)*jacobinv(ele_beam)*uvect(count-1)
                    epsnode(3)=epsnode(3)+ftz(j15)*sf(i15)*uvect(count)
                    epsnode(4)=epsnode(4)+(ftz(j15)*sf(i15)*uvect(count-1))+(ft(j15)*dsf(i15)*jacobinv(ele_beam)*uvect(count))
                    epsnode(5)=epsnode(5)+(ftz(j15)*sf(i15)*uvect(count-2))+(ftx(j15)*sf(i15)*uvect(count))
                    epsnode(6)=epsnode(6)+(ft(j15)*dsf(i15)*jacobinv(ele_beam)*uvect(count-2))+(ftx(j15)*sf(i15)*uvect(count-1)) 
                enddo
            enddo
            nodeval(cn)=nodeval(cn)+1
            
            !!!!!!!!!!!!!!! compute material stiffness and local cys inverse !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call matinv3(loc_e,loc_e_inv)
            call matinv3(loc_struc,loc_struc_inv)
            
            !!!!!!!!!!!!!!!!!!!! Global stress vector and tensor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            signode=matmul(trcmat,epsnode)
            signode_mat(1,1)=signode(1);    signode_mat(1,2)=signode(6);    signode_mat(1,3)=signode(5)
            signode_mat(2,1)=signode(6);    signode_mat(2,2)=signode(2);    signode_mat(2,3)=signode(4)
            signode_mat(3,1)=signode(5);    signode_mat(3,2)=signode(4);    signode_mat(3,3)=signode(3)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!! Local stress tensor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!! Local csys wrt to material csys!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            signode_Le_mat=matmul(transpose(loc_e_inv),matmul(signode_mat,loc_e_inv))
            signode_Le(1)=signode_Le_mat(1,1);  signode_Le(2)=signode_Le_mat(2,2);  signode_Le(3)=signode_Le_mat(3,3)
            signode_Le(4)=signode_Le_mat(2,3);  signode_Le(5)=signode_Le_mat(1,3);  signode_Le(6)=signode_Le_mat(1,2)
            !!!!!!!!!!!!!!!! Local csys wrt to structure csys !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            signode_Ls_mat=matmul(transpose(loc_struc_inv),matmul(signode_mat,loc_struc_inv))
            signode_Ls(1)=signode_Ls_mat(1,1);  signode_Ls(2)=signode_Ls_mat(2,2);    signode_Ls(3)=signode_Ls_mat(3,3)
            signode_Ls(4)=signode_Ls_mat(2,3);  signode_Ls(5)=signode_Ls_mat(1,3);    signode_Ls(6)=signode_Ls_mat(1,2)
            
            !!!!!!!!!!!!!!!!!! global displacements !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do i=1,3
                u_G(cn,i)=(((nodeval(cn)-1)*u_G(cn,i))+unode(i))/nodeval(cn)
            enddo
            do i=1,6
                eps_G(cn,i)=(((nodeval(cn)-1)*eps_G(cn,i))+epsnode(i))/nodeval(cn)
                sig_G(cn,i)=(((nodeval(cn)-1)*sig_G(cn,i))+signode(i))/nodeval(cn)
                sig_Le(cn,i)=(((nodeval(cn)-1)*sig_Le(cn,i))+signode_Le(i))/nodeval(cn)
                sig_Ls(cn,i)=(((nodeval(cn)-1)*sig_Ls(cn,i))+signode_Ls(i))/nodeval(cn)
            enddo
        enddo
    enddo

    return
    end subroutine   
    
    
end module
