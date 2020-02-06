module read_input_data
    contains
      
subroutine read_input   
use var_inputdata
    
    implicit none
!!!!!!!!!!!!!! internal variables !!!!!!!!!!!!!!!!!!!!!
    real*8,allocatable,dimension(:,:)   ::elas_constant                 ! elastic constants

    integer                             ::i,j,k,ij,ik,ijk
    character*1                         ::symbol,mark
     
    real*8,ALLOCATABLE,DIMENSION(:)     ::anod,xloc,yloc,zloc           
    integer,ALLOCATABLE,DIMENSION(:,:)  ::elem1con,elecon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    symbol='!'
999  format(A1)
!!!!!!!!!!!!!!!!!!!!!!!!!! Read inputs !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

101 read(1,999)mark
    if(mark.eq.symbol) go to 101
    read(1,*)nmat,nlayer
    allocate(matin(nlayer),angle(nlayer,2),elas_constant(nmat,9))
     
102  read(1,999)mark
     if(mark.eq.symbol) go to 102
     do k=1,nlayer                              
 	    read(1,*)matin(k),angle(k,1),angle(k,2)
     enddo
     
103  read(1,999)mark
     if(mark.eq.symbol) go to 103
     do i=1,nmat
        read(1,*)(elas_constant(i,j),j=1,9)
     enddo 
     
     call inp1(elas_constant)                   ! call inp1 subroutine for initial calculation of C matrix (defined in read_input.f90) 
     deallocate(elas_constant)
     
104  read(1,999)mark
     if(mark.eq.symbol) go to 104
     read(1,*)nexp,agp                          ! nexp: order of expansion; agp: number of gauss points (for cross-section integration)
     
105  read(1,999)mark
     if(mark.eq.symbol) go to 105
     read(1,*)ne,nne,beam_mesh
     
    SELECT CASE(nne)
        CASE (2)
        ngp=2                                   ! ngp: number of gauss points (for 1D integration along the beam)
        CASE (3)
        ngp=3
        CASE (4)
        ngp=4
    END SELECT
    
106  read(1,999)mark
     if(mark.eq.symbol) go to 106
     read(1,*)ylengthi,ylengthf
     
     call inp2                                  ! call inp2 subroutine for the distribution of beam nodes (defined in read_input.f90) 
   
107  read(1,999)mark
     if(mark.eq.symbol) go to 107
     read(1,*)phynormesh
     
108  read(1,999)mark
     if(mark.eq.symbol) go to 108
     read(1,*)ps_cond                           ! plane strain condition in the width direction
     
     mexp=(nexp+1)*(nexp+2)/2                   ! mexp: number of terms in the cross-section expansion
     tdof=3*mexp*tnodes                         ! tnodes: total number of beam nodes (calcluated in inp2 subroutine)
     
     write(11,*)'CS expansion order and terms'
     write(11,*)nexp,mexp 
     write(11,*)'beam nodes'
     write(11,*)tnodes
     write(11,*)'total dof'
     write(11,*)tdof 
     write(11,*) ' Node no    Y-cord'
     do i=1,tnodes
        write(11,10)i,y(i)
     enddo
10   format(2x,i4,5x,f21.18)
     
    !!!!!!!!!!!!!!!!!!!!!!!!!!!! area mesh input !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
     read(2,999)mark
     read(2,*)anode,ane

     allocate(cscord(anode,3),elemcon(ane,4+1),layerno(ane,2))
     allocate(anod(anode),xloc(anode),yloc(anode),zloc(anode))
     allocate (elecon(ane,4+5),elem1con(ane,4+1))
     cscord=0.0d0                               ! cscord array: x,y,z coordinates of the cross-section nodes (at y=0) 
     elemcon=0                                  ! elemcon array: element connectivity
     elecon=0                                   ! elecon array: element connectivity (as read from the the input file)
     layerno=0                                  ! layerno array: layer number for the corresponding element
     
     read(2,999)mark
     do i=1,anode
        read(2,*)anod(i),xloc(i),yloc(i),zloc(i)
     enddo
     cscord(1:anode,1)=anod(:)
     cscord(1:anode,2)=xloc(:)
     cscord(1:anode,3)=zloc(:)
    
     read(2,*)mark
     do i=1,ane
        read(2,*)(elecon(i,j),j=1,4+5)
     enddo
     elem1con(:,1)=elecon(:,1)
     elem1con(:,2:4+1)=elecon(:,1+5:4+5)
     layerno(1:ane,1)=elecon(:,1)
     layerno(1:ane,2)=elecon(:,4)
     call cwacw(elem1con)                   ! cwacw subroutine for getting anti-clocklwise element connectivity (defined in read_input.f90) 
     
     deallocate (elecon,elem1con,xloc,yloc,zloc,anod)
!!!!!!!!!!!!!!!!!!!!!!!!!! Load and BC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
109  read(4,999)mark
     if(mark.eq.symbol) go to 109
     read(4,*)np
     
     allocate(loadinp(7,np))
      
110  read(4,999)mark
     if(mark.eq.symbol) go to 110
     do i=1,np
        read(4,*)(loadinp(j,i),j=1,7)
     enddo
     
111  read(4,999)mark
     if(mark.eq.symbol) go to 111
     read(4,*)num_bcinp 
     allocate(node_bc(num_bcinp,3))
     
112  read(4,999)mark
     if(mark.eq.symbol) go to 112
     do i=1,num_bcinp
        read(4,*)node_bc(i,1),node_bc(i,2),node_bc(i,3)
     enddo
     
     allocate(cs_node_bc(num_bcinp,nexp+1))
     cs_node_bc=0
     cs_node_bc(:,1)=node_bc(:,2)
     cs_node_bc(:,2)=node_bc(:,3)
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
end subroutine


!----------------------------------------------------------------
!      subroutine inp2 FOR beam NODE COORDINATES
!----------------------------------------------------------------
subroutine inp2
    use var_inputdata
    use var_analysis
    implicit none
    Integer::ii,jj,el,tnds,kk
    real*8::tl,dbn,lpe,pi,y1,y2
    
    allocate(jacob(ne),jacobinv(ne))
    jacob=0.0d0                         ! jacob array: 1D jacobian of each beam element
    jacobinv=0.0d0                      ! jacobinv array: 1D jacobian inverse for each beam element
    pi=4.0*ATAN(1.0)                    ! pi (3.14159265 rad)
        
    tnodes=((nne-1)*ne)+1               ! tnodes: total number of beam nodes
    allocate(y(tnodes))                 ! y array: y-cordinate of each beam node
        
    y=0.0d0                             ! initialise y array
    y(1)=ylengthi
        
    SELECT CASE (beam_mesh)
    CASE ('CHEB','Cheb','cheb')                 ! chebyshev biased distribution
        tl=ylengthf-ylengthi
        do ii=2,tnodes-1
            y(ii)=(tl/2.0)-(0.5d0*tl*cos(pi*(ii-0.5d0)/tnodes))
        enddo
        y(tnodes)=ylengthf
            
    CASE ('UNIF','Unif','unif')                 ! uniformly distributed
        tl=ylengthf-ylengthi
        lpe=tl/ne            
        dbn=lpe/(nne-1)    
        do ii=2,tnodes-1
            y(ii)=y(ii-1)+dbn
        enddo
        y(tnodes)=ylengthf
        
    CASE ('MANU','Manu','manu')                 ! manual (read coordinates from a file)
        open(6,file="input/beam_yinp.dat")
        do ii=1,tnodes
            read(6,*)y(ii)
        enddo
        close(6)   
        
    END SELECT
    
    do kk=1,ne
        ii=((nne-1)*(kk-1))+1
        jj=((nne-1)*(kk))+1
        jacob(kk)=(y(jj)-y(ii))/2.0d0
        jacobinv(kk)=1.0d0/jacob(kk)
    enddo
    

    return
end subroutine
    
    
!------------------------------------------------------------------------------
!      Definition of subroutine inp1 Material Properties
!------------------------------------------------------------------------------
    subroutine inp1(elas_constant)
    use var_inputdata
    use mod_math
    Implicit none
    integer::i
    real*8::E1,E2,E3,G23,G13,G12,nu12,nu13,nu23,delta,nu21,nu31,nu32
    real*8,allocatable,dimension(:,:)::elas_constant
    allocate(cmat(6,6,nmat))
    cmat=0.0
    
    do i=1,nmat
        E1=elas_constant(i,1)                   ! Ex
        E2=elas_constant(i,2)                   ! Ey
        E3=elas_constant(i,3)                   ! Ez
        G23=elas_constant(i,4)                  ! Gyz
        G13=elas_constant(i,5)                  ! Gxz
        G12=elas_constant(i,6)                  ! Gxy
        nu12=elas_constant(i,7)                 ! NUxy
        nu13=elas_constant(i,8)                 ! NUxz
        nu23=elas_constant(i,9)                 ! NUyz

        nu21=nu12*E2/E1                         ! NUyx
	    nu31=nu13*E3/E1                         ! NUzx
	    nu32=nu23*E3/E2                         ! NUzy
    
	    delta=(1.0-(nu12*nu21)-(nu23*nu32)-(nu31*nu13)-(2.*nu21*nu13*nu32))
        cmat(1,1,i) = (1.0 - nu23*nu32)*E1/delta
        cmat(2,2,i) = (1.0 - nu13*nu31)*E2/delta
        cmat(3,3,i) = (1.0 - nu12*nu21)*E3/delta
        cmat(1,2,i) = (nu12 + nu32*nu13)*E2/delta
        cmat(2,1,i) = (nu21 + nu23*nu31)*E1/delta
        cmat(1,3,i) = (nu13 + nu12*nu23)*E3/delta
        cmat(3,1,i) = (nu31 + nu21*nu32)*E1/delta
        cmat(2,3,i) = (nu23 + nu21*nu13)*E3/delta
        cmat(3,2,i) = (nu32 + nu12*nu31)*E2/delta
    
        if (E1.eq.E2.and.E2.eq.E3) then                     ! isotropic material
            cmat(4,4,i) = (E1/(2.0d0*(1+nu12)))
            cmat(5,5,i) = (E1/(2.0d0*(1+nu12)))
            cmat(6,6,i) = (E1/(2.0d0*(1+nu12)))
        else                                                ! orthotropic material
            cmat(4,4,i) = G23
            cmat(5,5,i) = G13
            cmat(6,6,i) = G12
        endif
    enddo
    
    return
    end subroutine

!------------------------------------------------------------------------------
!      Definition of subroutine mat_stiff Material Properties
!------------------------------------------------------------------------------
    subroutine mat_stiff(ele,yloc)
    use var_inputdata
    use mod_math
    Implicit none
    integer::mat,lay,ele,i30,j30
    real*8::ang,yloc,roty
    real*8,allocatable,dimension(:)::V,N,P,T,Tcos,Tsin
    real*8,allocatable,dimension(:,:)::cmat1,Qrot,lmn,lmn1
    allocate(cmat1(6,6),Qrot(6,6),lmn(5,5),lmn1(3,3))
    allocate(V(3),N(3),T(3),P(3),Tcos(3),Tsin(3))
    
    lay=layerno(ele,2)                      ! layer number
    mat=matin(lay)                          ! material number
    ang=(2*(angle(lay,1)-angle(lay,2))*abs(yloc-0.5d0*(y(tnodes)-y(1)))/(y(tnodes)-y(1)))+angle(lay,2)  ! ang: fibre orientation angle at a particular ylocation (yloc)
    cmat1(:,:)=cmat(:,:,mat)                ! C matrix for material mat
    
    V(1)=0.0d0; V(2)=1.0d0; V(3)=0.0d0      ! y-direction (axial)
    N(1)=0.0d0; N(2)=0.0d0; N(3)=1.0d0      ! z-direction (normal)
    P(1)=V(2)*N(3)-N(2)*V(3); P(2)=V(3)*N(1)-N(3)*V(1); P(3)=V(1)*N(2)-N(1)*V(2)    ! x-direction (tangent)
    
    T(1)=N(2)*P(3)-P(2)*N(3); T(2)=N(3)*P(1)-P(3)*N(1); T(3)=N(1)*P(2)-P(1)*N(2)    ! axial
    Tsin(1)=T(1)*sind(-ang); Tsin(2)=T(2)*sind(-ang); Tsin(3)=T(3)*sind(-ang)       
    Tcos(1)=T(1)*cosd(-ang); Tcos(2)=T(2)*cosd(-ang); Tcos(3)=T(3)*cosd(-ang)
!!!!!!!!!!!!! ang is the fibre angle with respect to y axis (beam axis) (+ve input means anti-cs) !!!!!!!!!!!!!!!!!!!!!!!!
    !!! loc_e is the coordinate system such that y-axis is always along the fibre direction
    !!! It can be called as a material csys, Z is normal, y is along fibre, x is tangent. 
    !!! loc_e will be diff for laminates with diff fibre orientation
    !!! loc_struc is another local csys which is only depends on the structure orientation, which is used in post-processing
    !!!!!!!!!!!!!!!!!! local x, y and z coordinate !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    loc_e(1,1)=(Tcos(2)*N(3)-Tcos(3)*N(2))-Tsin(1)
    loc_e(1,2)=(Tcos(3)*N(1)-Tcos(1)*N(3))-Tsin(2)
    loc_e(1,3)=(Tcos(1)*N(2)-Tcos(2)*N(1))-Tsin(3)
    
    loc_e(2,1)=Tsin(2)*N(3)-Tsin(3)*N(2)+Tcos(1)
    loc_e(2,2)=Tsin(3)*N(1)-Tsin(1)*N(3)+Tcos(2)
    loc_e(2,3)=Tsin(1)*N(2)-Tsin(2)*N(1)+Tcos(3) 
    
    loc_e(3,1)=N(1)
    loc_e(3,2)=N(2)
    loc_e(3,3)=N(3)
    
    loc_struc(1,1)=P(1)
    loc_struc(1,2)=P(2)
    loc_struc(1,3)=P(3)
    
    loc_struc(2,1)=V(1)
    loc_struc(2,2)=V(2)
    loc_struc(2,3)=V(3) 
    
    loc_struc(3,1)=N(1)
    loc_struc(3,2)=N(2)
    loc_struc(3,3)=N(3)
    
    !!!!!!!!!!!!!!!!!! global x, y and z coordinate !!!!!!!!!!!!!!!!!!!!!
    glb_e(1,1)=1.0d0; glb_e(1,2)=0.0d0; glb_e(1,3)=0.0d0
    glb_e(2,1)=0.0d0; glb_e(2,2)=1.0d0; glb_e(2,3)=0.0d0
    glb_e(3,1)=0.0d0; glb_e(3,2)=0.0d0; glb_e(3,3)=1.0d0
    
    lmn1=matmul(loc_e,transpose(glb_e))
    lmn(1:3,1:3)=lmn1(:,:)
    lmn(1:3,4)=lmn1(1:3,1)
    lmn(1:3,5)=lmn1(1:3,2)
    lmn(4:5,1:5)=lmn(1:2,1:5)
    
    do i30=1,3
        do j30=1,3
            Qrot(i30,j30)=lmn(i30,j30)*lmn(i30,j30)
        enddo
    enddo
   
    do i30=1,3
        do j30=4,6
            Qrot(i30,j30)=lmn(i30,j30-2)*lmn(i30,j30-1)
        enddo
    enddo
   
    do i30=4,6
        do j30=1,3
            Qrot(i30,j30)=2*lmn(i30-2,j30)*lmn(i30-1,j30)
        enddo
    enddo
    
    do i30=4,6
        do j30=4,6
            Qrot(i30,j30)=lmn(i30-2,j30-2)*lmn(i30-1,j30-1)+lmn(i30-1,j30-2)*lmn(i30-2,j30-1)
        enddo
    enddo
        
    trcmat=matmul(transpose(Qrot),matmul(cmat1,Qrot))

!!!!!!!!!Plane strain condition!!!!!!!!!!!
    
    if (ps_cond.eq.1) then              ! refer appendix A in https://doi.org/10.1016/j.compositesb.2018.08.127
        trcmat(1,2:6)=0.0d0
        trcmat(5,1:4)=0.0d0
        trcmat(5,6)=0.0d0
        trcmat(6,1:5)=0.0d0
        trcmat(2:6,1)=0.0d0
        trcmat(1:4,5)=0.0d0
        trcmat(6,5)=0.0d0
        trcmat(1:5,6)=0.0d0
        trcmat(1:2,6)=0.0d0
        trcmat(6,1:2)=0.0d0
        trcmat(4,5)=0.0d0
        trcmat(5,4)=0.0d0
    endif
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    
    deallocate(cmat1,Qrot,lmn,lmn1,V,N,T,P,Tcos,Tsin)
    
    return
    end subroutine

!----------------------------------------------------------------
!      FUNCTION cwacw for getting anti-clocklwise connectivity
!----------------------------------------------------------------
 subroutine cwacw(elem1con)
  use var_inputdata
    implicit none
    Integer::i,ii,nodnm,jj
    real*8::cross,x(4),z(4)
    integer,allocatable::elem1con(:,:)
            
    do jj=1,ane
        do ii=1,4
            nodnm=elem1con(jj,ii+1)
            x(ii)=cscord(nodnm,2)
            z(ii)=cscord(nodnm,3)
        enddo
        
        i=1
12345    cross=((x(i+1)-x(i))*(z(i+2)-z(i+1)))-((x(i+2)-x(i+1))*(z(i+1)-z(i)))
        if (cross.eq.0.0) then
            i=i+1
            goto 12345
            
            ! clockwise
        else if (cross.lt.0.0) then
                elemcon(jj,1)=elem1con(jj,1)
                elemcon(jj,2)=elem1con(jj,2)
                elemcon(jj,3)=elem1con(jj,5)
                elemcon(jj,4)=elem1con(jj,4)
                elemcon(jj,5)=elem1con(jj,3)
            ! anti-clockwise
        else if (cross.gt.0.0) then                             
            elemcon(jj,:)=elem1con(jj,:)
        endif

     enddo
  return
end subroutine
    
end module
    