module solution
    contains
    
!-------------------------------------------------------
!   FUNCTION lin_static_TE FOR linear static analysis
!-------------------------------------------------------
    subroutine lin_static_TE
    use var_inputdata
    use var_analysis
    use solver
  
    pvect_res(:)=fvect(:)
    call solve_disp_sparse
    uvect(:)=uvect(:)+del_uvect(:)
    call uvect_BC(uvect)

    return
    end subroutine

    
!----------------------------------------------------------------
!      FUNCTION uvect_BC FOR enforcing uvect to zero to model BCs
!----------------------------------------------------------------
	subroutine uvect_BC(vector)    
    use var_inputdata
    use var_analysis
    
    implicit none
    integer::row_nz,i,j,beam_nd,row_nnz,ii,ij,j1
    real*8,allocatable,dimension(:)::vector
    
    do i=1,num_bcinp
        beam_nd=node_bc(i,1)
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Fixed BC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (node_bc(i,2).eq.0.and.node_bc(i,3).eq.0) then           
           if (beam_nd.eq.1) then
               do j=1,3*mexp
                    vector(j)=0.0d0
               enddo 
           elseif (beam_nd.eq.tnodes) then  
               do j=(tdof-3*mexp)+1,tdof
                  vector(j)=0.0d0
               enddo
           else
               do j=3*mexp*(beam_nd-1)+1,3*mexp*beam_nd
                  vector(j)=0.0d0
               enddo
           endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Simply-supported BC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        elseif (node_bc(i,2).eq.100.and.node_bc(i,3).eq.100) then
            if (beam_nd.eq.1) then                                                      ! first beam node
               do j=1,3*mexp
                   if (mod(j-1,3).eq.0.or.mod(j,3).eq.0) then                           !!!!!! x and z constrained !!!!!!
                       vector(j)=0.0d0
                   endif
               enddo 
           elseif (beam_nd.eq.tnodes) then                                              ! last beam node
               do j=(tdof-3*mexp)+1,tdof
                    if (mod(j-1,3).eq.0.or.mod(j,3).eq.0) then                          !!!!!! x and z constrained !!!!!!
                      vector(j)=0.0d0
                    endif                    
                enddo
           else
               do j=3*(beam_nd-1)*mexp+1,3*beam_nd*mexp
                   if (mod(j+1,3).eq.0) then                                           !!!!!! y constrained !!!!!! 
                      vector(j)=0.0d0
                   endif
               enddo
           endif
        endif
    enddo
              
    return
    end subroutine
    
end module