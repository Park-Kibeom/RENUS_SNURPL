subroutine mallocppr(npinl)
    use allocs
    use ppr
    use geom,   only : ng, nx, ny, nxy, nz, nxs,nxe,nxya
    implicit none
    
    integer                 :: npinl
    integer                 :: j,nupper,ndown
    
    npin    = npinl
    ncorn   = nxe(1)-nxs(1)+2
    ndown   = ncorn
    do j=1+1,ny
        nupper  = ndown
        ndown   = nxe(j)-nxs(j)+2
        ncorn   = ncorn+max(nupper,ndown)
    enddo
    ncorn = ncorn +ndown
    
    call dmalloc0(nodec,0,nx+2,0,ny+2)
    call dmalloc(nxcs,ny+1)
    call dmalloc(nxce,ny+1)
    call dmalloc0(lctox,0,ncorn)
    call dmalloc0(lctoy,0,ncorn)
    call dmalloc(lcc,8,ncorn)            !1-w,2-n,3-e,4-s,5-nw,6-ne,7-se,8-sw
    call dmalloc(lcn,4,ncorn)            !1-nw,2-ne,3-se,4-sw
    call dmalloc(lcnw,nxy)               !corners belonging to a node -nw
    call dmalloc(lcsw,nxy)               !corners belonging to a node -sw
    call dmalloc(lcne,nxy)               !corners belonging to a node -ne
    call dmalloc(lcse,nxy)               !corners belonging to a node -se

    call dmalloc(phicorn,ncorn,nz,ng)
    call dmalloc(phicorn0,ncorn,nz,ng)
    call dmalloc(powvalr,npin,npin,nxya)
    call dmalloc(powval1,npin,npin)
    call dmalloc(powval,npin,npin,ng)
    call dmalloc0(pppeak,1,nxya,0,nz)
    call dmalloc(powtot,ng)
    call dmalloc(avgjnetz,ng,nxy,nz)
    call dmalloc0(trlzcff,0,2,0,2,1,ng,1,nxy,1,nz)  

! buckling
    call dmalloc(kappa,nxy,nz,ng)
! pinprcoeffc
    call dmalloc0(qf2d,0,14,1,nxy,1,nz,1,ng) 
    call dmalloc0(qc2d,0,14,1,nxy,1,nz,1,ng)
    call dmalloc0(pc2d,0,14,1,nxy,1,nz,1,ng)
    call dmalloc(hc2d,8,nxy,nz,ng)
    call dmalloc(jcornx,4,nxy,nz,ng)
    call dmalloc(jcorny,4,nxy,nz,ng)
! pinprcoeff_ls
    call dmalloc(clsqf01,nxy,nz,ng)
    call dmalloc(clsqf02,nxy,nz,ng)
    call dmalloc(clsqf11,nxy,nz,ng)
    call dmalloc(clsqf12,nxy,nz,ng)
    call dmalloc(clsqf21,nxy,nz,ng)
    call dmalloc(clsqf22,nxy,nz,ng)
    call dmalloc(clsqf31,nxy,nz,ng)
    call dmalloc(clsqf32,nxy,nz,ng)
    call dmalloc(clsqf41,nxy,nz,ng)
    call dmalloc(clsqf42,nxy,nz,ng)
    call dmalloc(clsqfx1y1,nxy,nz,ng)
    call dmalloc(clsqf1221,nxy,nz,ng)
    call dmalloc(clsqf1331,nxy,nz,ng)
    call dmalloc(clsqfx2y2,nxy,nz,ng)
! pinprcoeff_gensol2d
    call dmalloc(cpc02,nxy,nz,ng) 
    call dmalloc(cpc04,nxy,nz,ng) 
    call dmalloc(cpc022,nxy,nz,ng)
    call dmalloc(cpc11,nxy,nz,ng) 
    call dmalloc(cpc12,nxy,nz,ng) 
    call dmalloc(cpc21,nxy,nz,ng) 
    call dmalloc(cpc22,nxy,nz,ng) 
    call dmalloc(chc6,nxy,nz,ng)  
    call dmalloc(chc13j,nxy,nz,ng)
    call dmalloc(chc13p,nxy,nz,ng)
    call dmalloc(chc57j,nxy,nz,ng)
    call dmalloc(chc57p,nxy,nz,ng)
    call dmalloc(chc8j,nxy,nz,ng) 
    call dmalloc(chc8p,nxy,nz,ng) 
    call dmalloc(chc24j,nxy,nz,ng)
    call dmalloc(chc24a,nxy,nz,ng)
! pinprcoeff_jcorner
    call dmalloc(cpjxh1,4,nxy,nz,ng)      
    call dmalloc(cpjxh2,4,nxy,nz,ng)      
    call dmalloc(cpjxh5,4,nxy,nz,ng)      
    call dmalloc(cpjxh6,4,nxy,nz,ng)      
    call dmalloc(cpjxh7,4,nxy,nz,ng)      
    call dmalloc(cpjxh8,4,nxy,nz,ng)      
    call dmalloc(cpjxp6,4,nxy,nz,ng)      
    call dmalloc(cpjxp7,4,nxy,nz,ng)      
    call dmalloc(cpjxp8,4,nxy,nz,ng)      
    call dmalloc(cpjxp9,4,nxy,nz,ng)      
    call dmalloc(cpjxp11,4,nxy,nz,ng)     
    call dmalloc(cpjxp12,4,nxy,nz,ng)     
    call dmalloc(cpjxp13,4,nxy,nz,ng)     
    call dmalloc(cpjxp14,4,nxy,nz,ng)     
    call dmalloc(cpjyh3,4,nxy,nz,ng)      
    call dmalloc(cpjyh4,4,nxy,nz,ng)      
    call dmalloc(cpjyh5,4,nxy,nz,ng)      
    call dmalloc(cpjyh6,4,nxy,nz,ng)      
    call dmalloc(cpjyh7,4,nxy,nz,ng)      
    call dmalloc(cpjyh8,4,nxy,nz,ng)      
    call dmalloc(cpjyp2,4,nxy,nz,ng)      
    call dmalloc(cpjyp3,4,nxy,nz,ng)      
    call dmalloc(cpjyp4,4,nxy,nz,ng)      
    call dmalloc(cpjyp9,4,nxy,nz,ng)      
    call dmalloc(cpjyp10,4,nxy,nz,ng)     
    call dmalloc(cpjyp12,4,nxy,nz,ng)           
    call dmalloc(cpjyp13,4,nxy,nz,ng)           
    call dmalloc(cpjyp14,4,nxy,nz,ng)           

    return    
end subroutine
