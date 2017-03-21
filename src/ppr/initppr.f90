subroutine initppr(iftran)
    use ppr
    use geom,   only : ng,nxy,nz,hmesh
    use xsec,   only : xstf,xsdf
    implicit none
    
    logical             :: iftran
    integer             :: m,l,k

!   initialized normalized buckling-kappa for hx=hy
!   initialize problem dependent parameters
    do k=1,nz
    do l=1,nxy
    do m=1,ng
        kappa(l,k,m)=0.5*hmesh(XDIR,l,k)*sqrt(xstf(m,l,k)/xsdf(m,l,k))
        call calcffls(m,l,k,kappa(l,k,m))
        call calcffsol2d(m,l,k,kappa(l,k,m))
        call calcffjcorn(m,l,k,kappa(l,k,m))
    enddo
    enddo
    enddo
    
    call caltrlz(iftran)
    call initcorn

    return
end subroutine