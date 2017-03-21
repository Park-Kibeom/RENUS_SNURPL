subroutine setlsmg(iftran)
!construct linear system to solve multi group problem.
    use const
    use cmfdmg
    use geom,       only : hmesh,volnode,lsfc
    use xsec,       only : xstf
    use tran,       only : rveloftm
    use geom,       only : ng   ! 2013_07_22 . scb
    implicit none
!
    logical      :: iftran

    integer                 :: l,k,m,idir
    real                    :: offdiag
    real                    :: area(4), areaz(2)

    do k=1,nz
        do l=1,nxy
!determine the area of surfaces at coarse meshes that is normal to directions
            area(1)=hmesh(YDIR,l,k)*hmesh(ZDIR,l,k)
            area(2)=area(1)
            area(3)=hmesh(XDIR,l,k)*hmesh(ZDIR,l,k)
            area(4)=area(3)

            areaz(1)=hmesh(XDIR,l,k)*hmesh(YDIR,l,k)
            areaz(2)=areaz(1)   

            do m=1,ng           
                diagf(l,k,m)=xstf(m,l,k)*volnode(l,k)
                ! west
                idir=1
                offdiag=dtilrf(m,lsfc(idir,l),k)*area(idir)
                ccrf(idir,l,k,m)=-offdiag
                diagf(l,k,m)=diagf(l,k,m)+offdiag 
                offdiag=dhatrf(m,lsfc(idir,l),k)*area(idir)
                ccrf(idir,l,k,m)=ccrf(idir,l,k,m)+offdiag
                diagf(l,k,m)=diagf(l,k,m)+offdiag 
                ! east
                idir=2
                offdiag=dtilrf(m,lsfc(idir,l),k)*area(idir)
                ccrf(idir,l,k,m)=-offdiag
                diagf(l,k,m)=diagf(l,k,m)+offdiag 
                offdiag=dhatrf(m,lsfc(idir,l),k)*area(idir)
                ccrf(idir,l,k,m)=ccrf(idir,l,k,m)-offdiag
                diagf(l,k,m)=diagf(l,k,m)-offdiag 
                ! north
                idir=3
                offdiag=dtilrf(m,lsfc(idir,l),k)*area(idir)
                ccrf(idir,l,k,m)=-offdiag
                diagf(l,k,m)=diagf(l,k,m)+offdiag 
                offdiag=dhatrf(m,lsfc(idir,l),k)*area(idir)
                ccrf(idir,l,k,m)=ccrf(idir,l,k,m)+offdiag
                diagf(l,k,m)=diagf(l,k,m)+offdiag 
                ! south
                idir=4
                offdiag=dtilrf(m,lsfc(idir,l),k)*area(idir)
                ccrf(idir,l,k,m)=-offdiag
                diagf(l,k,m)=diagf(l,k,m)+offdiag 
                offdiag=dhatrf(m,lsfc(idir,l),k)*area(idir)
                ccrf(idir,l,k,m)=ccrf(idir,l,k,m)-offdiag
                diagf(l,k,m)=diagf(l,k,m)-offdiag 
                ! bottom
                idir=1
                offdiag = dtilzf(m,l,k+idir-1)
                offdiag=offdiag*areaz(idir)
                cczf(idir,l,k,m)=-offdiag
                diagf(l,k,m)=diagf(l,k,m)+offdiag 
                offdiag=dhatzf(m,l,k+idir-1)*areaz(idir)
                cczf(idir,l,k,m)=cczf(idir,l,k,m)+offdiag
                diagf(l,k,m)=diagf(l,k,m)+offdiag 
                ! up
                idir=2
                offdiag = dtilzf(m,l,k+idir-1)
                offdiag=offdiag*areaz(idir)
                cczf(idir,l,k,m)=-offdiag
                diagf(l,k,m)=diagf(l,k,m)+offdiag 
                offdiag=dhatzf(m,l,k+idir-1)*areaz(idir)
                cczf(idir,l,k,m)=cczf(idir,l,k,m)-offdiag
                diagf(l,k,m)=diagf(l,k,m)-offdiag 
            enddo
        enddo
    enddo

    return
end subroutine
