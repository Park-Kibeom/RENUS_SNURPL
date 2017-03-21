subroutine setls2g(iftran)
!construct lineary system to solve 2-group problem.
!obtain coupling coeffs adn 2g matrix using d-tilda
    use const
    use cmfd2g
    use geom,       only : hmesh,volnode,lsfc    
    use xsec,       only : xschi,xsnf,xst,xss
    use xsec,       only : xsr   ! 2013_07_22 . scb
    use tran,       only : betap,betat,rvelotm, &
                            xschip ! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
    use tranxsec,   only : xschid
    !use cmfdmg,     only : dhatzf ! 2013_07_22 . scb for debugging 
    implicit none

    logical      :: iftran

    integer                 :: l,k,m2,mm,idir
    real                    :: vol,offdiag,wnp1
    real                    :: area(4), areaz(2), chi(ng2)

    do k=1,nz
        do l=1,nxy
            
!           determine the area of surfaces at coarse meshes that is normal to directions
            area(1)=hmesh(YDIR,l,k)*hmesh(ZDIR,l,k)
            area(2)=area(1)
            area(3)=hmesh(XDIR,l,k)*hmesh(ZDIR,l,k)
            area(4)=area(3)
            areaz(1)=hmesh(XDIR,l,k)*hmesh(YDIR,l,k)
            areaz(2)=areaz(1)

            chi(:)=0
            if(iftran) then 
                wnp1=betap(l,k)+betat(l,k)-ONE
!                chi(:)=(1-betat(l,k))*xschi(:,l,k)+wnp1*xschid(:,l,k)    ! added in ARTOS ver. 0.2 . 2012_08_06 by SCB    
                chi(:)=(1-betat(l,k))*xschip(:,l,k)+wnp1*xschid(:,l,k)    ! added in ARTOS ver. 0.2 . 2012_08_06 by SCB     ! mhtr
                continue
            else
                chi(:)=xschi(:,l,k)
            endif          

            vol=volnode(l,k)
            do m2=1,ng2
                mm=indm24(m2,m2)

                af(m2,l,k) = xsnf(m2,l,k) * vol
                !am(mm,l,k) = xst(m2,l,k) * vol - af(m2,l,k)*chi(m2)*reigvs
                am(mm,l,k) = xsr(m2,l,k) * vol - af(m2,l,k)*chi(m2)*reigvs   ! 2013_07_22 . scb
                !am(mm,l,k) = xsr(m2,l,k) * vol   ! 2013_07_22 . scb

                ! west
                idir=1
                offdiag=dtilr(m2,lsfc(idir,l),k)*area(idir)
                ccr(idir,m2,l,k)=-offdiag
                am(mm,l,k)=am(mm,l,k)+offdiag
                offdiag=dhatr(m2,lsfc(idir,l),k)*area(idir)
                ccr(idir,m2,l,k)=ccr(idir,m2,l,k)+offdiag
                am(mm,l,k)=am(mm,l,k)+offdiag

                ! east
                idir=2
                offdiag=dtilr(m2,lsfc(idir,l),k)*area(idir)
                ccr(idir,m2,l,k)=-offdiag
                am(mm,l,k)=am(mm,l,k)+offdiag
                offdiag=dhatr(m2,lsfc(idir,l),k)*area(idir)
                ccr(idir,m2,l,k)=ccr(idir,m2,l,k)-offdiag
                am(mm,l,k)=am(mm,l,k)-offdiag
                ! north
                idir=3
                offdiag=dtilr(m2,lsfc(idir,l),k)*area(idir)
                ccr(idir,m2,l,k)=-offdiag
                am(mm,l,k)=am(mm,l,k)+offdiag
                offdiag=dhatr(m2,lsfc(idir,l),k)*area(idir)
                ccr(idir,m2,l,k)=ccr(idir,m2,l,k)+offdiag
                am(mm,l,k)=am(mm,l,k)+offdiag
                ! south
                idir=4
                offdiag=dtilr(m2,lsfc(idir,l),k)*area(idir)
                ccr(idir,m2,l,k)=-offdiag
                am(mm,l,k)=am(mm,l,k)+offdiag
                offdiag=dhatr(m2,lsfc(idir,l),k)*area(idir)
                ccr(idir,m2,l,k)=ccr(idir,m2,l,k)-offdiag
                am(mm,l,k)=am(mm,l,k)-offdiag

                ! bottom
                idir=1
                offdiag = dtilz(m2,l,k+idir-1)
                offdiag=offdiag*areaz(idir)
                ccz(idir,m2,l,k)=-offdiag
                am(mm,l,k)=am(mm,l,k)+offdiag
                offdiag=dhatz(m2,l,k+idir-1)*areaz(idir)
                !offdiag=dhatzf(m2,l,k+idir-1)*areaz(idir)   ! 2013_07_22 . scb
                ccz(idir,m2,l,k)=ccz(idir,m2,l,k)+offdiag
                am(mm,l,k)=am(mm,l,k)+offdiag            
                ! up 
                idir=2
                offdiag = dtilz(m2,l,k+idir-1)
                offdiag=offdiag*areaz(idir)
                ccz(idir,m2,l,k)=-offdiag
                am(mm,l,k)=am(mm,l,k)+offdiag
                offdiag=dhatz(m2,l,k+idir-1)*areaz(idir)
                !offdiag=dhatzf(m2,l,k+idir-1)*areaz(idir)   ! 2013_07_22 . scb
                ccz(idir,m2,l,k)=ccz(idir,m2,l,k)-offdiag
                am(mm,l,k)=am(mm,l,k)-offdiag            
            enddo
            !
! 2013_07_22 . scb            
            am(indm24(1,2),l,k) = -xss(THERMAL,FAST,l,k)*vol - af(2,l,k)*chi(1)*reigvs
            am(indm24(2,1),l,k) = -xss(FAST,THERMAL,l,k)*vol - af(1,l,k)*chi(2)*reigvs
            !am(indm24(1,2),l,k) = -xss(THERMAL,FAST,l,k)*vol
            !am(indm24(2,1),l,k) = -xss(FAST,THERMAL,l,k)*vol
! added end            
        enddo
    enddo

    return    

end subroutine

! comment !!! fission & scattering terms are in LHS