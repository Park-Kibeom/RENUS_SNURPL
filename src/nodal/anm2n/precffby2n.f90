subroutine precffby2n
    use const
    use anm2n
    use geom,   only  : hmesh,ndir,ng,nxy,nz
    use xsec,   only  : xsadf,xsd
    implicit none
    
    real              :: amat(2,2),bmat(2,2),tmp,rh
    integer           :: sgn,idir,m,l,k
    
    do idir=1,ndir 
    do k=1,nz
    do l=1,nxy
!   construct partial linear systems for flux continuity
        amat(1,1)=xsadf(1,1,l,k)*fluxratio(1,l,k)
        amat(2,1)=xsadf(1,1,l,k)*fluxratio(2,l,k)
        amat(1,2)=xsadf(1,2,l,k)
        amat(2,2)=xsadf(1,2,l,k)

        fluxmat(1,1,l,k,idir) = amat(1,1)*snkp(l,k,idir)
        fluxmat(2,1,l,k,idir) = amat(2,1)*snmu(l,k,idir)
        fluxmat(1,2,l,k,idir) = amat(1,2)*snkp(l,k,idir)
        fluxmat(2,2,l,k,idir) = amat(2,2)*snmu(l,k,idir)

        bmat(1,1)             = amat(1,1)*cnkp(l,k,idir)
        bmat(2,1)             = amat(2,1)*cnmu(l,k,idir)
        bmat(1,2)             = amat(1,2)*cnkp(l,k,idir)
        bmat(2,2)             = amat(2,2)*cnmu(l,k,idir)

        fluxrhs(1,l,k,idir)   = bmat(1,1)*heven(1,l,k,idir) &
                               +bmat(2,1)*heven(2,l,k,idir) &
                               +xsadf(1,1,l,k)*(            &
                                    part0(1,l,k,idir)       &
                                   +part2(1,l,k,idir)       &
                               )
        fluxrhs(2,l,k,idir)   = bmat(1,2)*heven(1,l,k,idir) &
                               +bmat(2,2)*heven(2,l,k,idir) &
                               +xsadf(1,2,l,k)*(            &
                                    part0(2,l,k,idir)       &
                                   +part2(2,l,k,idir)       &
                               )
        
!       construct partial linear systems for current continuity
        sgn = 1
        if(kflag(l,k)) sgn = -1

        amat(1,1)=xsd(1,l,k)*fluxratio(1,l,k)
        amat(2,1)=xsd(1,l,k)*fluxratio(2,l,k)
        amat(1,2)=xsd(2,l,k)
        amat(2,2)=xsd(2,l,k)

        currmat(1,:,l,k,idir) = amat(1,:)*kp(l,k)*cnkp(l,k,idir) 
        currmat(2,:,l,k,idir) = amat(2,:)*mu(l,k)*cnmu(l,k,idir)

        bmat(1,:)             = amat(1,:)*kp(l,k)*snkp(l,k,idir)*sgn
        bmat(2,:)             = amat(2,:)*mu(l,k)*snmu(l,k,idir)

        currrhs(1,l,k,idir)   = bmat(1,1)*heven(1,l,k,idir) &
                               +bmat(2,1)*heven(2,l,k,idir) 
        currrhs(2,l,k,idir)   = bmat(1,2)*heven(1,l,k,idir) &
                               +bmat(2,2)*heven(2,l,k,idir) 

        rh = 1/hmesh(idir,l,k)
        do m=1,ng2
            tmp = -2*xsd(m,l,k)*rh
            curpartr(m,l,k,idir) = tmp*(part1(m,l,k,idir)   &
                                   +3*part2(m,l,k,idir))    
            curpartl(m,l,k,idir) = tmp*(part1(m,l,k,idir)   &
                                   -3*part2(m,l,k,idir))
        enddo
    enddo
    enddo
    enddo
end subroutine