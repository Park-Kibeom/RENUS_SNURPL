subroutine expflux13(k)
    use ppr
    use sfam,   only : jnet,phisfc,phif
    use geom,   only : ng, nxy, hmesh
    use xsec,   only : xsdf,xscdf
    implicit none
    
    integer,intent(in)      :: k
    integer                 :: m,l

    integer                 :: idir,io
    real                    :: alpha(2)                ! -2.*Diffusion Coeff/width
    real                    :: ps(2),pd(2)             ! suface phi sum, surface phi difference      
    real                    :: js(2),jd(2) 
    real                    :: val3,val4,val2,val1
    real                    :: p1,p2,p3,p4
    real                    :: temp,rxsdf
    
    do l=1,nxy
        do m=1,ng
            rxsdf=1/xsdf(m,l,k)

            do idir=1,NRDIR !i-direction
                alpha(idir)=-0.5*hmesh(idir,l,k)/xsdf(m,l,k)
                ps(idir)=phisfc(RIGHT,m,l,k,idir)+phisfc(LEFT,m,l,k,idir)
                pd(idir)=phisfc(RIGHT,m,l,k,idir)-phisfc(LEFT,m,l,k,idir)
                js(idir)=alpha(idir)*(jnet(RIGHT,m,l,k,idir)+jnet(LEFT,m,l,k,idir))
                jd(idir)=alpha(idir)*(jnet(RIGHT,m,l,k,idir)-jnet(LEFT,m,l,k,idir))
            enddo     
                                           
            val3=phicorn(lcse(l),k,m)/xscdf(4,m,l,k)            !  1(nw)-----2(ne)  
            val4=phicorn(lcsw(l),k,m)/xscdf(3,m,l,k)            !    |         | 
            val2=phicorn(lcne(l),k,m)/xscdf(2,m,l,k)            !    |   node  | 
            val1=phicorn(lcnw(l),k,m)/xscdf(1,m,l,k)            !    |         |
            p1=rfour*(val3-val4+val2-val1)                      !  4(sw)-----3(se)  
            p2=rfour*(val3-val4-val2+val1)
            p3=rfour*(val3+val4-val2-val1)
            p4=rfour*(val3+val4+val2+val1)

            qf2d(0,l,k,m)=phif(m,l,k)      
            do idir=1,NRDIR
                if(idir.eq.XDIR) then
                    io=5 ! i-order ! watch out for x and y terms
                else
                    io=1
                endif
                temp=ps(idir)-2*phif(m,l,k)
                qf2d(io+0,l,k,m)=r10*(-js(idir)+6*pd(idir))         !5- y^1, 1- x^1
                qf2d(io+1,l,k,m)=r14*(-jd(idir)+10*temp)          !6- y^2, 2- x^2
                qf2d(io+2,l,k,m)=r10*(+js(idir)-pd(idir))           !7- y^3, 3- x^3
                qf2d(io+3,l,k,m)=r14*(+jd(idir)-3*temp)           !8- y^4, 4- x^4
            enddo   
            qf2d(9,l,k,m)=p2                                  ! 9- x^1 * y^1
            qf2d(10,l,k,m)=p1-half*pd(1)                      !10- x^1 * y^2
            qf2d(11,l,k,m)=p3-half*pd(2)                      !11- x^2 * y^1
            qf2d(12,l,k,m)=phif(m,l,k)+p4-half*(ps(1)+ps(2))  !12- x^2 * y^2      
        enddo ! m
    enddo ! l
    return
end subroutine