subroutine calsol2d(k)
    use ppr
    use sfam,   only : jnet
    use geom,   only : ng,nxy,hmesh
    use xsec,   only : xstf, xsdf, xscdf
    implicit none
    
    integer                 :: k

    integer                 :: m,idir,l
    real,dimension(2)       :: alpha,jhr,jhl,jhs,jhd
    real                    :: p0yp1p1,p0ym1p1,px0p1p1, &
                               px0m1p1,pxyp1p1,pxym1m1, &
                               pxym1p1,pxyp1m1,         &
                               fm1m1,fm1p1,fp1p1,fp1m1, &
                               p1,p2,p3,p4
                               
    real                    :: cpc22qc2d12, chc24hc2d8,rxstf,temp
    real                    :: rh(NRDIR)
    
    
    do l=1,nxy
        rh(:) = 1./hmesh(XDIR:YDIR,l,k)
        do m=1,ng
            rxstf=1/xstf(m,l,k)

!           for pc2d(2,l,k,m) and pc2d(6,l,k,m)
            cpc22qc2d12=cpc22(l,k,m)*qc2d(12,l,k,m)

!           coefficients of particular solutions      
!           1-(0,1), 2-(0,2), 3-(0,3), 4-(0,4)
!           5-(1,0), 6-(2,0), 7-(3,0), 8-(4,0)
!           9-(1,1), 10-(1,2),11-(2,1),12-(2,2)
            pc2d(0,l,k,m)=cpc02(l,k,m)*(qc2d(2,l,k,m)+qc2d(6,l,k,m))        &
                         +cpc04(l,k,m)*(qc2d(4,l,k,m)+qc2d(8,l,k,m))        &
                         +cpc022(l,k,m)*qc2d(12,l,k,m)                      &
                         +rxstf*qc2d(0,l,k,m)      

            pc2d(1,l,k,m)=cpc11(l,k,m)*qc2d(3,l,k,m)                        &
                          +cpc12(l,k,m)*qc2d(11,l,k,m)                      &
                          +rxstf*qc2d(1,l,k,m)

            pc2d(2,l,k,m)=cpc21(l,k,m)*qc2d(4,l,k,m)                        &
                          +cpc22qc2d12                                      &
                          +rxstf*qc2d(2,l,k,m)

            pc2d(3,l,k,m)=rxstf*qc2d(3,l,k,m)

            pc2d(4,l,k,m)=rxstf*qc2d(4,l,k,m)

            pc2d(5,l,k,m)=cpc11(l,k,m)*qc2d(7,l,k,m)                        &
                          +cpc12(l,k,m)*qc2d(10,l,k,m)                      &
                          +rxstf*qc2d(5,l,k,m)

            pc2d(6,l,k,m)=cpc21(l,k,m)*qc2d(8,l,k,m)                        &
                          +cpc22qc2d12                                      &
                          +rxstf*qc2d(6,l,k,m)

            pc2d(7,l,k,m)=rxstf*qc2d(7,l,k,m)
            
            pc2d(8,l,k,m)=rxstf*qc2d(8,l,k,m)
            
            pc2d(9,l,k,m)=rxstf*qc2d(9,l,k,m)                               &
                          +60*xsdf(m,l,k)*rxstf*rxstf*rh(XDIR)*rh(XDIR)     &
                          *(qc2d(13,l,k,m)+qc2d(14,l,k,m))
                          
            pc2d(10,l,k,m)=rxstf*qc2d(10,l,k,m)

            pc2d(11,l,k,m)=rxstf*qc2d(11,l,k,m)

            pc2d(12,l,k,m)=rxstf*qc2d(12,l,k,m)

            if(term15) then
                pc2d(13,l,k,m)=rxstf*qc2d(13,l,k,m)
                pc2d(14,l,k,m)=rxstf*qc2d(14,l,k,m)
            endif

    !
    !       obtain corner points for homogenous boundary conditions
            p0yp1p1=+pc2d(1,l,k,m)+pc2d(2,l,k,m)+pc2d(3,l,k,m)+pc2d(4,l,k,m)
            p0ym1p1=-pc2d(1,l,k,m)+pc2d(2,l,k,m)-pc2d(3,l,k,m)+pc2d(4,l,k,m)
            px0p1p1=+pc2d(5,l,k,m)+pc2d(6,l,k,m)+pc2d(7,l,k,m)+pc2d(8,l,k,m)
            px0m1p1=-pc2d(5,l,k,m)+pc2d(6,l,k,m)-pc2d(7,l,k,m)+pc2d(8,l,k,m)

            pxyp1p1=+pc2d(9,l,k,m)+pc2d(10,l,k,m)+pc2d(11,l,k,m)+pc2d(12,l,k,m)+pc2d(13,l,k,m)+pc2d(14,l,k,m)
            pxym1m1=-pc2d(9,l,k,m)-pc2d(10,l,k,m)+pc2d(11,l,k,m)+pc2d(12,l,k,m)-pc2d(13,l,k,m)-pc2d(14,l,k,m)
            pxym1p1=-pc2d(9,l,k,m)+pc2d(10,l,k,m)-pc2d(11,l,k,m)+pc2d(12,l,k,m)-pc2d(13,l,k,m)-pc2d(14,l,k,m)
            pxyp1m1=+pc2d(9,l,k,m)-pc2d(10,l,k,m)-pc2d(11,l,k,m)+pc2d(12,l,k,m)+pc2d(13,l,k,m)+pc2d(14,l,k,m)

            fm1m1=phicorn(lcnw(l),k,m)/xscdf(1,m,l,k)-(pc2d(0,l,k,m)+p0ym1p1+px0m1p1+pxyp1m1) !x=-1,y=-1
            fp1m1=phicorn(lcne(l),k,m)/xscdf(2,m,l,k)-(pc2d(0,l,k,m)+p0ym1p1+px0p1p1+pxym1p1) !x=+1,y=-1
            fm1p1=phicorn(lcsw(l),k,m)/xscdf(3,m,l,k)-(pc2d(0,l,k,m)+p0yp1p1+px0m1p1+pxym1m1) !x=-1,y=+1
            fp1p1=phicorn(lcse(l),k,m)/xscdf(4,m,l,k)-(pc2d(0,l,k,m)+p0yp1p1+px0p1p1+pxyp1p1) !x=+1,y=+1
    !
            p1=0.25*(fp1p1-fm1p1+fp1m1-fm1m1)
            p2=0.25*(fp1p1-fm1p1-fp1m1+fm1m1)
            p3=0.25*(fp1p1+fm1p1-fp1m1-fm1m1)
            p4=0.25*(fp1p1+fm1p1+fp1m1+fm1m1)

    !       obtain jnet for homegeneous boundary condtions
            do idir=1,NRDIR
                alpha(idir)=-2*xsdf(m,l,k)*rh(idir)
            enddo
            jhr(1)=jnet(RIGHT,m,l,k,XDIR)-alpha(1)*(pc2d(5,l,k,m)+3*pc2d(6,l,k,m)+6*pc2d(7,l,k,m)+10*pc2d(8,l,k,m))
            jhl(1)=jnet(LEFT,m,l,k,XDIR)-alpha(1)*(pc2d(5,l,k,m)-3*pc2d(6,l,k,m)+6*pc2d(7,l,k,m)-10*pc2d(8,l,k,m))
            jhr(2)=jnet(RIGHT,m,l,k,YDIR)-alpha(2)*(pc2d(1,l,k,m)+3*pc2d(2,l,k,m)+6*pc2d(3,l,k,m)+10*pc2d(4,l,k,m))
            jhl(2)=jnet(LEFT,m,l,k,YDIR)-alpha(2)*(pc2d(1,l,k,m)-3*pc2d(2,l,k,m)+6*pc2d(3,l,k,m)-10*pc2d(4,l,k,m))

            temp=0.5/kappa(l,k,m)
            do idir=1,NRDIR
                jhs(idir)=temp*(jhr(idir)+jhl(idir))
                jhd(idir)=temp*(jhr(idir)-jhl(idir))
            enddo

    !       coefficients of homogeneous solutions    
            hc2d(6,l,k,m)=chc6(l,k,m)*p2
            hc2d(1,l,k,m)=chc13j(l,k,m)*jhs(1)+chc13p(l,k,m)*p1
            hc2d(3,l,k,m)=chc13j(l,k,m)*jhs(2)+chc13p(l,k,m)*p3
            hc2d(5,l,k,m)=chc57j(l,k,m)*jhs(1)+chc57p(l,k,m)*p1
            hc2d(7,l,k,m)=chc57j(l,k,m)*jhs(2)+chc57p(l,k,m)*p3
            hc2d(8,l,k,m)=chc8j(l,k,m)*(jhd(1)+jhd(2))+chc8p(l,k,m)*p4
            chc24hc2d8=chc24a(l,k,m)*hc2d(8,l,k,m)
            hc2d(2,l,k,m)=chc24j(l,k,m)*jhd(1)+chc24hc2d8
            hc2d(4,l,k,m)=chc24j(l,k,m)*jhd(2)+chc24hc2d8
        enddo
    enddo ! l
    return
end subroutine