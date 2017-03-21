subroutine calsol2drhs(k)
    use ppr
    use sfam,   only : reigv
    use geom,   only : ng,nxy
    use xsec,   only : xsnff,xssf,xssfs,xssfe,xschif
    implicit none

    integer,intent(in)  :: k
    
    integer             :: i,m,ms,l
    real                :: psi2d,sc2d(ng)
    
    do l=1,nxy
        do i=0,14
            psi2d=0.0
            sc2d=0.0
            do m=1,ng
                psi2d=psi2d+reigv*xsnff(m,l,k)*qf2d(i,l,k,m)
                do ms=xssfs(m,l,k),xssfe(m,l,k)
                    sc2d(m)=sc2d(m)+xssf(ms,m,l,k)*qf2d(i,l,k,ms)
                enddo
            enddo
            do m=1,ng
                qc2d(i,l,k,m)=xschif(m,l,k)*psi2d+sc2d(m)
            enddo
        enddo 

!       1-(0,1), 2-(0,2), 3-(0,3), 4-(0,4)
!       5-(1,0), 6-(2,0), 7-(3,0), 8-(4,0)
!       9-(1,1), 10-(1,2),11-(2,1),12-(2,2)
        do m=1,ng
            qc2d(0,l,k,m)=qc2d(0,l,k,m)-trlzcff(0,0,m,l,k)

            qc2d(1,l,k,m)=qc2d(1,l,k,m)-trlzcff(0,1,m,l,k)
            qc2d(2,l,k,m)=qc2d(2,l,k,m)-trlzcff(0,2,m,l,k)

            qc2d(5,l,k,m)=qc2d(5,l,k,m)-trlzcff(1,0,m,l,k)
            qc2d(6,l,k,m)=qc2d(6,l,k,m)-trlzcff(2,0,m,l,k)

            qc2d(9,l,k,m)=qc2d(9,l,k,m)-trlzcff(1,1,m,l,k)
            qc2d(10,l,k,m)=qc2d(10,l,k,m)-trlzcff(1,2,m,l,k)
            qc2d(11,l,k,m)=qc2d(11,l,k,m)-trlzcff(2,1,m,l,k)
            qc2d(12,l,k,m)=qc2d(12,l,k,m)-trlzcff(2,2,m,l,k)
        enddo ! m      
    enddo ! l
    return
end subroutine