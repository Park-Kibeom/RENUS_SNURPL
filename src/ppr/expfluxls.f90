subroutine expfluxls(k)
    use ppr
    use sfam,   only : phif
    use geom,   only : nxy,ng
    implicit none
    
    integer             :: k

    integer             :: m,l
    real                :: c22hc2d8,c42hc2d8

    do l=1,nxy
    do m=1,ng
        qf2d(0,l,k,m)   =clsqf01(l,k,m)*(hc2d(2,l,k,m)+hc2d(4,l,k,m))   &
                        +clsqf02(l,k,m)*hc2d(8,l,k,m)                   &
                        +pc2d(0,l,k,m)

        qf2d(0,l,k,m)   =phif(m,l,k)
                  
        qf2d(1,l,k,m)   =clsqf11(l,k,m)*hc2d(3,l,k,m)                   &
                        +clsqf12(l,k,m)*hc2d(7,l,k,m)                   &
                        +pc2d(1,l,k,m)

        qf2d(5,l,k,m)   =clsqf11(l,k,m)*hc2d(1,l,k,m)                   &
                        +clsqf12(l,k,m)*hc2d(5,l,k,m)                   &
                        +pc2d(5,l,k,m)

        c22hc2d8        =clsqf22(l,k,m)*hc2d(8,l,k,m)
        
        qf2d(2,l,k,m)   =clsqf21(l,k,m)*hc2d(4,l,k,m)                   &
                        +c22hc2d8+pc2d(2,l,k,m)
        
        qf2d(6,l,k,m)   =clsqf21(l,k,m)*hc2d(2,l,k,m)                   &
                        +c22hc2d8+pc2d(6,l,k,m)
        
        qf2d(3,l,k,m)   =clsqf31(l,k,m)*hc2d(3,l,k,m)                   &
                        +clsqf32(l,k,m)*hc2d(7,l,k,m)                   &
                        +pc2d(3,l,k,m)
                      
        qf2d(7,l,k,m)   =clsqf31(l,k,m)*hc2d(1,l,k,m)                   &
                        +clsqf32(l,k,m)*hc2d(5,l,k,m)                   &
                        +pc2d(7,l,k,m)

        c42hc2d8        =clsqf42(l,k,m)*hc2d(8,l,k,m)
        
        qf2d(4,l,k,m)   =clsqf41(l,k,m)*hc2d(4,l,k,m)                   &
                        +c42hc2d8+pc2d(4,l,k,m)

        qf2d(8,l,k,m)   =clsqf41(l,k,m)*hc2d(2,l,k,m)                   &
                        +c42hc2d8+pc2d(8,l,k,m)
!       cross terms
        qf2d(9,l,k,m)=clsqfx1y1(l,k,m)*hc2d(6,l,k,m)+pc2d(9,l,k,m)
        qf2d(10,l,k,m)=clsqf1221(l,k,m)*hc2d(5,l,k,m)+pc2d(10,l,k,m)
        qf2d(11,l,k,m)=clsqf1221(l,k,m)*hc2d(7,l,k,m)+pc2d(11,l,k,m)
        qf2d(12,l,k,m)=clsqfx2y2(l,k,m)*hc2d(8,l,k,m)+pc2d(12,l,k,m)
        qf2d(13,l,k,m)=clsqf1331(l,k,m)*hc2d(6,l,k,m)+pc2d(13,l,k,m)
        qf2d(14,l,k,m)=clsqf1331(l,k,m)*hc2d(6,l,k,m)+pc2d(14,l,k,m)
    enddo
    enddo
    return
end
