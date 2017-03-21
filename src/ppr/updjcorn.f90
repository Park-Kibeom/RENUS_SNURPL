subroutine updjcorn(k)
    use ppr
    use geom,   only : ng,nxy,hmesh
    use xsec,   only : xsdf

    integer                 :: k
    real                    :: alpha
    
    do m=1,ng
    do l=1,nxy
        alpha=-2.*xsdf(m,l,k)/hmesh(XDIR,l,k)

        do ic=1,4
            jcornx(ic,l,k,m)=alpha*(																	&
                               cpjxh1(ic,l,k,m)*hc2d(1,l,k,m)         &
                              +cpjxh2(ic,l,k,m)*hc2d(2,l,k,m)         &
                              +cpjxh5(ic,l,k,m)*hc2d(5,l,k,m)         &
                              +cpjxh6(ic,l,k,m)*hc2d(6,l,k,m)         &
                              +cpjxh7(ic,l,k,m)*hc2d(7,l,k,m)         &
                              +cpjxh8(ic,l,k,m)*hc2d(8,l,k,m)         &
                              +pc2d(5,l,k,m)                          &
                              +cpjxp6(ic,l,k,m)*pc2d(6,l,k,m)         &
                              +cpjxp7(ic,l,k,m)*pc2d(7,l,k,m)         &
                              +cpjxp8(ic,l,k,m)*pc2d(8,l,k,m)         &  
                              +cpjxp9(ic,l,k,m)*pc2d(9,l,k,m)         &
                              +pc2d(10,l,k,m)                         &
                              +cpjxp11(ic,l,k,m)*pc2d(11,l,k,m)       &
                              +cpjxp12(ic,l,k,m)*pc2d(12,l,k,m)       &
                              +cpjxp13(ic,l,k,m)*pc2d(13,l,k,m)       &
                              +cpjxp14(ic,l,k,m)*pc2d(14,l,k,m))

            jcorny(ic,l,k,m)=alpha*(																	&
                               cpjyh3(ic,l,k,m)*hc2d(3,l,k,m)         &
                              +cpjyh4(ic,l,k,m)*hc2d(4,l,k,m)         &
                              +cpjyh5(ic,l,k,m)*hc2d(5,l,k,m)         &
                              +cpjyh6(ic,l,k,m)*hc2d(6,l,k,m)         &
                              +cpjyh7(ic,l,k,m)*hc2d(7,l,k,m)         &
                              +cpjyh8(ic,l,k,m)*hc2d(8,l,k,m)         &
                              +pc2d(1,l,k,m)                          &
                              +cpjyp2(ic,l,k,m)*pc2d(2,l,k,m)         &
                              +cpjyp3(ic,l,k,m)*pc2d(3,l,k,m)         &
                              +cpjyp4(ic,l,k,m)*pc2d(4,l,k,m)         &  
                              +cpjyp9(ic,l,k,m)*pc2d(9,l,k,m)         &
                              +cpjyp10(ic,l,k,m)*pc2d(10,l,k,m)       &
                              +pc2d(11,l,k,m)                         &
                              +cpjyp12(ic,l,k,m)*pc2d(12,l,k,m)       &
                              +cpjyp13(ic,l,k,m)*pc2d(13,l,k,m)       &
                              +cpjyp14(ic,l,k,m)*pc2d(14,l,k,m))
        enddo
    enddo
    enddo    
        !
    return
end
