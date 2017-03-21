subroutine caltrlcusp1n(iftran,l,krod,rodfrac,hfine,ntfine,trlfine)
    use const
    use geom,       only    :  ng,hmesh,nxs,nxe,nodel,volnode
    use tran,       only    :  srctrf,rveloftm,betap,pinvf
    use tranxsec,   only    :  xschifd
    use sfam,       only    :  phif,psi,jnet
    implicit none
    
    integer                 :: l,krod,ntfine
    logical                 :: iftran
    real                    :: hfine(ntfine)
    real                    :: rodfrac
    real                    :: trlfine(ng,ntfine)
    
    integer                 :: k5,k,kc,kl,kr,kfine,kreg,kregfine,idir,j,i,m
    real                    :: avgjnet(2)
    real                    :: trlcff(ng,0:2,-1:1) ! ng, order, (lower, center, upper)
    real                    :: avgtrl3(ng, 0:2), hmesh3(0:2) ! 0-center, 1-left, 2-right
    real                    :: spsrctr,rvol   
    real                    :: rhx,rhy,rhz,zfinel,zfiner,hdelta

    spsrctr=0
    avgjnet=0
    trlcff=0
! obtain transverse leakage
    do k5=-1,1
        k=krod+k5
        do m=1,ng
            avgjnet(XDIR)=(jnet(RIGHT,m,l,k,XDIR)-jnet(LEFT,m,l,k,XDIR))/hmesh(XDIR,l,k)
            avgjnet(YDIR)=(jnet(RIGHT,m,l,k,YDIR)-jnet(LEFT,m,l,k,YDIR))/hmesh(YDIR,l,k)

            if(iftran) then
                rvol=1/volnode(l,k)
                spsrctr=srctrf(m,l, k)*rvol-(rveloftm(m,l,k)+pinvf(m,l,k))*phif(m,l,k)      &
                                -xschifd(m,l,k)*(1-betap(l,k))*psi(l,k)*rvol
            endif            

            trlcff(m,0,k5)=avgjnet(XDIR) + avgjnet(YDIR) - spsrctr         
        enddo
    enddo   !k  
    
    kc=0;kl=-1;kr=+1
    hmesh3(CENTER)=hmesh(ZDIR,l,krod+kc)
    hmesh3(LEFT)=hmesh(ZDIR,l,krod+kl)
    hmesh3(RIGHT)=hmesh(ZDIR,l,krod+kr)
    do m=1,ng
        avgtrl3=0;

        avgtrl3(:,CENTER)=trlcff(:,0,kc)
        avgtrl3(:,LEFT)=trlcff(:,0,kl)
        avgtrl3(:,RIGHT)=trlcff(:,0,kr)

        call trlcffbyintg(avgtrl3, hmesh3, ng,               &
                trlcff(:,1,kc),trlcff(:,2,kc))
    enddo    
    
    kfine = 0
    zfinel = -1
    rhz = 1/hmesh(ZDIR,l,krod)
    k=0
    do kfine=1,ntfine
        hdelta= 2*hfine(kfine)*rhz
        zfiner = zfinel + hdelta
        do m=1,ng
            trlfine(m, kfine) = ((trlcff(m,0,k)-trlcff(m,2,k)*0.5)*(zfiner-zfinel)  &
                              + trlcff(m,1,k)*0.5*(zfiner**2-zfinel**2)             &
                              + trlcff(m,2,k)*0.5*(zfiner**3-zfinel**3))/hdelta
        enddo
        zfinel=zfiner
    enddo
end subroutine