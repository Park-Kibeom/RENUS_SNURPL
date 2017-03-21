subroutine caltrlcusp3n(iftran,l,krod,rodfrac)
    use const
    use geom,       only    :  ng,hmesh,nxe,nodel,volnode  ! nxs - doesnt used, 06.01.16_alc
    use tran,       only    :  srctrf,rveloftm,betap,pinvf
    use tranxsec,   only    :  xschifd
    use sfam,       only    :  phif,psi,jnet
    use decusping3n
    implicit none
    
    integer                 :: l,krod
    logical                 :: iftran
    real                    :: rodfrac
    
    integer                 :: k5,k,kc,kl,kr,kfine,kreg,kregfine,idir,j,i,m
    real                    :: avgjnet(2)
    real                    :: trlcff(ng,0:2,-2:2) ! ng, order, (lower, center, upper)
    real                    :: avgtrl3(ng, 0:2), hmesh3(0:2) ! 0-center, 1-left, 2-right
    real                    :: spsrctr,rvol   
    real                    :: rhx,rhy,rhz,zfinel,zfiner,hdelta

    spsrctr=0
    avgjnet=0
    trlcff=0
! obtain transverse leakage
    do k5=-2,2
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
    
    do k=-1,1
        kc=k;kl=k-1;kr=k+1
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
    enddo  
    
    kfine = 0
    do kreg=1,4
        k = kreg-2
        if(k.ge.1) k=k-1
        zfinel = -1
        if(kreg.eq.3) zfinel=1-2*rodfrac
        rhz = 1/hmesh(ZDIR,l,krod+k)
        do kregfine=1,nfine(kreg)
            kfine = kfine + 1 
            hdelta= 2*hfine(kreg)*rhz
            zfiner = zfinel + hdelta
            do m=1,ng
                trlfine(m, kfine) = ((trlcff(m,0,k)-trlcff(m,2,k)*0.5)*(zfiner-zfinel)  &
                                  + trlcff(m,1,k)*0.5*(zfiner**2-zfinel**2)             &
                                  + trlcff(m,2,k)*0.5*(zfiner**3-zfinel**3))/hdelta
                                  
            enddo
            zfinel=zfiner
        enddo
    enddo
      
end subroutine