! 2014_10_06 . scb
  subroutine updsrcdtcff_sp3senm(idir,l,k,m,i)    
    use const
    use bdf,        only : bdfcoef,bdfordersave,dtcff2
    use tranxsec,   only : rvelof
    use xsec,       only : xstrf,xsdf,xstf,xsdf2,xstfsp3
    use geom,       only : hmesh
    use sp3senm,    only : phishp
    use senm1d,     only : tdmd
    
    integer                 :: idir,l,k,m,im
    
    real                    :: srcdtcff(0:2,2), coef(2), rh
    real                    :: dtcff0(0:2,1:2)
    integer :: nbdf
    
    nbdf=bdfordersave
    rh=1.d0/hmesh(idir,l,k)
    
    do im=1,2
      dtcff0(0:2,im) = dtcff2(:,im,m,l,k,idir,0)    
        
      dtcff0(:,im)=0.d0
      do icff=0,2
        do iorder=1,nbdf
          dtcff0(icff,im) = dtcff0(icff,im) + bdfcoef(iorder)*dtcff2(icff,im,m,l,k,idir,iorder)
        enddo
        dtcff0(icff,im) = rvelof(m,l,k)*dtcff0(icff,im)
      enddo
    enddo
        
    dtcff0(0,1) = dtcff0(0,1) + 2.d0/3.d0*rh*(3.d0*phishp(1,2,XDIR,l,k,m)+10.d0*phishp(1,4,XDIR,l,k,m) &
                                                +6.d0*phishp(2,2,XDIR,l,k,m)+20.d0*phishp(2,4,XDIR,l,k,m))
    dtcff0(1,1) = dtcff0(1,1) + 10.d0*rh*(phishp(1,3,XDIR,l,k,m)+2*phishp(2,3,XDIR,l,k,m))
    dtcff0(2,1) = dtcff0(2,1) + 70.d0/3.d0*rh*(phishp(1,4,XDIR,l,k,m)+2*phishp(2,4,XDIR,l,k,m)) 

    
    dtcff0(0,2) = dtcff0(0,2) + 6.d0/7.d0*rh*(3.d0*phishp(2,2,XDIR,l,k,m)+10.d0*phishp(2,4,XDIR,l,k,m))
    dtcff0(1,2) = dtcff0(1,2) + 90.d0/7.d0*rh*phishp(2,3,XDIR,l,k,m)
    dtcff0(2,2) = dtcff0(2,2) + 30.d0*rh*phishp(2,4,XDIR,l,k,m)

    coef(1)=-1.d0/(xstrf(m,l,k) + bdfcoef(0)*rvelof(m,l,k))
    coef(2)=-1.d0/(xstfsp3(m,l,k) + bdfcoef(0)*rvelof(m,l,k))
    do im=1,2
      do icff=0,2
        dtcff0(icff,im) = coef(im)*dtcff0(icff,im)
      enddo
      dtcff2(:,im,m,l,k,idir,0) = dtcff0(0:2,im)  
    enddo
          
    coef(1)=2.d0*rvelof(m,l,k)*rh*xsdf(m,l,k)
    coef(2)=2.d0*rvelof(m,l,k)*rh*xsdf2(m,l,k)    
    do im=1,2   
      srcdtcff(:,im)=0.d0
      do icff=0,2
        do iorder=0,nbdf
          srcdtcff(icff,im) = srcdtcff(icff,im) +  bdfcoef(iorder)*dtcff2(icff,im,m,l,k,idir,iorder)
        enddo
        srcdtcff(icff,im) = coef(im)*srcdtcff(icff,im)
      enddo    
    enddo
    
    do im=1,2
      do icff=0,2
        tdmd(im,icff,i,m)=srcdtcff(icff,im)
      enddo
    enddo    

  end subroutine
  
