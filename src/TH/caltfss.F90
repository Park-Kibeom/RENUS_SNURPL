    subroutine caltfss(k,lchan,tcool1,htcoef1,qf)
! steady-state fuel temperature calculation at the specified channel (k,lchan)      
      use param
      use allocs

      include 'global.h'
      include 'thgeom.inc'
      include 'thfuel.inc'
      include 'thcntl.inc'
      include 'thfdbk.inc'
      include 'geom.h'   ! 2015_08_03 . scb
      
      logical, save :: first=TRUE
      real :: kconv,kconv1,kconv2,kconv4,kconv4b,kfi,kmr,kml,kmrb,kmlb
      real :: krat   ! 2015_08_03 . scb
      real, pointer, save, dimension(:) ::  &
                              kf  & !thermal conductivity of uo2 at mesh pointss
                             ,kfm & ! kf at the middle points
                             ,x   & ! fuel temp
                             ,xd,ad,al,au,b
      
      if(first) then
        call dmalloc(kf,nrp5)
        call dmalloc(kfm,nrp5)
        call dmalloc(x,nrp5)
        call dmalloc(xd,nrp5)
        call dmalloc(ad,nrp5)
        call dmalloc(al,nrp5)
        call dmalloc(au,nrp5)
        call dmalloc(b,nrp5)
        first=FALSE
      endif
      
      errtf=CINF
      itr=0
      
! 2015_08_03 . scb
      l=lchantol(lchan)
      la=ltola(l)
      iat=iassytyp(la)
      krat=kratio(iat)
! added end      
      do while(errtf.gt.epstf)
         itr=itr+1
! update coefficients
         do i=1,nrp4
            x(i)=tfuel(i,k,lchan)+CKELVIN
         enddo
         do i=1,nrp1
            kf(i)=fkf(x(i))
            kf(i)=krat*kf(i)  ! 2015_08_03 . scb
         enddo
! coefficient at the middle points
         do i=1,nr
            kfm(i)=0.5*(kf(i)+kf(i+1))
         enddo
         do i=nrp2,nrp4
            kf(i)=fkc(x(i))
         enddo
         m=nrp3
         kmr=0.5*(kf(m)+kf(m+1))
         kml=0.5*(kf(m)+kf(m-1))      
! setup linear system
         do i=1,nrp4
            x(i)=tfuel(i,k,lchan)
            xd(i)=x(i)
         enddo
         qfd2=qf*delr2
         m=1
         ad(m)=4*kf(m)
         au(m)=-4*kf(m)
         b(m)=qfd2
         i=m
         do m=2,nr
            ri=1/float(i)
            ad(m)=kfm(i)+kfm(m)+0.5*(kfm(m)-kfm(i))*ri
            al(m)=-kfm(i)*(1-0.5*ri)
            au(m)=-kfm(m)*(1+0.5*ri)
            b(m)=qfd2
            i=m
         enddo
         m=nrp1
         alpha=kgap*(1-kf(m-1)/kf(m))
         ad(m)=2*(kf(m)+kgap*(1+0.5/nr))+alpha
         al(m)=-2*kf(m)
         au(m)=-2*kgap*(1+0.5/nr)-alpha
         b(m)=qfd2
         m=nrp2
         alpha=2*kgap2*(kf(m+1)/kf(m)-1)
         ad(m)=8*kf(m)+kgap4-alpha
         al(m)=-kgap4+alpha
         au(m)=-8*kf(m)
         b(m)=0
         m=nrp3
         ad(m)=4*(kmr+kml)+tworm*(kmr-kml)
         al(m)=-kml*(4-tworm)
         au(m)=-kmr*(4+tworm)
         b(m)=0
         m=nrp4
         kconv1=htcoef1*tw
         alpha=2*kconv1*(1-kf(m-1)/kf(m))
         kconv=htcoef1*tw*(4+tw/rw)
         ad(m)=8*kf(m)+kconv+alpha
         al(m)=-8*kf(m)
         b(m)=(kconv+alpha)*tcool1
         
! solve the tridiagonal system by gauss elimination
         im1=1
         do i=2,nrp4
            aldi=al(i)/ad(im1)
            ad(i)=ad(i)-aldi*au(im1)
            b(i)=b(i)-aldi*b(im1)
            im1=i
         enddo
         i=nrp4
         ip1=nrp4
         x(i)=b(i)/ad(i)
         do i=nrp3,1,-1
            x(i)=(b(i)-au(i)*x(ip1))/ad(i)
            ip1=i
            if(x(i).lt.0)  x(i)=1    ! 2015_07_10 . scb fix
         enddo
         errtf=0
         do i=1,nrp4
            tfuel(i,k,lchan)=x(i)
            errtf=max(errtf,abs(x(i)-xd(i)))
         enddo                 
      enddo

!#ifdef DBG
!      if(lchan.eq.1)                                                        &
!     &     print '(3i5,3f10.2,f10.4)',k,lchan,itr,x(1),x(nrp4),tcool1,relp(k,lchan)
!#endif
      tfuel(nrp5,k,lchan)=ftfavg(x,nr,delr2)
      tdoplold=tdopl(k,lchan)
      if(ifeffdopl) then !txk undocumented option for doppler temperature?
         tdopl(k,lchan)=sqrt(wfcl*x(1)+wfsurf*x(nrp1)+CKELVIN)                        ! NEA w=0.64
!         tdopl(k,lchan)=sqrt(wfcl*tfuel(nrp5,k,lchan)+wfsurf*x(nrp1)+CKELVIN)        !  Studsvik w=0.78  
      else
         tdopl(k,lchan)=sqrt(tfuel(nrp5,k,lchan)+CKELVIN)
      endif
      tdoplmax=max(tdoplmax,abs(1-tdoplold/tdopl(k,lchan)))      
      
      return
      
    end subroutine