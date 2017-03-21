    subroutine caltftr(k,lchan,tcool1,htcoef1,qf,dr2odt,tw2odt,tfupd)
! steady-state fuel temperature calculation at the specified channel (k,lchan)      
      use param
      use allocs

      include 'global.h'
      include 'thgeom.inc'
      include 'thfuel.inc'
      include 'thcntl.inc'
      include 'thfdbk.inc'
      include 'geom.h'   ! 2015_08_03 . scb
      
      logical, intent(in) :: tfupd
      
      logical, save :: first=TRUE
      real :: kconv,kconv1,kconv2,kconv4,kconv4b,kfi,kmr,kml,kmrb,kmlb
      real :: krat   ! 2015_08_03 . scb
      real, pointer, save, dimension(:) :: &
                              kf   & !thermal conductivity of uo2 at mesh pointss
                             ,kfb  & !thermal conductivity of uo2 at mesh pointss
                             ,kfm  & ! kf at the middle points
                             ,kfmp & !thermal conductivity of uo2 at mesh pointss
                             ,x    & ! fuel temp
                             ,rhocp ,ad ,al ,au ,b
      
      if(first) then
        call dmalloc(kf,nrp5)
        call dmalloc(kfb,nrp5)
        call dmalloc(kfm,nrp5)
        call dmalloc(kfmp,nrp5)
        call dmalloc(x,nrp5)
        call dmalloc(rhocp,nrp5)
        call dmalloc(ad,nrp5)
        call dmalloc(al,nrp5)
        call dmalloc(au,nrp5)
        call dmalloc(b,nrp5)
        first=FALSE
      endif      
!
! 2015_08_03 . scb
      l=lchantol(lchan)
      la=ltola(l)
      iat=iassytyp(la)
      krat=kratio(iat)
! added end      

! update coefficients
      tdopld=tdopl(k,lchan)
      do i=1,nrp4
         x(i)=tfuel(i,k,lchan)+CKELVIN
      enddo
      do i=1,nrp1
         kfi=fkf(x(i))
         kfi=krat*kfi  ! 2015_08_03 . scb
         kf(i)=thetaf*kfi
         kfb(i)=thetafb*kfi
         rhocp(i)=frhocpf(x(i))*dr2odt
      enddo
! coefficient at the middle points
      do i=1,nr
         kfm(i)=0.5*(kf(i)+kf(i+1))
         kfmp(i)=0.5*(kfb(i)+kfb(i+1))
      enddo
      do i=nrp2,nrp4
         kfi=fkc(x(i))
         kf(i)=thetaf*kfi
         kfb(i)=thetafb*kfi
         rhocp(i)=frhocpc(x(i))*tw2odt
      enddo
      m=nrp3
      kmr=0.5*(kf(m)+kf(m+1))
      kml=0.5*(kf(m)+kf(m-1))
      kmrb=0.5*(kfb(m)+kfb(m+1))
      kmlb=0.5*(kfb(m)+kfb(m-1))
! setup linear system
      do i=1,nrp4
         x(i)=tfuel(i,k,lchan)
      enddo
      qfd2=qf*delr2
      m=1
      ad(m)=rhocp(m)+4*kf(m)
      au(m)=-4*kf(m)
      b(m)=qfd2+x(m)*rhocp(m)-4*kfb(m)*x(m)+4*kfb(m)*x(m+1)
      i=m
      do m=2,nr
         ri=1/float(i)
         ad(m)=rhocp(m)+(kfm(i)+kfm(m)+0.5*(kfm(m)-kfm(i))*ri)
         al(m)=-kfm(i)*(1-0.5*ri)
         au(m)=-kfm(m)*(1+0.5*ri)
         b(m)=qfd2+x(m)*rhocp(m)-(kfmp(i)+kfmp(m)+0.5*(kfmp(m)-kfmp(i))*ri)*x(m)  &
             +kfmp(i)*(1-0.5*ri)*x(m-1)+kfmp(m)*(1+0.5*ri)*x(m+1)
         i=m
      enddo
      m=nrp1
      alpha=kgap*(1-kf(m-1)/kf(m))
      alphab=alpha/thetaf*thetafb
      ad(m)=rhocp(m)+2*(kf(m)+kgap*(1+0.5/nr))+alpha
      al(m)=-2*kf(m)
      au(m)=-2*kgap*(1+0.5/nr)-alpha
      b(m)=qfd2+x(m)*rhocp(m)-(2*(kfb(m)+kgapb*(1+0.5/nr))+alphab)*x(m)  &
          +2*kfb(m)*x(m-1)+(2*kgapb*(1+0.5/nr)+alphab)*x(m+1)
      m=nrp2
      alpha=2*kgap2*(kf(m+1)/kf(m)-1)
      alphab=alpha/thetaf*thetafb
      ad(m)=rhocp(m)+8*kf(m)+kgap4-alpha
      al(m)=-kgap4+alpha
      au(m)=-8*kf(m)
      b(m)=x(m)*rhocp(m)-(8*kfb(m)+kgap4b)*x(m)+kgap4b*x(m-1)  &
          +8*kfb(m)*x(m+1)+x(m)*alphab-x(m-1)*alphab
      m=nrp3
      ad(m)=rhocp(m)+4*(kmr+kml)+tworm*(kmr-kml)
      al(m)=-kml*(4-tworm)
      au(m)=-kmr*(4+tworm)
      b(m)=x(m)*rhocp(m)-(4*(kmrb+kmlb)+tworm*(kmrb-kmlb))*x(m)  &
          +kmlb*(4-tworm)*x(m-1)+kmrb*(4+tworm)*x(m+1)
      m=nrp4
      kconv1=thetaf*htcoef1*tw
      alpha=2*kconv1*(1-kf(m-1)/kf(m))
      alphab=alpha/thetaf*thetafb
      kconv=htcoef1*tw*(4+tw/rw)
      kconv4=thetaf*kconv
      kconv4b=thetafb*kconv
      ad(m)=rhocp(m)+8*kf(m)+kconv4+alpha
      al(m)=-8*kf(m)
      b(m)=x(m)*rhocp(m)-(8*kfb(m)+kconv4b+alphab)*x(m)+8*kfb(m)*x(m-1)  &
          +(kconv+alpha/thetaf)*tcool1


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
      enddo
      if(tfupd) then
         do i=1,nrp4
            tfuel(i,k,lchan)=x(i)
         enddo
      endif
      
      tdoplold=tdopl(k,lchan)
      tfuel(nrp5,k,lchan)=ftfavg(x,nr,delr2)
      if(ifeffdopl) then !txk undocumented option for doppler temperature?
         tdopl(k,lchan)=sqrt(wfcl*x(1)+wfsurf*x(nrp1)+CKELVIN)
      else
         tdopl(k,lchan)=sqrt(tfuel(nrp5,k,lchan)+CKELVIN)
      endif
      tdoplmax=max(tdoplmax,abs(1-tdoplold/tdopl(k,lchan)))
!
      return
      
    end subroutine
