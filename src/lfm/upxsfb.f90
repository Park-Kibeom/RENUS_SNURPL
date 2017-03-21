subroutine upxsfb
      use param
      use MASTERXSL, ONLY : LFM_coeff, NUFS,FISS,CAPT,TRAN,SCA,ABSO,KAPA,NTYPE1,NTYPE2
      include 'global.h'
      include 'ffdm.h'
      include 'times.h'
      include 'xsec.h'
!      include 'fbxs.h'
      include 'geom.h'
      include 'files.h'
      include 'xesm.h'  ! 2013_09_27 . scb
      logical, save :: first=.TRUE.
!      
      integer :: i,j
      integer :: negative
      real :: lf, lt, sumv
      real :: zeros=1.0
      logical :: off
      
      allocate(lkg(nxy,ng))
      allocate(philkg(nxy,ng))
      off=0
      lkg=0;       philkg=0
      do k=1,nz
        ka=ktoka(k)
        do l=1,nxy
          la=ltola(l)
          iat=iassytyp(la)
          ia=latoia(la)
          ja=latoja(la)
          icomp=icompnumz(ka,iat)
          do m=1,ng
            philkg(la,m)=philkg(la,m)+phif(m,l,k)
            do i=1,2 !left/right
              do j=1,2 !x,y   consider radial leakage
!              do j=1,3 !x,y,z
                if(i.eq.1)then
                !lkg(la,m) =lkg(la,m) -jnet(i,m,l,k,j)*rhmesh(j,l,k)
                lkg(l,m) =lkg(l,m) -jnet(i,m,l,k,j)*rhmesh(j,l,k)
                else
                !lkg(la,m) =lkg(la,m) +jnet(i,m,l,k,j)*rhmesh(j,l,k)
                lkg(l,m) =lkg(l,m) +jnet(i,m,l,k,j)*rhmesh(j,l,k)
                endif
              enddo
            enddo
            lkg(l,m)=lkg(l,m)/phif(m,l,k)
          enddo
        enddo
      enddo
      do la=1, nxya
        do m=1, ng
          !lkg(la,m)=lkg(la,m)/philkg(la,m)
        enddo
      enddo
      
      negative=0;
      do k=1,nz
        ka=ktoka(k)
        do l=1,nxy
          la=ltola(l)
          iat=iassytyp(la)
          ia=latoia(la)
          ja=latoja(la)
          icomp=icompnumz(ka,iat)
          sumv=0
          off=0
          !if( esigd(1,icomp).eq.0 )then
          !    cycle
          !endif          
          do m=1,ng
            
            lf=lkg(l,1)  !*0.0
            lt=lkg(l,2) !*0.0
            !---D
            if(off)then
                !xsdf(m,l,k)=xsdf(m,l,k)
            else
                xsdf(m,l,k)=xsdf(m,l,k)*(1+ LFM_coeff(icomp,m,TRAN,1)*lf + LFM_coeff(icomp,m,TRAN,2)*lt)
            endif
            xstrf(m,l,k)=1.0/3/xsdf(m,l,k)
            xsdf2(m,l,k)=9.d0/7.d0*xsdf(m,l,k)   ! 2013_07_15 . scb
            if( xsdf(m,l,k) .lt. 0) negative=negative+1
            
            !---sig_A
            if(off)then
                !xsaf(m,l,k)=xsaf(m,l,k)
            else
                xsaf(m,l,k)=xsaf(m,l,k)*(1+ LFM_coeff(icomp,m,ABSO,1)*lf + LFM_coeff(icomp,m,ABSO,2)*lt)
            endif
            
            if( xsaf(m,l,k) .lt. 0) negative=negative+1
            !---sig_nf            
            if(off)then
                !xsnff(m,l,k)=xsnff(m,l,k)
            else
                xsnff(m,l,k)=xsnff(m,l,k)*(1+ LFM_coeff(icomp,m,NUFS,1)*lf + LFM_coeff(icomp,m,NUFS,2)*lt)
                xsff(m,l,k) =xsff(m,l,k) *(1+ LFM_coeff(icomp,m,KAPA,1)*lf + LFM_coeff(icomp,m,KAPA,2)*lt)
                xskpf(m,l,k)=xskpf(m,l,k)*(1+ LFM_coeff(icomp,m,KAPA,1)*lf + LFM_coeff(icomp,m,KAPA,2)*lt)
            endif
            
            if( xsnff(m,l,k) .lt. 0) negative=negative+1
            !
            !xsff(m,l,k)=xsnff(m,l,k)/2.3  ! 2013_10_02 . scb      
            
            !if(off)then
            !xskpf(m,l,k)=esigkf(m,icomp)
            !else
            !xskpf(m,l,k)=asigkf(m,icomp)*ratio + csigkf(m,icomp)*fratio         + esigkf(m,icomp)
            !endif
            
            if( xskpf(m,l,k) .lt. 0) negative=negative+1
            
            !---sig_s
            if(off)then
                !if(m.eq.1)then
                !  xssf(m,m+1,l,k) =esigs(m,icomp)
                !  if( xssf(m,m+1,l,k) .lt. 0) negative=negative+1
                !else
                !  xssf(m,m-1,l,k) =esigs(m,icomp)         
                !  if( xssf(m,m-1,l,k) .lt. 0) negative=negative+1       
                !endif
            else
                if(m.eq.1)then
                  !xssf(m,m+1,l,k) = asigs(m,icomp)*ratio + csigs(m,icomp)*fratio + esigs(m,icomp)
                  xssf(m,m+1,l,k)=xssf(m,m+1,l,k)*(1+ LFM_coeff(icomp,m,SCA,1)*lf + LFM_coeff(icomp,m,SCA,2)*lt)
                  if( xssf(m,m+1,l,k) .lt. 0) negative=negative+1
                else
                  !xssf(m,m-1,l,k) = asigs(m,icomp)*ratio + csigs(m,icomp)*fratio + esigs(m,icomp)         
                  xssf(m,m-1,l,k)=xssf(m,m-1,l,k)*(1+ LFM_coeff(icomp,m,SCA,1)*lf + LFM_coeff(icomp,m,SCA,2)*lt)
                  if( xssf(m,m-1,l,k) .lt. 0) negative=negative+1       
                endif
            endif
            
            
            xbeta(m,l,k,XDIR)=xsdf(m,l,k)*rhmesh(XDIR,l,k)
            xbeta(m,l,k,YDIR)=xsdf(m,l,k)*rhmesh(YDIR,l,k)
	        xbeta(m,l,k,ZDIR)=xsdf(m,l,k)*rhmesh(ZDIR,l,k)
	        xbeta2(m,l,k,XDIR)=xsdf2(m,l,k)*rhmesh(XDIR,l,k)
            xbeta2(m,l,k,YDIR)=xsdf2(m,l,k)*rhmesh(YDIR,l,k)
	        xbeta2(m,l,k,ZDIR)=xsdf2(m,l,k)*rhmesh(ZDIR,l,k)
          
            
                if(m.eq.1)then
          write(10,'(2i5,100e17.5)') l,m, xsdf(m,l,k),  xsaf(m,l,k),  xsnff(m,l,k), xsff(m,l,k), xskpf(m,l,k),  xssf(m,m+1,l,k)
                else  
          write(10,'(2i5,100e17.5)') l,m, xsdf(m,l,k),  xsaf(m,l,k),  xsnff(m,l,k), xsff(m,l,k), xskpf(m,l,k),  xssf(m,m-1,l,k)
                endif
            
          enddo !group
        enddo !nxy
      enddo !nz
      if( negative .ne. 0 )then
          write(10,*) '!!!'
      endif
      
      
      do k=1,nz
        do l=1,nxy
          do m=1,ng
            summ=0
            do md=1,ng
              if(m .ge. xssfs(md,l,k) .and. m.le.xssfe(md,l,k)) then
                summ=summ+xssf(m,md,l,k)
              endif
            enddo
            if(usesiga) then
              xstf(m,l,k)=xsaf(m,l,k)+summ
            else
              xsaf(m,l,k)=xstf(m,l,k)-summ
            endif
          enddo
        enddo
      enddo
      
1000  continue
      do l=1,nxy
          write(10,*) lkg(l,:)
      enddo
     write(10,*) '---'
      
          
      

      return
end subroutine

