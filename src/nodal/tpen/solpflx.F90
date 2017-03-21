! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
    subroutine solpflx

! point flux solver using CPB relation

	    use const
	    use geomhex
	    use tpen
	    use xsec
	    use sfam, only : phi, phif
	    implicit none

      real,allocatable,save,dimension(:) :: psrc,pcoef,psrctem, &
                                            xsdmwt(:,:)
	    real :: bsflx(6), radfl, xsdmwtb, wts1, wts2, wts3, wts4, &
	            adfl, reflratl, dat1, dat2, dat3, pbdvl, dathflx
	    integer :: inpt(6), iz, ig, ip, ih, it, nn, iter
	    logical,save :: first=TRUE
      real :: phipt(6)   ! 2014_05_15 . scb adf hex

      if(ng.ne.2) then
         call solpflxmg
         return
      endif
      
      if(first) then
        allocate(psrc(ncorn),pcoef(ncorn),psrctem(ncorn))
        allocate(xsdmwt(3,nassy))     
        psrc=0.d0
        pcoef=0.d0
        psrctem=0.d0
        xsdmwt=0.d0
        first=.FALSE.
      endif
!
      do iz=1,nz
        do ig=1,ng2
          do ip=1,nxpnt
            psrc(ip)=0.d0
            pflxt(ip)=pflx(ig,ip,iz)
            pcoef(ip)=0.d0
          enddo
          do ih=1,nassy
            radfl=1.
            xsdmwtb=wtass(ih)*xsd(ig,ih,iz)
            xsdmwt(1,ih)=chlval(1)*xsdmwtb*radfl
            xsdmwt(2,ih)=chlval(2)*xsdmwtb*radfl
            xsdmwt(3,ih)=chlval(3)*xsdmwtb*radfl
            wts1=xsdmwtb*chlval(4)
            wts2=xsdmwtb*chlval(5)
            wts3=xsdmwtb*chlval(6)
            wts4=xsdmwtb*chlval(7)
            do it=1,6
              nn=neignd(it,ih)
              if(nn.eq.0) then
!                adfl=xsadf(ig,ih,iz)
                adfl=1.
                reflratl=(adfl-2*alxr)/(adfl+2*alxr)
                bsflx(it)=2.*(1.+reflratl)*cnto(ig,it,ih,iz)
              else
                bsflx(it)=2.*(cnto(ig,it,ih,iz)  &
                         +cnto(ig,neigjin(it,ih),nn,iz))
              endif
              bsflx(it)=bsflx(it)/xsadf(it,ig,ih,iz)    ! 2014_05_15 . scb adf hex
            enddo
#ifdef DEBUG
	      if(ig.eq.1) print '(a10,10F9.4)', 'BSFlux', bsflx(1)
#endif
!            dathflx=wts4*hflx(ig,ih,iz)*xsadf(ig,ih,iz)
            dathflx=wts4*hflx(ig,ih,iz)
            do it=1,3
              nn=neigpt(it,ih)
              dat1=bsflx(it)+bsflx(mp5(it))
              dat2=bsflx(mp2(it))+bsflx(mp3(it))
              dat3=wts3*(bsflx(mp1(it))+bsflx(mp4(it)))
              pbdvl=pbdv(ig,nn)*codpnt(nn)*wtass(ih)
              psrc(nn)=psrc(nn)+(wts1*dat1+wts2*dat2+dat3+dathflx)*radfl
              !pcoef(nn)=pcoef(nn)+(xsdmwtb+pbdvl)*radfl
              pcoef(nn)=pcoef(nn)+xsdmwtb/xscdf(it,ig,ih,iz)+pbdvl   ! 2014_05_15 . scb adf hex
              nn=neigpt(mp3(it),ih)
              pbdvl=pbdv(ig,nn)*codpnt(nn)*wtass(ih)
              psrc(nn)=psrc(nn)+(wts1*dat2+wts2*dat1+dat3+dathflx)*radfl
              !pcoef(nn)=pcoef(nn)+(xsdmwtb+pbdvl)*radfl
              pcoef(nn)=pcoef(nn)+xsdmwtb/xscdf(mp3(it),ig,ih,iz)+pbdvl    ! 2014_05_15 . scb adf hex
            enddo
          enddo

          do iter=1, nitrpfc
            do ip=1,nxpnt
              psrctem(ip)=psrc(ip)
            enddo
            do ih=1,nassy
              do it=1,6
                inpt(it)=neigpt(it,ih)
                phipt(it)=pflxt(inpt(it))/xscdf(it,ig,ih,iz)    ! 2014_05_15 . scb adf hex
              enddo
              do it=1,3
! 2014_05_15 . scb adf hex                
                !dat1=pflxt(inpt(mp1(it)))+pflxt(inpt(mp5(it)))
                !dat2=pflxt(inpt(mp2(it)))+pflxt(inpt(mp4(it)))
                !psrctem(inpt(it))=psrctem(inpt(it))              &
                !             +xsdmwt(3,ih)*pflxt(inpt(mp3(it)))  &
                !             +xsdmwt(1,ih)*dat1+xsdmwt(2,ih)*dat2
                !psrctem(inpt(mp3(it)))=psrctem(inpt(mp3(it)))    &
                !             +xsdmwt(3,ih)*pflxt(inpt(it))       &
                !             +xsdmwt(1,ih)*dat2+xsdmwt(2,ih)*dat1
                dat1=phipt(mp1(it))+phipt(mp5(it))
                dat2=phipt(mp2(it))+phipt(mp4(it))
                psrctem(inpt(it))=psrctem(inpt(it))              &
                             +xsdmwt(3,ih)*phipt(mp3(it))  &
                             +xsdmwt(1,ih)*dat1+xsdmwt(2,ih)*dat2
                psrctem(inpt(mp3(it)))=psrctem(inpt(mp3(it)))    &
                             +xsdmwt(3,ih)*phipt(it)       &
                             +xsdmwt(1,ih)*dat2+xsdmwt(2,ih)*dat1
! added end                
              enddo
            enddo
            do ip=1,nxpnt
              pflxt(ip)=0.7*psrctem(ip)/pcoef(ip)+0.3*pflxt(ip)
            enddo
          enddo
          do ip=1,nxpnt
            pflx(ig,ip,iz)=pflxt(ip)
          enddo
        enddo
      enddo

      return
    end
