! added in ARTOS ver. 0.2 . 2012_07_24 by SCB        
    subroutine solpflx_sp3

! point flux solver using CPB relation

	    use const
	    use geomhex
	    use tpen_sp3
	    use xsec, only : xsdf2, xsdf
	    implicit none

! 2013_05_22 . scb      
      !real,allocatable,save,dimension(:,:) :: psrc,pcoef,psrctem,pflxt, &
      !                                        xsdmwt(:,:,:)
      real,pointer,dimension(:,:) :: psrc,pcoef,psrctem,pflxt,xsdmwt(:,:,:)
      common / rsthexsp3 / psrc,pcoef,psrctem,pflxt,xsdmwt
! added end      
      
      real :: xsdmwtb(2), wts1(2), wts2(2), wts3(2), wts4(2), &
               bsflx(2,6), dathflx(2)
	    real :: radfl, temp, adfl, reflratl, dat1, dat2, dat3, pbdvl
	    integer :: inpt(6), iz, ig, ip, ih, it, nn, iter, im
      real :: jsum(2)
	    logical,save :: first=TRUE
      

      !if(first) then
      if(first .and. .not.associated(psrc)) then  ! 2013_05_22 . scb
        allocate(psrc(2,ncorn),pcoef(2,ncorn),psrctem(2,ncorn),pflxt(2,ncorn))
        allocate(xsdmwt(2,3,nassy))     
        psrc=0.d0
        pcoef=0.d0
        psrctem=0.d0
        pflxt=0.d0
        xsdmwt=0.d0
        first=.FALSE.
      endif

      do iz=1,nz
         do ig=1,ng
            do ip=1,nxpnt   
               do im=1,2
                  psrc(im,ip)=0.
                  pflxt(im,ip)=pflx2(im,ig,ip,iz)
                  pcoef(im,ip)=0
               enddo ! im  
            enddo ! ip
            do ih=1,nassy
               xsdmwtb(1)=wtass(ih)*xsdf(ig,ih,iz)
               xsdmwtb(2)=wtass(ih)*xsdf2(ig,ih,iz)
               do im=1,2
                  temp=xsdmwtb(im)
                  xsdmwt(im,1,ih)=chlval(1)*temp
                  xsdmwt(im,2,ih)=chlval(2)*temp
                  xsdmwt(im,3,ih)=chlval(3)*temp
                  wts1(im)=temp*chlval(4)
                  wts2(im)=temp*chlval(5)
                  wts3(im)=temp*chlval(6)
                  wts4(im)=temp*chlval(7)
               enddo ! im
               do it=1,6
                  nn=neignd(it,ih)
                  if(nn.eq.0) then                
                     if(alxr.eq.0) then
                        jsum(1)=2.*cnto2(1,ig,it,ih,iz)
                        jsum(2)=2.*cnto2(2,ig,it,ih,iz)
                     else
                        jsum(1)=cnto2(1,ig,it,ih,iz)
                        jsum(2)=cnto2(2,ig,it,ih,iz)  
                     endif
                  else
                     jsum(1)=cnto2(1,ig,it,ih,iz)+cnto2(1,ig,neigjin(it,ih),nn,iz)   
                     jsum(2)=cnto2(2,ig,it,ih,iz)+cnto2(2,ig,neigjin(it,ih),nn,iz)     
                  endif

                  bsflx(1,it)=(56.*jsum(1)+24.*jsum(2))/25.
                  bsflx(2,it)=(8.*jsum(1)+32.*jsum(2))/25. 
               enddo ! it
               do im=1,2
                  dathflx(im)=wts4(im)*hflxf2(im,ig,ih,iz)
               enddo

               do it=1,3          
                  do im=1,2
                     nn=neigpt(it,ih)
                     dat1=bsflx(im,it)+bsflx(im,mp5(it))
                     dat2=bsflx(im,mp2(it))+bsflx(im,mp3(it))
                     dat3=wts3(im)*(bsflx(im,mp1(it))+bsflx(im,mp4(it)))
                     pbdvl=pbdv(ig,nn)*codpnt(nn)*wtass(ih)
                     psrc(im,nn)=psrc(im,nn)+(wts1(im)*dat1+wts2(im)*dat2+dat3+dathflx(im))
                     pcoef(im,nn)=pcoef(im,nn)+(xsdmwtb(im)+pbdvl)
                     nn=neigpt(mp3(it),ih)
                     pbdvl=pbdv(ig,nn)*codpnt(nn)*wtass(ih)
                     psrc(im,nn)=psrc(im,nn)+(wts1(im)*dat2+wts2(im)*dat1+dat3+dathflx(im))
                     pcoef(im,nn)=pcoef(im,nn)+(xsdmwtb(im)+pbdvl)
                  enddo ! im
               enddo ! it
            enddo ! ih

            do iter=1, nitrpfc
               do ip=1,nxpnt
                  do im=1,2
                     psrctem(im,ip)=psrc(im,ip)
                  enddo
               enddo
               do ih=1,nassy
                  do it=1,6
                     inpt(it)=neigpt(it,ih)
                  enddo
                  do it=1,3
                     do im=1,2
                        dat1=pflxt(im,inpt(mp1(it)))+pflxt(im,inpt(mp5(it)))
                        dat2=pflxt(im,inpt(mp2(it)))+pflxt(im,inpt(mp4(it)))
                        psrctem(im,inpt(it))=psrctem(im,inpt(it))      &
                             +xsdmwt(im,3,ih)*pflxt(im,inpt(mp3(it)))  &
                             +xsdmwt(im,1,ih)*dat1+xsdmwt(im,2,ih)*dat2
                        psrctem(im,inpt(mp3(it)))=psrctem(im,inpt(mp3(it)))    &
                             +xsdmwt(im,3,ih)*pflxt(im,inpt(it))       &
                             +xsdmwt(im,1,ih)*dat2+xsdmwt(im,2,ih)*dat1
                     enddo
                  enddo ! it
               enddo ! ih
               do ip=1,nxpnt
                  do im=1,2
                     pflxt(im,ip)=0.7*psrctem(im,ip)/pcoef(im,ip)+0.3*pflxt(im,ip)
                  enddo
               enddo
            enddo ! iter
            do ip=1,nxpnt
               do im=1,2
                  pflx2(im,ig,ip,iz)=pflxt(im,ip)
               enddo
            enddo ! ip
         enddo ! ig
      enddo ! iz

      return
    end subroutine
