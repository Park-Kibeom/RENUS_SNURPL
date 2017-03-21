! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
   subroutine upddhathex2g_sp3

#define sp3_nodal_correction   ! 2014_03_12 . scb      

! Update CMFD coefficient(dhat and beta) from Nodal Solution 

      use const
      use geomhex
      use tpen, only : fohflx, fhflx, focnto, fcnto, focntzo, fcntzo
      use tpen_sp3
      use cmfdhex2g
      use xsec
      use sfam, only : fphi, phi, fphi2

      implicit none

! variables for nem solver in z-direction
! atleak : transverse leakage for axial direction
! atleaku,atleakl : transverse leakage for upper and lower node
! xsdu,xsdl : xsd(diffusion coefficients) for upper and lower node

      integer :: mvz(2)
      data mvz/2,1/

      integer :: iz, ih, ig, ig2, it, isfc, ind, &
	              is1, is2, id1, id2, iz1, iz2, im
      real :: flx1, flx2, dfdm, w, fdflxs, adfl, albl,  &
	           reflratfl, hzr, xsd0
      real :: reflratzf(ng)
      real :: sumflx, sumcnt, dc, sumflx2, sumcntz
      real :: curo1, curo2, jm(2), jp(2), sflx(2)

      integer :: m, l, k, md, mbeg, mend, m2
      real :: xss2t, xssup



      if(ng.ne.2) then
         do iz=1,nz 
            do ih=1,nassy
               do ig=1,ng2
                  sumflx=0
                  do ig2=mgb(ig),mge(ig)
                     sumflx=sumflx+hflxf2(1,ig2,ih,iz)
                  enddo               
                  hflx2(1,ig,ih,iz)=sumflx
#ifdef sp3_nodal_correction   ! 2014_03_12 . scb                      
                  phi(ig,ih,iz)=hflx2(1,ig,ih,iz)   ! 2014_03_10 . scb
#endif                  

                  sumflx2=0
                  do ig2=mgb(ig),mge(ig)
                     sumflx2=sumflx2+hflxf2(2,ig2,ih,iz)
                  enddo               
                  hflx2(2,ig,ih,iz)=sumflx2

                  do ig2=mgb(ig),mge(ig)
                     fohflx(ig2,ih,iz)=fhflx(ig2,ih,iz)
                     fhflx(ig2,ih,iz)=hflxf2(1,ig2,ih,iz)/sumflx
                     fphi(ig2,ih,iz)=fhflx(ig2,ih,iz)
                     fphi2(ig2,ih,iz)=hflxf2(2,ig2,ih,iz)/sumflx2
                  enddo
               enddo ! ig
               do it=1,6
                  do ig=1,ng2
                     sumcnt=0
                     do ig2=mgb(ig),mge(ig)
                        sumcnt=sumcnt+cnto2(1,ig2,it,ih,iz)
                     enddo
                     do ig2=mgb(ig),mge(ig)
                        focnto(ig2,it,ih,iz)=fcnto(ig2,it,ih,iz)
                        fcnto(ig2,it,ih,iz)=cnto2(1,ig2,it,ih,iz)/sumcnt
                     enddo
                  enddo
               enddo ! it
            enddo
         enddo
         
         if(nz.gt.1) then
            do iz=1,nz
               do ih=1,nassy
                  do it=1,2
                     do ig=1,ng2
                        sumcntz=0
                        do ig2=mgb(ig),mge(ig)
                           sumcntz=sumcntz+cntzo2(1,ig2,it,ih,iz)
                        enddo
                        do ig2=mgb(ig),mge(ig)
                           focntzo(ig2,it,ih,iz)=fcntzo(ig2,it,ih,iz)
                           fcntzo(ig2,it,ih,iz)=cntzo2(1,ig2,it,ih,iz)/sumcntz
                        enddo
                     enddo
                  enddo
               enddo ! ih
            enddo ! iz
         endif
         
! 2014_03_10 . scb   
#ifdef sp3_nodal_correction   ! 2014_03_12 . scb    
         call colxs(fphi)
         call colxs_sp3(fphi) 
         call upddtilhex2g_sp3
#endif         
! added end         
      endif

! update D(diffusion coefficient) for CMFD
      do isfc=1,nxsfc
         ind=neigsnd(3,isfc)
         if(ind.eq.12) then
            is1=neigsnd(1,isfc)
            is2=neigsnd(2,isfc)
            id1=neigsnd(4,isfc)
            id2=neigsnd(5,isfc)
            do iz=1,nz
               do ig=1,ng2
                  dc=xsd(ig,is1,iz)
                  do im=1,2
                     curo1=0.
                     curo2=0.
                     do ig2=mgb(ig),mge(ig)
                        curo1=curo1+cnto2(im,ig2,id1,is1,iz)
                        curo2=curo2+cnto2(im,ig2,id2,is2,iz)
                     enddo
                     jp(im)=curo1+curo2
                     jm(im)=curo1-curo2
                  enddo
                  sflx(1)=(40.*jp(1)-40.*jp(2))/25.
                  sflx(2)=(8.*jp(1)+32.*jp(2))/25.  
                               
                  flx1=hflx2(1,ig,is1,iz)
                  flx2=hflx2(1,ig,is2,iz)

                  dfdm=dfd(ig,isfc,iz)
                  w=dfdm/dc/2.
                  fdflxs=(1.-w)*flx1+w*flx2 
                  dhat(ig,isfc,iz)=-(rt3*hside*jm(1)+dfdm*(flx2-flx1))/(flx2+flx1)
                  betaphis(ig,isfc,iz)=2.*(sflx(1)-fdflxs)/(flx1+flx2)


                  flx1=hflx2(2,ig,is1,iz)
                  flx2=hflx2(2,ig,is2,iz)

                  dc=xsd2(ig,is1,iz)
                  dfdm=dfd2(ig,isfc,iz)
                  w=dfdm/dc/2.
                  fdflxs=(1.-w)*flx1+w*flx2 

                  dhat2(ig,isfc,iz)=-(rt3*hside*jm(2)+dfdm*(flx2-flx1))/(flx2+flx1)
                  betaphis2(ig,isfc,iz)=2.*(sflx(2)-fdflxs)/(flx1+flx2)
               enddo ! ig
            enddo ! iz
         else
            is1=neigsnd(ind,isfc)
            id1=neigsnd(ind+3,isfc)
            do iz=1,nz 
               do ig=1,ng2
                  do im=1,2
                     curo1=0.
                     do ig2=mgb(ig),mge(ig)
                        curo1=curo1+cnto2(im,ig2,id1,is1,iz)
                     enddo
                     if(alxr.eq.0) then
                        jp(im)=2.*curo1
                        jm(im)=0.
                     else
                        jp(im)=curo1
                        jm(im)=curo1
                     endif
                  enddo

                  sflx(1)=(40.*jp(1)-40.*jp(2))/25.
                  sflx(2)=(8.*jp(1)+32.*jp(2))/25. 
                  
                  flx1=hflx2(1,ig,is1,iz)
                  dfdm=dfd(ig,isfc,iz)
                  dhat(ig,isfc,iz)=(dfdm*flx1-rt3*hside*jm(1))/flx1
                  if(ind.eq.2) dhat(ig,isfc,iz)=-dhat(ig,isfc,iz)
                  betaphis(ig,isfc,iz)=sflx(1)/flx1

                  flx1=hflx2(2,ig,is1,iz)
                  dfdm=dfd2(ig,isfc,iz)
                  dhat2(ig,isfc,iz)=(dfdm*flx1-rt3*hside*jm(2))/flx1
                  if(ind.eq.2) dhat2(ig,isfc,iz)=-dhat2(ig,isfc,iz)
                  betaphis2(ig,isfc,iz)=sflx(2)/flx1
               enddo ! ig
            enddo ! iz
         endif ! ind
      enddo ! isfc

      if(nz.eq.1) return
! update D in z-direction for CMFD
      do iz=1,nz+1 
         ind=neigsndz(3,iz)
         if(ind.eq.12) then
            iz1=neigsndz(1,iz)
            iz2=neigsndz(2,iz)
            id1=neigsndz(4,iz)
            id2=neigsndz(5,iz)
            hzr=hz(iz2)/hz(iz1)
            do ih=1,nassy 
               do ig=1,ng2
                  dc=xsd(ig,ih,iz1)
                  do im=1,2
                     curo1=0.
                     curo2=0.
                     do ig2=mgb(ig),mge(ig)
                        curo1=curo1+cntzo2(im,ig2,id1,ih,iz1)
                        curo2=curo2+cntzo2(im,ig2,id2,ih,iz2)
                     enddo
                     jp(im)=curo1+curo2
                     jm(im)=curo1-curo2
                  enddo
                  sflx(1)=(40.*jp(1)-40.*jp(2))/25.
                  sflx(2)=(8.*jp(1)+32.*jp(2))/25. 
                                                                                 
                  flx1=hflx2(1,ig,ih,iz1)  
                  flx2=hflx2(1,ig,ih,iz2) 

                  dfdm=dfdz(ig,ih,iz)
                  w=dfdm/dc/(1.+hzr)
                  fdflxs=(1.-w)*flx1+w*flx2 
                  dhatz(ig,ih,iz)=-((hz(iz1)+hz(iz2))*jm(1)/2. &
                                 +dfdm*(flx2-flx1))/(flx2+flx1)
                  betaphisz(ig,ih,iz)=2.*(sflx(1)-fdflxs)/(flx1+flx2)


                  flx1=hflx2(2,ig,ih,iz1)  
                  flx2=hflx2(2,ig,ih,iz2) 

                  dc=xsd2(ig,ih,iz1)
                  dfdm=dfdz2(ig,ih,iz)
                  w=dfdm/dc/(1.+hzr)
                  fdflxs=(1.-w)*flx1+w*flx2 
                  betaphisz2(ig,ih,iz)=2.*(sflx(2)-fdflxs)/(flx1+flx2)
               enddo
            enddo
         else
            iz1=neigsndz(ind,iz)
            id1=neigsndz(ind+3,iz)
            do ih=1,nassy 
               do ig=1,ng2
                  do im=1,2
                     curo1=0.
                     do ig2=mgb(ig),mge(ig)
                        curo1=curo1+cntzo2(im,ig2,id1,ih,iz1)
                     enddo                    
                     if(iz.eq.1) then
                        if(alzl.eq.0) then
                           jp(im)=2.*curo1
                           jm(im)=0.
                        else
                           jp(im)=curo1
                           jm(im)=curo1
                        endif
                     else
                        if(alzr.eq.0) then
                           jp(im)=2.*curo1
                           jm(im)=0.
                        else
                           jp(im)=curo1
                           jm(im)=curo1
                        endif
                     endif                                                                                                                                      
                  enddo

                  sflx(1)=(40.*jp(1)-40.*jp(2))/25.
                  sflx(2)=(8.*jp(1)+32.*jp(2))/25.

                  flx1=hflx2(1,ig,ih,iz1)  
                  dfdm=dfdz(ig,ih,iz)
                  dhatz(ig,ih,iz)=(dfdm*flx1-hz(iz1)*jm(1))/flx1
                  if(ind.eq.2) dhatz(ig,ih,iz)=-dhatz(ig,ih,iz)
                  betaphisz(ig,ih,iz)=sflx(1)/flx1

                  flx1=hflx2(2,ig,ih,iz1)  
                  betaphisz2(ig,ih,iz)=sflx(2)/flx1
               enddo
            enddo
         endif
      enddo

      return
   end subroutine
