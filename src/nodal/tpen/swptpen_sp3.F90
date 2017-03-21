! added in ARTOS ver. 0.2 . 2012_07_24 by SCB        
   subroutine swptpen_sp3(iftran)

      use const
      use geomhex
      use tpen_sp3
      use xsec
      use sfam, only : reigv, phi_p2
      
      implicit none

#define paral_sp3_tpen

#ifdef paral_sp3_tpen
! boundary conditions for TPEN solver in r-direction
      real :: cnti(2,ng,6),pbflx(2,ng,6)

      integer :: it, ig, ih, iz, nn, im, ig2, m, iter, iin
      real :: asrc(2,6), xsrc(2,6), ysrc(2,6)
      real :: asrcd(2,6), xsrcd(2,6), ysrcd(2,6)
      real :: reflr, reflz

      real :: areavol, sum, rhz
      real :: srczn(2,ng,6)

      real :: srcbal(2,ng,6), srcmomx(2,ng,6), srcmomy(2,ng,6), srczm
      logical,save :: first=TRUE
      
      logical :: iftran   ! 2012_12_05 . scb


! axial source for TPEN cal.
      if(nz.eq.1) goto 100   ! 2012_12_18 . scb 
! redo this nz.eq.1 statements because of some bugs... ! 2013_09_24 . scb
      iz=1
      rhz=1./hz(iz)
      do ih=1,nassy
         do ig=1,ng
            do im=1,2
               srcz2(im,ig,ih,iz)= &
                   +((reflratzb-1.)*cntzo2(im,ig,1,ih,iz) &
                   +cntzo2(im,ig,1,ih,iz+1)-cntzo2(im,ig,2,ih,iz))*rhz
            enddo
         enddo
      enddo
      do iz=2,nz-1
         rhz=1./hz(iz)
         do ih=1,nassy
            do ig=1,ng
               do im=1,2
                 ! incoming - outgoing !! 2013_08_12 . scb comment
                  srcz2(im,ig,ih,iz)= &
                      +(cntzo2(im,ig,2,ih,iz-1)+cntzo2(im,ig,1,ih,iz+1) &
                      -cntzo2(im,ig,1,ih,iz)-cntzo2(im,ig,2,ih,iz))*rhz
               enddo
            enddo
         enddo
      enddo
      iz=nz
      rhz=1./hz(iz)
      do ih=1,nassy
         do ig=1,ng
            do im=1,2
               srcz2(im,ig,ih,iz)= &
                   +(cntzo2(im,ig,2,ih,iz-1)-cntzo2(im,ig,1,ih,iz) &
                   +(reflratzt-1.)*cntzo2(im,ig,2,ih,iz))*rhz
            enddo
         enddo
      enddo  
      
! prepare boundary condition(incoming J(cnti) and point flux(pbflx))
100 continue
      
! 2012_12_07 . scb      
      if(iftran) then
        do iz=1,nz
          do ih=1,nassy
            do ig=1,ng
              do im=1,2
                srcz2(im,ig,ih,iz)=srcz2(im,ig,ih,iz)+sefve2(im,ig,ih,iz)
              enddo
            enddo
          enddo
        enddo        
      endif
! added end    
    
!!$OMP PARALLEL DO  private(iz,ih,it,ig,nn,reflr,im,cnti,srczn,pbflx,srczm,srcbal,srcmomx,srcmomy,asrc,xsrc,ysrc,asrcd,xsrcd,ysrcd)    
      do iz=1,nz
         do ih=1,nassy
            ! preparation for radial solver(cnti,srczn,pbflx)
            !if(nassy.gt.1) then   ! 2012_12_18 . scb
               do it=1,6
                  nn=neignd(it,ih)
                  if(nn.eq.0) then
                     do ig=1,ng
                        ! reflrhom: reflective=1, vacuum=0
                        reflr=(1-2*alxr)/(1.+2*alxr)      
                        do im=1,2
                           cnti(im,ig,it)=cnto2(im,ig,it,ih,iz)*reflr
                        enddo

                        if(ifbcref) then
                           do im=1,2
                              srczn(im,ig,it)=srcz2(im,ig,ih,iz)
                           enddo
                        else 
                           do im=1,2
                              srczn(im,ig,it)=0.
                              !srczn(im,ig,it)=-srcz2(im,ig,ih,iz) ! fix_tpen
                           enddo
                        endif
                     enddo ! ig
                  else
                     do ig=1,ng
                        do im=1,2
                           cnti(im,ig,it)=cnto2(im,ig,neigjin(it,ih),nn,iz)
                           srczn(im,ig,it)=srcz2(im,ig,nn,iz)
                        enddo
                     enddo ! ig
                  endif ! nn
               enddo ! it

               do it=1,6 
                  nn=neigpt(it,ih)
                  do ig=1,ng
                     pbflx(1,ig,it)=pflx2(1,ig,nn,iz)
                     pbflx(2,ig,it)=pflx2(2,ig,nn,iz)
                  enddo
               enddo
            !endif   ! 2012_12_18 . scb

            if(nz.gt.1) then   ! 2012_12_18 . scb
               do it=1,6
                  do ig=1,ng
                     do im=1,2
                        srczm=srcz2(im,ig,ih,iz) 
                        srcbal(im,ig,it)=srczm                                 &
                           +(83.*srczn(im,ig,it)+17.*srczn(im,ig,mp1(it))        &
                           -37.*srczn(im,ig,mp2(it))-43.*srczn(im,ig,mp3(it))   &
                           -37.*srczn(im,ig,mp4(it))+17.*srczn(im,ig,mp5(it)))/540. 
                        srcmomx(im,ig,it)=(-60.*srczm                          &
                           +59.*srczn(im,ig,it)+14.*srczn(im,ig,mp1(it))        &
                           -10.*srczn(im,ig,mp2(it))-7.*srczn(im,ig,mp3(it))    & 
                           -10.*srczn(im,ig,mp4(it))+14.*srczn(im,ig,mp5(it)))/3240. 
                        srcmomy(im,ig,it)=-(srczn(im,ig,mp1(it))-srczn(im,ig,mp5(it)))/40.  &  !m-g
                           -(srczn(im,ig,mp2(it))-srczn(im,ig,mp4(it)))/360. 
                     enddo         
                  enddo
               enddo
   ! 2012_12_18 . scb              
            else
               srcbal=0.
               srcmomx=0.
               srcmomy=0.
            endif

            ! response matrix solver          
            do it=1,6
               asrcd(1,it)=0.
               xsrcd(1,it)=0.
               ysrcd(1,it)=0.
               do m=1,ng
                  asrcd(1,it)=asrcd(1,it)+xsnff(m,ih,iz)*(aflx2(1,m,it,ih,iz)-2.*aflx2(2,m,it,ih,iz))
                  xsrcd(1,it)=xsrcd(1,it)+xsnff(m,ih,iz)*(xmom2(1,m,it,ih,iz)-2.*xmom2(2,m,it,ih,iz))
                  ysrcd(1,it)=ysrcd(1,it)+xsnff(m,ih,iz)*(ymom2(1,m,it,ih,iz)-2.*ymom2(2,m,it,ih,iz))
               enddo
            enddo

            do iter=1,1
               do ig=1,ng
                  do it=1,6
                     asrc(1,it)=reigv*xschif(ig,ih,iz)*asrcd(1,it)
                     xsrc(1,it)=reigv*xschif(ig,ih,iz)*xsrcd(1,it)
                     ysrc(1,it)=reigv*xschif(ig,ih,iz)*ysrcd(1,it)
                     do ig2=xssfs(ig,ih,iz),xssfe(ig,ih,iz)
                        asrc(1,it)=asrc(1,it)+xssf(ig2,ig,ih,iz)*(aflx2(1,ig2,it,ih,iz)-2.*aflx2(2,ig2,it,ih,iz))
                        xsrc(1,it)=xsrc(1,it)+xssf(ig2,ig,ih,iz)*(xmom2(1,ig2,it,ih,iz)-2.*xmom2(2,ig2,it,ih,iz))
                        ysrc(1,it)=ysrc(1,it)+xssf(ig2,ig,ih,iz)*(ymom2(1,ig2,it,ih,iz)-2.*ymom2(2,ig2,it,ih,iz))
                     enddo
                     asrc(2,it)=-2./3.*asrc(1,it) +srcbal(2,ig,it)
                     xsrc(2,it)=-2./3.*xsrc(1,it) +srcmomx(2,ig,it)
                     ysrc(2,it)=-2./3.*ysrc(1,it) +srcmomy(2,ig,it)

                     asrc(1,it)=asrc(1,it) +srcbal(1,ig,it)
                     xsrc(1,it)=xsrc(1,it) +srcmomx(1,ig,it)
                     ysrc(1,it)=ysrc(1,it) +srcmomy(1,ig,it)
                  enddo

                  ! solve the problem by TPEN
                  call onetpen_sp3( &
                     !xstf(ig,ih,iz),xstrf(ig,ih,iz), &
                     xstf(ig,ih,iz),xstfsp3(ig,ih,iz), &   ! 2014_10_07 . scb
                     xsdf(ig,ih,iz),xsdf2(ig,ih,iz), &
                     hside,pbflx(:,ig,:),cnti(:,ig,:),                                         &
                     aflx2(:,ig,:,ih,iz),xmom2(:,ig,:,ih,iz),ymom2(:,ig,:,ih,iz), &
                     cnto2(:,ig,:,ih,iz),hflxf2(:,ig,ih,iz),   &
                     asrc,xsrc,ysrc,       &
                     sfli2(:,ig,:,ih,iz))
               enddo ! ig
            enddo ! iter     
         enddo ! ih
      enddo ! iz
!!$OMP END PARALLEL DO        
   
      return
   end subroutine
   
   
   
   
#else
! boundary conditions for TPEN solver in r-direction
      real :: cnti(2,ng,6),pbflx(2,ng,6)

      integer :: it, ig, ih, iz, nn, im, ig2, m, iter, iin
      real :: asrc(2,6), xsrc(2,6), ysrc(2,6)
      real :: asrcd(2,6), xsrcd(2,6), ysrcd(2,6)
      real :: reflr, reflz

      real :: areavol, sum, rhz
      real :: srczn(2,ng,6)

      real :: srcbal(2,ng,6), srcmomx(2,ng,6), srcmomy(2,ng,6), srczm
      logical,save :: first=TRUE
      
      logical :: iftran   ! 2012_12_05 . scb


! axial source for TPEN cal.
      iz=1
      rhz=1./hz(iz)
      do ih=1,nassy
         do ig=1,ng
            do im=1,2
               srcz2(im,ig,ih,iz)= &
                   +((reflratzb-1.)*cntzo2(im,ig,1,ih,iz) &
                   +cntzo2(im,ig,1,ih,iz+1)-cntzo2(im,ig,2,ih,iz))*rhz
            enddo
         enddo
      enddo
      do iz=2,nz-1
         rhz=1./hz(iz)
         do ih=1,nassy
            do ig=1,ng
               do im=1,2
                  srcz2(im,ig,ih,iz)= &
                      +(cntzo2(im,ig,2,ih,iz-1)+cntzo2(im,ig,1,ih,iz+1) &
                      -cntzo2(im,ig,1,ih,iz)-cntzo2(im,ig,2,ih,iz))*rhz
               enddo
            enddo
         enddo
      enddo
      iz=nz
      rhz=1./hz(iz)
      do ih=1,nassy
         do ig=1,ng
            do im=1,2
               srcz2(im,ig,ih,iz)= &
                   +(cntzo2(im,ig,2,ih,iz-1)-cntzo2(im,ig,1,ih,iz) &
                   +(reflratzt-1.)*cntzo2(im,ig,2,ih,iz))*rhz
            enddo
         enddo
      enddo
      
! 2012_12_07 . scb      
      if(iftran) then
        do iz=1,nz
          do ih=1,nassy
            do ig=1,ng
              do im=1,2
                srcz2(im,ig,ih,iz)=srcz2(im,ig,ih,iz)+sefve2(im,ig,ih,iz)
              enddo
            enddo
          enddo
        enddo        
      endif
! added end      
      
! prepare boundary condition(incoming J(cnti) and point flux(pbflx))
100   continue
      do iz=1,nz
         do ih=1,nassy
            ! preparation for radial solver(cnti,srczn,pbflx)
            !if(nassy.gt.1) then   ! 2012_12_18 . scb
               do it=1,6
                  nn=neignd(it,ih)
                  if(nn.eq.0) then
                     do ig=1,ng
                        ! reflrhom: reflective=1, vacuum=0
                        reflr=(1-2*alxr)/(1.+2*alxr)      
                        do im=1,2
                           cnti(im,ig,it)=cnto2(im,ig,it,ih,iz)*reflr
                        enddo

                        if(ifbcref) then
                           do im=1,2
                              srczn(im,ig,it)=srcz2(im,ig,ih,iz)
                           enddo
                        else 
                           do im=1,2
                              srczn(im,ig,it)=0.
                              !srczn(im,ig,it)=-srcz2(im,ig,ih,iz) ! fix_tpen
                           enddo
                        endif
                     enddo ! ig
                  else
                     do ig=1,ng
                        do im=1,2
                           cnti(im,ig,it)=cnto2(im,ig,neigjin(it,ih),nn,iz)
                           srczn(im,ig,it)=srcz2(im,ig,nn,iz)
                        enddo
                     enddo ! ig
                  endif ! nn
               enddo ! it

               do it=1,6 
                  nn=neigpt(it,ih)
                  do ig=1,ng
                     pbflx(1,ig,it)=pflx2(1,ig,nn,iz)
                     pbflx(2,ig,it)=pflx2(2,ig,nn,iz)
                  enddo
               enddo
            !endif   ! 2012_12_18 . scb

            !if(nz.gt.1) then   ! 2012_12_18 . scb
               do it=1,6
                  do ig=1,ng
                     do im=1,2
                        srczm=srcz2(im,ig,ih,iz) 
                        srcbal(im,ig,it)=srczm                                 &
                           +(83.*srczn(im,ig,it)+17.*srczn(im,ig,mp1(it))        &
                           -37.*srczn(im,ig,mp2(it))-43.*srczn(im,ig,mp3(it))   &
                           -37.*srczn(im,ig,mp4(it))+17.*srczn(im,ig,mp5(it)))/540. 
                        srcmomx(im,ig,it)=(-60.*srczm                          &
                           +59.*srczn(im,ig,it)+14.*srczn(im,ig,mp1(it))        &
                           -10.*srczn(im,ig,mp2(it))-7.*srczn(im,ig,mp3(it))    & 
                           -10.*srczn(im,ig,mp4(it))+14.*srczn(im,ig,mp5(it)))/3240. 
                        srcmomy(im,ig,it)=-(srczn(im,ig,mp1(it))-srczn(im,ig,mp5(it)))/40.  &  !m-g
                           -(srczn(im,ig,mp2(it))-srczn(im,ig,mp4(it)))/360. 
                     enddo         
                  enddo
               enddo
   ! 2012_12_18 . scb              
            !else
            !   srcbal=0.
            !   srcmomx=0.
            !   srcmomy=0.
            !endif

            ! response matrix solver          
            do it=1,6
               asrcd(1,it)=0.
               xsrcd(1,it)=0.
               ysrcd(1,it)=0.
               do m=1,ng
                  asrcd(1,it)=asrcd(1,it)+xsnff(m,ih,iz)*(aflx2(1,m,it,ih,iz)-2.*aflx2(2,m,it,ih,iz))
                  xsrcd(1,it)=xsrcd(1,it)+xsnff(m,ih,iz)*(xmom2(1,m,it,ih,iz)-2.*xmom2(2,m,it,ih,iz))
                  ysrcd(1,it)=ysrcd(1,it)+xsnff(m,ih,iz)*(ymom2(1,m,it,ih,iz)-2.*ymom2(2,m,it,ih,iz))
               enddo
            enddo

            do iter=1,1
               do ig=1,ng
                  do it=1,6
                     asrc(1,it)=reigv*xschif(ig,ih,iz)*asrcd(1,it)
                     xsrc(1,it)=reigv*xschif(ig,ih,iz)*xsrcd(1,it)
                     ysrc(1,it)=reigv*xschif(ig,ih,iz)*ysrcd(1,it)
                     do ig2=xssfs(ig,ih,iz),xssfe(ig,ih,iz)
                        asrc(1,it)=asrc(1,it)+xssf(ig2,ig,ih,iz)*(aflx2(1,ig2,it,ih,iz)-2.*aflx2(2,ig2,it,ih,iz))
                        xsrc(1,it)=xsrc(1,it)+xssf(ig2,ig,ih,iz)*(xmom2(1,ig2,it,ih,iz)-2.*xmom2(2,ig2,it,ih,iz))
                        ysrc(1,it)=ysrc(1,it)+xssf(ig2,ig,ih,iz)*(ymom2(1,ig2,it,ih,iz)-2.*ymom2(2,ig2,it,ih,iz))
                     enddo
                     asrc(2,it)=-2./3.*asrc(1,it) +srcbal(2,ig,it)
                     xsrc(2,it)=-2./3.*xsrc(1,it) +srcmomx(2,ig,it)
                     ysrc(2,it)=-2./3.*ysrc(1,it) +srcmomy(2,ig,it)

                     asrc(1,it)=asrc(1,it) +srcbal(1,ig,it)
                     xsrc(1,it)=xsrc(1,it) +srcmomx(1,ig,it)
                     ysrc(1,it)=ysrc(1,it) +srcmomy(1,ig,it)
                  enddo

                  ! solve the problem by TPEN
                  call onetpen_sp3( &
                     !xstf(ig,ih,iz),xstrf(ig,ih,iz), &
                     xstf(ig,ih,iz),xstfsp3(ig,ih,iz), &    ! 2014_10_07 . scb
                     xsdf(ig,ih,iz),xsdf2(ig,ih,iz), &
                     hside,pbflx(:,ig,:),cnti(:,ig,:),                                         &
                     aflx2(:,ig,:,ih,iz),xmom2(:,ig,:,ih,iz),ymom2(:,ig,:,ih,iz), &
                     cnto2(:,ig,:,ih,iz),hflxf2(:,ig,ih,iz),   &
                     asrc,xsrc,ysrc,       &
                     sfli2(:,ig,:,ih,iz))
               enddo ! ig
            enddo ! iter     
         enddo ! ih
      enddo ! iz
   
      return
    end subroutine
#endif