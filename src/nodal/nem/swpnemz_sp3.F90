! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
  !subroutine swpnemz_sp3(iftran)
  subroutine swpnemz_sp3_old(iftran)

      use const
      use geomhex
      use tpen_sp3
      use cmfdhex2g, only : phisz
      use xsec
      use sfam, only : reigv, phi_p2
      implicit none

#define paral_sp3_nem
#define scb_correction_tl_sp3   ! 2013_07_10 . scb
#define paral_dbg   ! 2014_03_23 . scb

#ifdef paral_sp3_nem
      logical :: iftran  ! 2012_12_05 . scb

      integer :: iz, ih, ig, im, ig2, it
      real :: hzm, hzl, hzu, atleakm(2,ng), atleakl(2,ng), atleaku(2,ng)
      real :: dcm(2,ng), dcu(2,ng), dcl(2,ng)   ! 2012_12_18 . scb
      real :: cntz0i(2,ng), cntz1i(2,ng) 
      integer :: neigdn, neigup, nn

      real :: ss0d, ss1d, ss2d
      real :: ss0, ss1, ss2
      real :: q0(2), q1(2), q2(2)
      integer :: iter

      real :: sum, reflratl
      real :: dc(2), tl0(2,ng), tl1(2,ng), alp1(2,ng), alp2(2,ng)
      real :: sbal(2,ng), smm1(2,ng), smm2(2,ng)

      integer :: ifrom, ito, istp
! 2013_07_10 . scb      
#ifdef scb_correction_tl_sp3
      real :: ratl, ratr, difsl(2,ng), difsr(2,ng)
#endif
! added end

! 2014_03_23 . scb for paral dbg
      logical,save :: first=.true.
      integer :: iodbg
      
      real :: cdum1(2), cdum2(2)
! added end      

! average transverse leakage cal. for z-direction 
      do ih=1,nassy
        do iz=1,nz
          do ig=1,ng  
            do im=1,2
              sum=0
              do it=1,6
                nn=neignd(it,ih)
                if(nn.eq.0) then
                  reflratl=(1.d0-2.d0*alxr)/(1.d0+2.d0*alxr)
                  sum=sum+(1.d0-reflratl)*cnto2(im,ig,it,ih,iz)
                else
                  ! outgoing - incoming current !!  ! 2013_08_12 . scb comment
                  sum=sum+cnto2(im,ig,it,ih,iz)-cnto2(im,ig,neigjin(it,ih),nn,iz)
                endif
              enddo
              atleak2(im,ig,ih,iz)=sum*2.d0*rt3/9.d0/hside 
! 2012_12_05 . scb              
	            if(iftran) then
                !sefve2(im,ig,ih,iz)=0.d0   ! 2014_02_12 . scb for dbg 
	              atleak2(im,ig,ih,iz)=atleak2(im,ig,ih,iz)-sefve2(im,ig,ih,iz)
              endif                  
! added end              
              hflxf2d(im,ig,ih,iz)=hflxf2(im,ig,ih,iz)
              zmom1d(im,ig,ih,iz)=zmom1(im,ig,ih,iz)
              zmom2d(im,ig,ih,iz)=zmom2(im,ig,ih,iz)
            enddo 
          enddo   
        enddo   
      enddo  

!$OMP PARALLEL DO  private(iz,ih,hzm,hzu,hzl,neigdn,neigup,ig,im,atleakm,atleakl,atleaku,cntz0i,cntz1i,dc,tl0,tl1,alp1,alp2,sbal,smm1,smm2,ss0d,ss1d,ss2d,ss0,ss1,ss2,ig2,q0,q1,q2,ratl, ratr, difsl, difsr,iodbg)
      do iz=1,nz
         do ih=1,nassy
            ! preparation for axial solver(cnti,srczn,pbflx)
            hzm=hz(iz) 
            hzu=hzm
            hzl=hzm 
            !if(nz.gt.1) then   ! 2012_12_18 . scb
               neigdn=neigz(1,iz) 
               neigup=neigz(2,iz)
               do ig=1,ng
                  do im=1,2
                     atleakm(im,ig)=atleak2(im,ig,ih,iz) 
                     atleakl(im,ig)=atleakm(im,ig)
                     atleaku(im,ig)=atleakm(im,ig)
                  enddo 
               enddo
               if(neigdn.ne.0) then 
                  hzl=hz(neigdn) 
                  do ig=1,ng
                     do im=1,2
                        cntz0i(im,ig)=cntzo2(im,ig,2,ih,neigdn)
                        atleakl(im,ig)=atleak2(im,ig,ih,neigdn)
                     enddo
                  enddo                       
               else
                  do ig=1,ng
                     do im=1,2
                        cntz0i(im,ig)=reflratzb*cntzo2(im,ig,1,ih,iz)
#ifndef scb_correction_tl_sp3     ! 2013_08_12 . scb
                        if(.not.ifbcrefzb) atleakl(im,ig)=-atleak2(im,ig,ih,iz)   ! 2012_12_14 . scb 
#else                        
                        if(.not.ifbcrefzb) atleakl(im,ig)=0.d0   ! 2012_12_14 . scb
#endif
                     enddo
                  enddo
! 2013_07_10 . scb                  
#ifdef scb_correction_tl_sp3                        
                  if(.not.ifbcrefzb) hzl=0.d0
#endif
! added end
               endif
               if(neigup.ne.0) then 
                  hzu=hz(neigup) 
                  do ig=1,ng
                     do im=1,2
                        cntz1i(im,ig)=cntzo2(im,ig,1,ih,neigup)
                        if(iz.eq.2) write(324,*) im,ig,ih,neigup,cntzo2(im,ig,1,ih,neigup)
                        atleaku(im,ig)=atleak2(im,ig,ih,neigup)
                     enddo  
                  enddo
               else
                  do ig=1,ng
                     do im=1,2
                        cntz1i(im,ig)=reflratzt*cntzo2(im,ig,2,ih,iz)
                        if(iz.eq.2) write(324,*) im,ig,ih,neigup,cntzo2(im,ig,2,ih,neigup)
#ifndef scb_correction_tl_sp3        ! 2013_08_12 . scb
                        if(.not.ifbcrefzt) atleaku(im,ig)=-atleak2(im,ig,ih,iz)    ! 2012_12_14 . scb
#else                        
                        if(.not.ifbcrefzt) atleaku(im,ig)=0.d0   ! 2012_12_14 . scb
#endif                        
                     enddo
                  enddo
! 2013_07_10 . scb                  
#ifdef scb_correction_tl_sp3                        
                  if(.not.ifbcrefzt) hzu=0.d0
#endif
! added end
               endif
#ifdef paral_dbg
              !iodbg=100*iz
              !if(iz.eq.2) then
              !  write(iodbg,*) cntz1i(:,:) 
              !  stop
              !endif              
#endif                        
            !endif   ! 2012_12_18 . scb

! Transverse Leakage and coeff's at interface 
!#ifdef DEBUG      ! 2012_12_18 . scb
#ifdef scb_correction_tl_sp3
            ratl=hzl/hzm
            ratr=hzu/hzm
            
            difsl=atleakl-atleakm
            difsr=atleaku-atleakm
            
            if(hzl.lt.1.e-10) then
              alp1(:,:)=(atleakm(:,:)+difsr(:,:)/((1+ratr)*(1+2*ratr))) / (1+1/(1+2*ratr))
              alp2(:,:)=(-atleakm(:,:)+difsr(:,:)/(1+ratr)) / (2*(1+ratr))
            elseif(hzu.lt.1.e-10) then
              alp1(:,:)=-(atleakm(:,:)+difsl(:,:)/((1+ratl)*(1+2*ratl))) / (1+1/(1+2*ratl))
              alp2(:,:)=(-atleakm(:,:)+difsl(:,:)/(1+ratl)) / (2*(1+ratl))     
            else
              alp1(:,:)=(-difsl(:,:)/((1+ratl)*(1+2*ratl))+difsr(:,:)/((1+ratr)*(1+2*ratr))) / (1/(1+2*ratl)+1/(1+2*ratr))
              alp2(:,:)=(difsl(:,:)/(1+ratl)+difsr(:,:)/(1+ratr)) / (2*(1+ratl+ratr))        
            endif
            do ig=1,ng
               do im=1,2          
                  sbal(im,ig)=-atleakm(im,ig)
                  smm1(im,ig)=-alp1(im,ig)/3.d0
                  smm2(im,ig)=-alp2(im,ig)/5.d0
               enddo
            enddo
#else
            do ig=1,ng
               dc(1)=xsdf(ig,ih,iz)
               dc(2)=xsdf2(ig,ih,iz)

               do im=1,2
                  tl0(im,ig)=(atleakl(im,ig)/hzl+atleakm(im,ig)/hzm) &
                        /(dc(im)/hzl+dc(im)/hzm)
                  tl1(im,ig)=(atleaku(im,ig)/hzu+atleakm(im,ig)/hzm) &
                        /(dc(im)/hzu+dc(im)/hzm)

                  alp1(im,ig)=(tl1(im,ig)-tl0(im,ig))/2.
                  alp2(im,ig)=atleakm(im,ig)/dc(im)-(tl1(im,ig)+tl0(im,ig))/2.

                  sbal(im,ig)=-atleakm(im,ig)
                  smm1(im,ig)=-alp1(im,ig)*dc(im)/3.
                  smm2(im,ig)=-alp2(im,ig)*dc(im)/5.
               enddo
            enddo
#endif            
! 2012_12_18 . scb
!#else
!            do ig=1,ng
!               do im=1,2          
!                  alp1(im,ig)=0.25d0*(-atleakl(im,ig)+atleaku(im,ig))
!                  alp2(im,ig)=1.d0/12.d0*(-atleakl(im,ig)+2.d0*atleakm(im,ig)-atleaku(im,ig))
!
!                  sbal(im,ig)=-atleakm(im,ig)
!                  smm1(im,ig)=-alp1(im,ig)/3.d0
!                  smm2(im,ig)=-alp2(im,ig)/5.d0
!               enddo
!            enddo
!#endif

            ss0d=0.d0
            ss1d=0.d0
            ss2d=0.d0
            do ig=1,ng
               ss0d=ss0d+xsnff(ig,ih,iz)*(hflxf2(1,ig,ih,iz)-2.*hflxf2(2,ig,ih,iz))
               ss1d=ss1d+xsnff(ig,ih,iz)*(zmom1(1,ig,ih,iz)-2.*zmom1(2,ig,ih,iz))
               ss2d=ss2d+xsnff(ig,ih,iz)*(zmom2(1,ig,ih,iz)-2.*zmom2(2,ig,ih,iz))
            enddo

            do iter=1,1
               do ig=1,ng
                  ss0=reigv*xschif(ig,ih,iz)*ss0d
                  ss1=reigv*xschif(ig,ih,iz)*ss1d
                  ss2=reigv*xschif(ig,ih,iz)*ss2d
                  do ig2=xssfs(ig,ih,iz),xssfe(ig,ih,iz)
                     ss0=ss0+xssf(ig2,ig,ih,iz)*(hflxf2d(1,ig2,ih,iz)-2.*hflxf2d(2,ig2,ih,iz))
                     ss1=ss1+xssf(ig2,ig,ih,iz)*(zmom1(1,ig2,ih,iz)-2.*zmom1(2,ig2,ih,iz))
                     ss2=ss2+xssf(ig2,ig,ih,iz)*(zmom2(1,ig2,ih,iz)-2.*zmom2(2,ig2,ih,iz))
                  enddo
                  q0(1)=ss0 +sbal(1,ig)
                  q1(1)=ss1 +smm1(1,ig)
                  q2(1)=ss2 +smm2(1,ig)

                  q0(2)=-2.d0/3.d0*ss0 +sbal(2,ig)
                  q1(2)=-2.d0/3.d0*ss1 +smm1(2,ig)           
                  q2(2)=-2.d0/3.d0*ss2 +smm2(2,ig)

#ifdef paral_dbg
                  iodbg=100*iz+10*ig
                  !write(iodbg,*) hzm, xsdf(ig,ih,iz), xsdf2(ig,ih,iz),xstf(ig,ih,iz),xstrf(ig,ih,iz)
                  !write(iodbg,*) cntz0i(:,ig), cntz1i(:,ig),q0,q1,q2
                  !write(iodbg,*) cntzo2(:,ig,1,ih,iz),cntzo2(:,ig,2,ih,iz)
                  !write(iodbg,*) hflxf2(:,ig,ih,iz),zmom1(:,ig,ih,iz),zmom2(:,ig,ih,iz)
                  !write(iodbg,*) phisz(:,ig,1,ih,iz),phisz(:,ig,2,ih,iz)
                  !write(iodbg,*) cntz1i(:,ig)
                  write(iodbg,*) ig, ih, iz, cntzo2                           
#endif                  
                  ! solve the problem by NEM
                  !cdum1=cntzo2(:,ig,1,ih,iz)
                  !cdum2=cntzo2(:,ig,2,ih,iz)
                  call onenemz_sp3(hzm, xsdf(ig,ih,iz),xsdf2(ig,ih,iz), &
                                 xstf(ig,ih,iz),xstrf(ig,ih,iz), &
                                 cntz0i(:,ig), cntz1i(:,ig), &
                                 q0,q1,q2, &
                                 cntzo2(:,ig,1,ih,iz),cntzo2(:,ig,2,ih,iz), &   ! 2014_03_23 . scb
                                 !cdum1,cdum2, &
                                 hflxf2(:,ig,ih,iz),zmom1(:,ig,ih,iz),zmom2(:,ig,ih,iz), &
                                 phisz(:,ig,1,ih,iz),phisz(:,ig,2,ih,iz))
                  
                  !cntzo2(:,ig,1,ih,iz)=cdum1
                  !cntzo2(:,ig,2,ih,iz)=cdum2
#ifdef paral_dbg
                  iodbg=100*iz+10*ig+1
                  !write(iodbg,*) hzm, xsdf(ig,ih,iz), xsdf2(ig,ih,iz),xstf(ig,ih,iz),xstrf(ig,ih,iz)
                  !write(iodbg,*) cntz0i(:,ig), cntz1i(:,ig),q0,q1,q2
                  !write(iodbg,*) cntzo2(:,ig,1,ih,iz),cntzo2(:,ig,2,ih,iz)
                  !write(iodbg,*) hflxf2(:,ig,ih,iz),zmom1(:,ig,ih,iz),zmom2(:,ig,ih,iz)
                  !write(iodbg,*) phisz(:,ig,1,ih,iz),phisz(:,ig,2,ih,iz)
                  !write(iodbg,*) cntz1i(:,ig) 
                  write(iodbg,*) ig, ih, iz, cntzo2
                  !stop
#endif        
               enddo ! ig
            enddo
         enddo ! ih
      enddo ! iz
!$OMP END PARALLEL DO      

      !stop

      return
   end subroutine


#else
!
!
!      logical :: iftran  ! 2012_12_05 . scb
!
!      integer :: iz, ih, ig, im, ig2, it
!      real :: hzm, hzl, hzu, atleakm(2,ng), atleakl(2,ng), atleaku(2,ng)
!      real :: dcm(2,ng), dcu(2,ng), dcl(2,ng)   ! 2012_12_18 . scb
!      real :: cntz0i(2,ng), cntz1i(2,ng) 
!      integer :: neigdn, neigup, nn
!
!      real :: ss0d, ss1d, ss2d
!      real :: ss0, ss1, ss2
!      real :: q0(2), q1(2), q2(2)
!      integer :: iter
!
!      real :: sum, reflratl
!      real :: dc(2), tl0(2,ng), tl1(2,ng), alp1(2,ng), alp2(2,ng)
!      real :: sbal(2,ng), smm1(2,ng), smm2(2,ng)
!
!      integer :: ifrom, ito, istp
!
!! average transverse leakage cal. for z-direction 
!      do ih=1,nassy
!        do iz=1,nz
!          do ig=1,ng  
!            do im=1,2
!              sum=0
!              do it=1,6
!                nn=neignd(it,ih)
!                if(nn.eq.0) then
!                  reflratl=(1.d0-2.d0*alxr)/(1.d0+2.d0*alxr)
!                  sum=sum+(1.d0-reflratl)*cnto2(im,ig,it,ih,iz)
!                else
!                  sum=sum+cnto2(im,ig,it,ih,iz)-cnto2(im,ig,neigjin(it,ih),nn,iz)
!                endif
!              enddo
!              atleak2(im,ig,ih,iz)=sum*2.d0*rt3/9.d0/hside 
!! 2012_12_05 . scb              
!	            if(iftran) then
!	              atleak2(im,ig,ih,iz)=atleak2(im,ig,ih,iz)-sefve2(im,ig,ih,iz)
!              endif                  
!! added end              
!              hflxf2d(im,ig,ih,iz)=hflxf2(im,ig,ih,iz)
!              zmom1d(im,ig,ih,iz)=zmom1(im,ig,ih,iz)
!              zmom2d(im,ig,ih,iz)=zmom2(im,ig,ih,iz)
!            enddo 
!          enddo   
!        enddo   
!      enddo  
!
!      do iz=1,nz
!         do ih=1,nassy
!            ! preparation for axial solver(cnti,srczn,pbflx)
!            hzm=hz(iz) 
!            hzu=hzm
!            hzl=hzm 
!            !if(nz.gt.1) then   ! 2012_12_18 . scb
!               neigdn=neigz(1,iz) 
!               neigup=neigz(2,iz)
!               do ig=1,ng
!                  do im=1,2
!                     atleakm(im,ig)=atleak2(im,ig,ih,iz) 
!                     atleakl(im,ig)=atleakm(im,ig)
!                     atleaku(im,ig)=atleakm(im,ig)
!                  enddo 
!               enddo
!               if(neigdn.ne.0) then 
!                  hzl=hz(neigdn) 
!                  do ig=1,ng
!                     do im=1,2
!                        cntz0i(im,ig)=cntzo2(im,ig,2,ih,neigdn)
!                        atleakl(im,ig)=atleak2(im,ig,ih,neigdn)
!                     enddo
!                  enddo                       
!               else
!                  do ig=1,ng
!                     do im=1,2
!                        cntz0i(im,ig)=reflratzb*cntzo2(im,ig,1,ih,iz)
!                        if(.not.ifbcrefzb) atleakl(im,ig)=-atleak2(im,ig,ih,iz)   ! 2012_12_14 . scb
!                        !if(.not.ifbcrefzb) atleakl(im,ig)=0.d0   ! 2012_12_14 . scb
!                     enddo
!                  enddo
!               endif
!               if(neigup.ne.0) then 
!                  hzu=hz(neigup) 
!                  do ig=1,ng
!                     do im=1,2
!                        cntz1i(im,ig)=cntzo2(im,ig,1,ih,neigup)
!                        atleaku(im,ig)=atleak2(im,ig,ih,neigup)
!                     enddo  
!                  enddo
!               else
!                  do ig=1,ng
!                     do im=1,2
!                        cntz1i(im,ig)=reflratzt*cntzo2(im,ig,2,ih,iz)
!                        if(.not.ifbcrefzt) atleaku(im,ig)=-atleak2(im,ig,ih,iz)    ! 2012_12_14 . scb
!                        !if(.not.ifbcrefzt) atleaku(im,ig)=0.d0   ! 2012_12_14 . scb
!                     enddo
!                  enddo
!               endif
!            !endif   ! 2012_12_18 . scb
!
!! Transverse Leakage and coeff's at interface 
!!#ifdef DEBUG      ! 2012_12_18 . scb
!            do ig=1,ng
!               dc(1)=xsdf(ig,ih,iz)
!               dc(2)=xsdf2(ig,ih,iz)
!
!               do im=1,2
!                  tl0(im,ig)=(atleakl(im,ig)/hzl+atleakm(im,ig)/hzm) &
!                        /(dc(im)/hzl+dc(im)/hzm)
!                  tl1(im,ig)=(atleaku(im,ig)/hzu+atleakm(im,ig)/hzm) &
!                        /(dc(im)/hzu+dc(im)/hzm)
!
!                  alp1(im,ig)=(tl1(im,ig)-tl0(im,ig))/2.
!                  alp2(im,ig)=atleakm(im,ig)/dc(im)-(tl1(im,ig)+tl0(im,ig))/2.
!
!                  sbal(im,ig)=-atleakm(im,ig)
!                  smm1(im,ig)=-alp1(im,ig)*dc(im)/3.
!                  smm2(im,ig)=-alp2(im,ig)*dc(im)/5.
!               enddo
!            enddo
!! 2012_12_18 . scb
!!#else
!!            do ig=1,ng
!!               do im=1,2          
!!                  alp1(im,ig)=0.25d0*(-atleakl(im,ig)+atleaku(im,ig))
!!                  alp2(im,ig)=1.d0/12.d0*(-atleakl(im,ig)+2.d0*atleakm(im,ig)-atleaku(im,ig))
!!
!!                  sbal(im,ig)=-atleakm(im,ig)
!!                  smm1(im,ig)=-alp1(im,ig)/3.d0
!!                  smm2(im,ig)=-alp2(im,ig)/5.d0
!!               enddo
!!            enddo
!!#endif
!
!            ss0d=0.d0
!            ss1d=0.d0
!            ss2d=0.d0
!            do ig=1,ng
!               ss0d=ss0d+xsnff(ig,ih,iz)*(hflxf2(1,ig,ih,iz)-2.*hflxf2(2,ig,ih,iz))
!               ss1d=ss1d+xsnff(ig,ih,iz)*(zmom1(1,ig,ih,iz)-2.*zmom1(2,ig,ih,iz))
!               ss2d=ss2d+xsnff(ig,ih,iz)*(zmom2(1,ig,ih,iz)-2.*zmom2(2,ig,ih,iz))
!            enddo
!
!            do iter=1,1
!               do ig=1,ng
!                  ss0=reigv*xschif(ig,ih,iz)*ss0d
!                  ss1=reigv*xschif(ig,ih,iz)*ss1d
!                  ss2=reigv*xschif(ig,ih,iz)*ss2d
!                  do ig2=xssfs(ig,ih,iz),xssfe(ig,ih,iz)
!                     ss0=ss0+xssf(ig2,ig,ih,iz)*(hflxf2d(1,ig2,ih,iz)-2.*hflxf2d(2,ig2,ih,iz))
!                     ss1=ss1+xssf(ig2,ig,ih,iz)*(zmom1(1,ig2,ih,iz)-2.*zmom1(2,ig2,ih,iz))
!                     ss2=ss2+xssf(ig2,ig,ih,iz)*(zmom2(1,ig2,ih,iz)-2.*zmom2(2,ig2,ih,iz))
!                  enddo
!                  q0(1)=ss0 +sbal(1,ig)
!                  q1(1)=ss1 +smm1(1,ig)
!                  q2(1)=ss2 +smm2(1,ig)
!
!                  q0(2)=-2.d0/3.d0*ss0 +sbal(2,ig)
!                  q1(2)=-2.d0/3.d0*ss1 +smm1(2,ig)           
!                  q2(2)=-2.d0/3.d0*ss2 +smm2(2,ig)
!
!                  ! solve the problem by NEM
!                  call onenemz_sp3(hzm, xsdf(ig,ih,iz),xsdf2(ig,ih,iz), &
!                                 xstf(ig,ih,iz),xstrf(ig,ih,iz), &
!                                 cntz0i(:,ig), cntz1i(:,ig), &
!                                 q0,q1,q2, &
!                                 cntzo2(:,ig,1,ih,iz),cntzo2(:,ig,2,ih,iz), &
!                                 hflxf2(:,ig,ih,iz),zmom1(:,ig,ih,iz),zmom2(:,ig,ih,iz), &
!                                 phisz(:,ig,1,ih,iz),phisz(:,ig,2,ih,iz))
!               enddo ! ig
!            enddo
!         enddo ! ih
!      enddo ! iz
!
!      return
!      end subroutine
#endif      