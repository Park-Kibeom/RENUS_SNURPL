! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
   subroutine swpnemz_sp3(iftran)
   !subroutine swpnemz_sp3_new(iftran)

      use const
      use geomhex
      use tpen_sp3
      use cmfdhex2g, only : phisz
      use xsec
      use sfam, only : reigv, phi_p2
      implicit none

!#define paral_dbg   ! 2014_03_23 . scb

      logical :: iftran  ! 2012_12_05 . scb

      integer :: iz, ih, ig, im, ig2, it
      real :: hzm, hzl, hzu, atleakm(2,ng), atleakl(2,ng), atleaku(2,ng)
      !real :: cntz0i(2,ng), cntz1i(2,ng) 
      real :: cntz0i(2,ng,nassy,nz), cntz1i(2,ng,nassy,nz)  ! 2014_04_01 . scb
      integer :: neigdn, neigup, nn

      real :: ss0d, ss1d, ss2d
      real :: ss0, ss1, ss2
      real :: q0(2), q1(2), q2(2)
      integer :: iter

      real :: sum, reflratl
      real :: alp1(2,ng), alp2(2,ng)
      real :: sbal(2,ng,nassy,nz), smm1(2,ng,nassy,nz), smm2(2,ng,nassy,nz)

      integer :: ifrom, ito, istp
! 2013_07_10 . scb      
      real :: ratl, ratr, difsl(2,ng), difsr(2,ng)
! added end

! 2014_03_23 . scb for paral dbg
      logical,save :: first=.true.
      integer :: iodbg
      
      real :: cdum1(2), cdum2(2)
! added end      
      integer :: ig1,ih1,iz1

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

!!$OMP PARALLEL DO  private(iz,ih,hzm,hzu,hzl,neigdn,neigup,ig,im,atleakm,atleakl,atleaku,alp1,alp2,ratl,ratr,difsl,difsr,iodbg)
      do iz=1,nz
         do ih=1,nassy
            ! preparation for axial solver(cnti,srczn,pbflx)
            hzm=hz(iz) 
            hzu=hzm
            hzl=hzm 
            
            neigdn=neigz(1,iz) 
            neigup=neigz(2,iz)

            atleakm=atleak2(:,:,ih,iz) 
            atleakl=atleakm
            atleaku=atleakm

            if(neigdn.ne.0) then 
              hzl=hz(neigdn) 
              atleakl=atleak2(:,:,ih,neigdn)               
            else
              if(.not.ifbcrefzb) then
                atleakl=0.d0   ! 2012_12_14 . scb       
                hzl=0.d0       ! 2013_07_10 . scb  
              endif              
            endif
            if(neigup.ne.0) then 
              hzu=hz(neigup) 
              atleaku=atleak2(:,:,ih,neigup)
            else
              if(.not.ifbcrefzt) then
                atleaku=0.d0   ! 2012_12_14 . scb
                hzu=0.d0       ! 2013_07_10 . scb  
              endif              
            endif

! Transverse Leakage and coeff's at interface 
!#ifdef DEBUG      ! 2012_12_18 . scb
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
                  sbal(im,ig,ih,iz)=-atleakm(im,ig)
                  smm1(im,ig,ih,iz)=-alp1(im,ig)/3.d0
                  smm2(im,ig,ih,iz)=-alp2(im,ig)/5.d0
               enddo
            enddo
         enddo
      enddo
!!$OMP END PARALLEL DO      
      
      cntz0i=0.d0
      cntz1i=0.d0
!!$OMP PARALLEL DO  private(iz,ih,ss0d,ss1d,ss2d,ig,ss0,ss1,ss2,ig2,q0,q1,q2,iodbg,cntz0i,cntz1i,ig1,ih1,iz1)
      do iz=1,nz
         do ih=1,nassy      
            if(neigz(1,iz).ne.0) then 
              cntz0i(:,:,ih,iz)=cntzo2(:,:,2,ih,neigz(1,iz))    
!              write(401,*)  cntzo2(:,:,2,ih,neigz(1,iz))    
            else
              cntz0i(:,:,ih,iz)=reflratzb*cntzo2(:,:,1,ih,iz)
!              write(401,*)  reflratzb*cntzo2(:,:,1,ih,iz)
            endif
            
            if(neigz(2,iz).ne.0) then 
              cntz1i(:,:,ih,iz)=cntzo2(:,:,1,ih,neigz(2,iz))
!              write(401,*)  cntzo2(:,:,1,ih,neigz(2,iz))
            else
              cntz1i(:,:,ih,iz)=reflratzt*cntzo2(:,:,2,ih,iz)
!              write(401,*)  cntzo2(:,:,1,ih,neigz(2,iz))
            endif
           
            ss0d=0.d0
            ss1d=0.d0
            ss2d=0.d0
            do ig=1,ng
               ss0d=ss0d+xsnff(ig,ih,iz)*(hflxf2(1,ig,ih,iz)-2.*hflxf2(2,ig,ih,iz))
               ss1d=ss1d+xsnff(ig,ih,iz)*(zmom1(1,ig,ih,iz)-2.*zmom1(2,ig,ih,iz))
               ss2d=ss2d+xsnff(ig,ih,iz)*(zmom2(1,ig,ih,iz)-2.*zmom2(2,ig,ih,iz))
            enddo

            do ig=1,ng
              ss0=reigv*xschif(ig,ih,iz)*ss0d
              ss1=reigv*xschif(ig,ih,iz)*ss1d
              ss2=reigv*xschif(ig,ih,iz)*ss2d
              do ig2=xssfs(ig,ih,iz),xssfe(ig,ih,iz)
                  ss0=ss0+xssf(ig2,ig,ih,iz)*(hflxf2d(1,ig2,ih,iz)-2.*hflxf2d(2,ig2,ih,iz))
                  ss1=ss1+xssf(ig2,ig,ih,iz)*(zmom1(1,ig2,ih,iz)-2.*zmom1(2,ig2,ih,iz))
                  ss2=ss2+xssf(ig2,ig,ih,iz)*(zmom2(1,ig2,ih,iz)-2.*zmom2(2,ig2,ih,iz))
              enddo
              q0(1)=ss0 +sbal(1,ig,ih,iz)
              q1(1)=ss1 +smm1(1,ig,ih,iz)
              q2(1)=ss2 +smm2(1,ig,ih,iz)

              q0(2)=-2.d0/3.d0*ss0 +sbal(2,ig,ih,iz)
              q1(2)=-2.d0/3.d0*ss1 +smm1(2,ig,ih,iz)           
              q2(2)=-2.d0/3.d0*ss2 +smm2(2,ig,ih,iz)

#ifdef paral_dbg
              !do ig1=1,ng
              !  do ih1=1,nassy
              !    do iz1=1,nz
              !      write(401,*) ig1,ih1,iz1
              !      write(401,*) 'cntz0i',cntz0i(:,ig1,ih1,iz1)
              !      write(401,*) 'cntz1i',cntz1i(:,ig1,ih1,iz1)
              !      write(401,*) 'cntzo21',cntzo2(:,ig1,1,ih1,iz1)
              !      write(401,*) 'cntzo22',cntzo2(:,ig1,2,ih1,iz1)       
              !    enddo
              !  enddo
              !enddo
              !write(1401,*) 'cntz0i',cntz0i
              !write(1401,*) 'cntz1i',cntz1i
              !write(1401,*) 'cntzo21',cntzo2
              !write(1401,*) 'cntzo22',cntzo2        
              !iodbg=401
              !write(iodbg,*) 1,2,3
              !!stop
#endif                  
              ! solve the problem by NEM
              call onenemz_sp3(hz(iz) , xsdf(ig,ih,iz),xsdf2(ig,ih,iz), &
                              !xstf(ig,ih,iz),xstrf(ig,ih,iz), &
                              xstf(ig,ih,iz),xstfsp3(ig,ih,iz), &   ! 2014_10_07 . scb
                              cntz0i(:,ig,ih,iz), cntz1i(:,ig,ih,iz), &
                              q0,q1,q2, &
                              cntzo2(:,ig,1,ih,iz),cntzo2(:,ig,2,ih,iz), &   ! 2014_03_23 . scb
                              !cdum1,cdum2, &
                              hflxf2(:,ig,ih,iz),zmom1(:,ig,ih,iz),zmom2(:,ig,ih,iz), &
                              phisz(:,ig,1,ih,iz),phisz(:,ig,2,ih,iz))
              
! 2015_01_13 . scb for moment shape
!#define sol_shape              
#ifdef sol_shape
              !write(113,*) 'ig, ih, iz', ig, ih, iz
              !write(113,*) 'hz(iz)', hz(iz)
              !write(113,*) 'hflxf2(:,ig,ih,iz)', hflxf2(:,ig,ih,iz)
              !write(113,*) 'zmom1(:,ig,ih,iz)', zmom1(:,ig,ih,iz)
              !write(113,*) 'zmom2(:,ig,ih,iz)', zmom2(:,ig,ih,iz)
              !write(113,*) 'phisl', phisz(:,ig,1,ih,iz)
              !write(113,*) 'phisr', phisz(:,ig,2,ih,iz)
              write(113,*) ig, ih, iz
              write(113,*) hz(iz)
              write(113,*) hflxf2(1,ig,ih,iz)-2.*hflxf2(2,ig,ih,iz), hflxf2(2,ig,ih,iz)
              write(113,*) zmom1(1,ig,ih,iz)-2*zmom1(2,ig,ih,iz), zmom1(2,ig,ih,iz)
              write(113,*) zmom2(1,ig,ih,iz)-2*zmom2(2,ig,ih,iz), zmom2(2,ig,ih,iz)
              write(113,*) phisz(1,ig,1,ih,iz)-2*phisz(2,ig,1,ih,iz), phisz(2,ig,1,ih,iz)
              write(113,*) phisz(1,ig,2,ih,iz)-2*phisz(2,ig,2,ih,iz), phisz(2,ig,2,ih,iz)
              write(113,*) xsnff(ig,ih,iz), xsdf(ig,ih,iz)
              
              write(115,*) ig, ih, iz
              write(115,*) 'jinl',cntz0i(:,ig,ih,iz)
              write(115,*) 'jinr', cntz1i(:,ig,ih,iz)
              write(115,*) 'joutl', cntzo2(:,ig,1,ih,iz)
              write(115,*) 'joutr', cntzo2(:,ig,2,ih,iz)
              write(115,*) 'jnetl', cntz0i(:,ig,ih,iz)-cntzo2(:,ig,1,ih,iz)
              write(115,*) 'jnetr', cntzo2(:,ig,2,ih,iz)-cntz1i(:,ig,ih,iz)
#endif  
! added end
        !pause
#ifdef paral_dbg
              !do ig1=1,ng
              !  do ih1=1,nassy
              !    do iz1=1,nz
              !      write(401,*) ig1,ih1,iz1
              !      !write(401,*) 'cntz0i',cntz0i(:,ig1,ih1,iz1)
              !      !write(401,*) 'cntz1i',cntz1i(:,ig1,ih1,iz1)
              !      write(401,*) 'cntzo21',cntzo2(:,ig1,1,ih1,iz1)
              !      write(401,*) 'cntzo22',cntzo2(:,ig1,2,ih1,iz1)       
              !    enddo
              !  enddo
              !enddo
#endif        
            enddo ! ig
         enddo ! ih
      enddo ! iz
!!$OMP END PARALLEL DO      
       
#ifdef sol_shape
      close(113)
      close(115)
#endif            

#ifdef paral_dbg
      stop
#endif      

      return
   end subroutine

