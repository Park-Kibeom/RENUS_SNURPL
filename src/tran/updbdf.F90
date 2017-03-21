! 2014_09_18 . scb : use initial value as previous solutions...
    !subroutine updbdf(istep, bdforder, autodt, deltm, ng, nxy, nz)
    subroutine updbdf(istep, bdforder, autodt, deltm, ng, nxy, nz,ifsp3) ! 2013_08_12 . scb
    
      use sfam,     only : phi, phif, psi
      use tran,     only : psit, psitd, betap, betapd
      use bdf
      use geomhex,  only : ifhexsp3       ! 2012_12_06 . scb
      use tpen_sp3, only : hflx2, hflxf2  ! 2012_12_06 . scb
      use sp3senm,  only : phishp       ! 2013_08_12 . scb
      use const             ! 2014_09_23 . scb
      use sfam_cntl,  only : ifrect   ! 2014_09_26 . scb
      use sp3senm,  only : phi2    ! 2014_11_13 . scb
      
      integer :: istep, bdforder, ng, nxy, nz
      logical :: autodt
      real*8 :: deltm
      
      !integer,parameter :: ng2=2
      
      logical :: ifsp3 ! 2013_08_12 . scb
      
      bdforder=bdfordersave
      
      if(autodt) then
        automain=.TRUE.
        if(deltmupd) then
          deltm=deltmsave
          deltmupd=.FALSE.
        endif
      endif
      
      if(updneed) then              
        do i=bdforder,2,(-1)
          do iz=1,nz
            do ixy=1,nxy
              do ig=1,ng
                phifbdf(ig,ixy,iz,i)=phifbdf(ig,ixy,iz,i-1)
              enddo
            enddo
          enddo          
        enddo         
        do iz=1,nz
          do ixy=1,nxy
            do ig=1,ng
              phifbdf(ig,ixy,iz,1)=phif(ig,ixy,iz)
            enddo
          enddo
        enddo        
           
        if(ng.ne.ng2) then
          do i=bdforder,2,(-1)
            do iz=1,nz
              do ixy=1,nxy
                do ig=1,ng2
                  phibdf(ig,ixy,iz,i)=phibdf(ig,ixy,iz,i-1)
                enddo
              enddo
            enddo          
            !phifbdf(:,:,:,i)=phifbdf(:,:,:,i-1)
          enddo
          do iz=1,nz
            do ixy=1,nxy
              do ig=1,ng2
                phibdf(ig,ixy,iz,1)=phi(ig,ixy,iz)
              enddo
            enddo
          enddo       
          !phifbdf(:,:,:,1)=phif(:,:,:)   
        endif         

! 2012_12_06 . scb        
        if(ifsp3) then
          do i=bdforder,2,(-1)
            do iz=1,nz
              do ixy=1,nxy
                do ig=1,ng
                  phifbdf2(ig,ixy,iz,i)=phifbdf2(ig,ixy,iz,i-1)
                enddo
              enddo
            enddo        
            !phibdf2(:,:,:,i)=phibdf2(:,:,:,i-1)
          enddo                   
          
          if(ifhexsp3) then
            do iz=1,nz
              do ixy=1,nxy
                do ig=1,ng
                  phifbdf2(ig,ixy,iz,1)=hflxf2(2,ig,ixy,iz)         
                enddo
              enddo
            enddo
          else
! 2014_11_13 . scb            
!! 2013_08_13 . scb            
            do iz=1,nz
              do ixy=1,nxy
                do ig=1,ng
                  phifbdf2(ig,ixy,iz,1)=0.d0
                  do idir=1,3
                    phifbdf2(ig,ixy,iz,1)=phifbdf2(ig,ixy,iz,1)+phishp(2,0,idir,ixy,iz,ig)      
                  enddo
                  phifbdf2(ig,ixy,iz,1)=phifbdf2(ig,ixy,iz,1)/3.d0
                enddo
              enddo
            enddo
!! added end    
            !do iz=1,nz
            !  do ixy=1,nxy
            !    do ig=1,ng
            !      phifbdf2(ig,ixy,iz,1)=phi2(ixy,iz,ig)
            !    enddo
            !  enddo
            !enddo
! added end            
          endif
                               
          if(ng.ne.ng2) then
            do i=bdforder,2,(-1)
              do iz=1,nz
                do ixy=1,nxy
                  do ig=1,ng2
                    phibdf2(ig,ixy,iz,i)=phibdf2(ig,ixy,iz,i-1)
                  enddo
                enddo
              enddo   
              !phifbdf2(:,:,:,i)=phifbdf2(:,:,:,i-1)
            enddo        
            
            if(ifhexsp3) then
              do iz=1,nz
                do ixy=1,nxy
                  do ig=1,ng2
                    phibdf2(ig,ixy,iz,1)=hflx2(2,ig,ixy,iz)         
                  enddo
                enddo
              enddo
            endif                        
          endif
        endif
! added end        
        
        if(automain) then
          do i=1,nz
            do j=1,nxy
              psitd(j,i)=psit(j,i)
              psit(j,i)=psi(j,i)
              betapd(j,i)=betap(j,i)
            enddo
          enddo        
        endif

        do i=bdforder,2,(-1)
          deltmarray(i)=deltmarray(i-1)
        enddo
        
! 2014_09_22 . scb        
        !if(ifrect .and. .not. ifsp3) then
        if(flagsrcdt) then
          do iorder=bdforder,1,(-1)  
            iorder1=iorder-1
            do k=1,nz
              do l=1,nxy
                do m=1,ng
                  srcdtbdf(m,l,k,XDIR,iorder)=srcdtbdf(m,l,k,XDIR,iorder1)
                  srcdtbdf(m,l,k,YDIR,iorder)=srcdtbdf(m,l,k,YDIR,iorder1)
                  srcdtbdf(m,l,k,ZDIR,iorder)=srcdtbdf(m,l,k,ZDIR,iorder1)
                enddo
              enddo
            enddo 
          enddo            
! 2014_09_26 . scb            
        elseif(flagsrcdt2) then
! 2014_10_06 . scb            
          if(ifsp3) then                        
            do iorder=bdforder,1,(-1)  
              iorder1=iorder-1
              do k=1,nz
                do l=1,nxy
                  do m=1,ng
                    do im=1,2
                      dtcff2(:,im,m,l,k,XDIR,iorder)=dtcff2(:,im,m,l,k,XDIR,iorder1)
                      dtcff2(:,im,m,l,k,YDIR,iorder)=dtcff2(:,im,m,l,k,YDIR,iorder1)
                      dtcff2(:,im,m,l,k,ZDIR,iorder)=dtcff2(:,im,m,l,k,ZDIR,iorder1)            
                    enddo
                  enddo
                enddo 
              enddo     
            enddo              
! added end                
          else              
            do iorder=bdforder,1,(-1)  
              iorder1=iorder-1
              do k=1,nz
                do l=1,nxy
                  do m=1,ng
                    dtcff(:,m,l,k,XDIR,iorder)=dtcff(:,m,l,k,XDIR,iorder1)
                    dtcff(:,m,l,k,YDIR,iorder)=dtcff(:,m,l,k,YDIR,iorder1)
                    dtcff(:,m,l,k,ZDIR,iorder)=dtcff(:,m,l,k,ZDIR,iorder1)            
                  enddo
                enddo
              enddo 
            enddo     
          endif            
        endif
! added end          
        !endif        
! added end        
        
      else
        do i=1,nz
          do j=1,nxy
            psi(j,i)=psit(j,i)
            betap(j,i)=betapd(j,i)
          enddo
        enddo
        
        do iz=1,nz
          do ixy=1,nxy
            do ig=1,ng
              phif(ig,ixy,iz)=phifbdf(ig,ixy,iz,1)
            enddo
          enddo
        enddo  
        
        if(ng.ne.ng2) then
          do iz=1,nz
            do ixy=1,nxy
              do ig=1,ng2
                phi(ig,ixy,iz)=phibdf(ig,ixy,iz,1)
              enddo
            enddo
          enddo  
        endif        
          
        !phi(:,:,:)=phibdf(:,:,:,1)
        !if(ng.ne.ng2) phif(:,:,:)=phifbdf(:,:,:,1)
        
! 2012_12_06 . scb
        if(ifhexsp3) then
          do iz=1,nz
            do ixy=1,nxy
              do ig=1,ng
                !hflx2(ig,ixy,iz,1)=phibdf2(2,ig,ixy,iz)         
                hflxf2(2,ig,ixy,iz)=phifbdf2(ig,ixy,iz,1)   ! 2012_12_26 . scb (debugged)
              enddo
            enddo
          enddo
          
          if(ng.ne.ng2) then
            do iz=1,nz
              do ixy=1,nxy
                do ig=1,ng2
                  hflx2(2,ig,ixy,iz)=phibdf2(ig,ixy,iz,1)         
                enddo
              enddo
            enddo
          endif
        endif        
! added end
        
        updneed=.TRUE.
      endif
      deltmarray(1)=deltm
      
      call calbdfcoef(bdforder)
      !call calbdfcoef_const(bdforder)
      
      imaxiter=imaxiter+1
      
    end subroutine