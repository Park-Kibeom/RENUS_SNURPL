    subroutine readmoxbchxs(ifile,icbase)
      use param

      include 'global.h'
      include 'files.h'
      include 'cards.h'      
      include 'cntl.h'
      include 'geom.h'
      include 'xsec.h'
      include 'moxbch.inc'
      
      integer, intent(in) :: ifile
      logical :: ifnumeric, ifsigkf0,ifdetass
      integer :: iblock

      icomp=icbase
      !print *, icomp
      pnts(:,:,icomp)=0
      iblock=0
      ncomppnt=0
      do while (TRUE)
        read(ifile,'(a512)',end=1000) oneline
        write(io8,'(a)') trim(oneline)
        if(probe.eq.BANG .or. probe.eq.AST .or. .not.ifnumeric(oneline)) cycle
        backspace(ifile)
        
        !print *, icomp
        select case(iblock)
          case(0) ! point conf. 
            read(ifile, *) (npnt(iprop,icomp), iprop=1,nprop)
            do iprop=1,nprop
              if(npnt(iprop,icomp).eq.0) cycle
              ncomppnt=max(npnt(IDM,icomp),1)*max(npnt(IPPM,icomp),1)*max(npnt(ITF,icomp),1)*max(npnt(ITM,icomp),1)
              do while(true)
                read(ifile,'(a512)',end=1000) oneline            
                if(probe.ne.BANG .and. probe.ne.AST) exit
              enddo  

              read(oneline, *)(pnts(ipnt,iprop,icomp),ipnt=1,npnt(iprop,icomp))
! fuel temp. uses sqrt interpolation
              if(iprop.eq.3) pnts(:,iprop,icomp)=sqrt(pnts(:,iprop,icomp))              
            enddo
            iblock=iblock+1
          case(1) ! xsctr
            do m=1,ng
              do while(true)
                read(ifile,'(a512)',end=1000) oneline            
                if(probe.ne.BANG .and. probe.ne.AST) exit
              enddo  
              backspace(ifile)
              if(ncomppnt.eq.3) then
                read(ifile,*) (sigbchtr(3*i-2,m,icomp), i=1,3)
                do i=1,3
                  sigbchtr(3*i-1,m,icomp)=sigbchtr(3*i-2,m,icomp)
                  sigbchtr(3*i,m,icomp)=sigbchtr(3*i-1,m,icomp)
                  
                  sigbchtr(3*i-2+9:3*i+9,m,icomp)=sigbchtr(3*i-2:3*i,m,icomp)
                  sigbchtr(3*i-2+18:3*i+18,m,icomp)=sigbchtr(3*i-2:3*i,m,icomp)
                enddo
              else
                read(ifile,*) sigbchtr(1:ncomppnt,m,icomp)
              endif
            enddo
            iblock=iblock+1
          case(2) ! xsca
            do m=1,ng
              do while(true)
                read(ifile,'(a512)',end=1000) oneline            
                if(probe.ne.BANG .and. probe.ne.AST) exit
              enddo  
              backspace(ifile)
              if(ncomppnt.eq.3) then
                read(ifile,*) (sigbcha(3*i-2,m,icomp), i=1,3)
                do i=1,3
                  sigbcha(3*i-1,m,icomp)=sigbcha(3*i-2,m,icomp)
                  sigbcha(3*i,m,icomp)=sigbcha(3*i-1,m,icomp)

                  sigbcha(3*i-2+9:3*i+9,m,icomp)=sigbcha(3*i-2:3*i,m,icomp)
                  sigbcha(3*i-2+18:3*i+18,m,icomp)=sigbcha(3*i-2:3*i,m,icomp)
                enddo
              else
                read(ifile,*) sigbcha(1:ncomppnt,m,icomp)                  
              endif
            enddo
            iblock=iblock+1
          case(3) ! xscnf
            do m=1,ng
              do while(true)
                read(ifile,'(a512)',end=1000) oneline            
                if(probe.ne.BANG .and. probe.ne.AST) exit
              enddo  
              backspace(ifile)
              if(ncomppnt.eq.3) then
                read(ifile,*) (sigbchnf(3*i-2,m,icomp), i=1,3)
                do i=1,3
                  sigbchnf(3*i-1,m,icomp)=sigbchnf(3*i-2,m,icomp)
                  sigbchnf(3*i,m,icomp)=sigbchnf(3*i-1,m,icomp)

                  sigbchnf(3*i-2+9:3*i+9,m,icomp)=sigbchnf(3*i-2:3*i,m,icomp)
                  sigbchnf(3*i-2+18:3*i+18,m,icomp)=sigbchnf(3*i-2:3*i,m,icomp)
                enddo
              else
                read(ifile,*) sigbchnf(1:ncomppnt,m,icomp)
              endif
            enddo
            iblock=iblock+1
          case(4) ! xsckf
            do m=1,ng
              do while(true)
                read(ifile,'(a512)',end=1000) oneline            
                if(probe.ne.BANG .and. probe.ne.AST) exit
              enddo  
              backspace(ifile)
              if(ncomppnt.eq.3) then
                read(ifile,*) (sigbchkf(3*i-2,m,icomp), i=1,3)
                do i=1,3
                  sigbchkf(3*i-1,m,icomp)=sigbchkf(3*i-2,m,icomp)
                  sigbchkf(3*i,m,icomp)=sigbchkf(3*i-1,m,icomp)

                  sigbchkf(3*i-2+9:3*i+9,m,icomp)=sigbchkf(3*i-2:3*i,m,icomp)
                  sigbchkf(3*i-2+18:3*i+18,m,icomp)=sigbchkf(3*i-2:3*i,m,icomp)
                enddo
              else
                read(ifile,*) sigbchkf(:,m,icomp)
              endif
            enddo
            iblock=iblock+1
          case(5) ! sigbchs
            do ms=1,ng
              do md=1,ng
                do while(true)
                  read(ifile,'(a512)',end=1000) oneline
                  if(probe.ne.BANG .and. probe.ne.AST) exit
                enddo  
                backspace(ifile)
              if(ncomppnt.eq.3) then
                read(ifile,*) (sigbchs(3*i-2,ms,md,icomp), i=1,3)
                do i=1,3
                  sigbchs(3*i-1,ms,md,icomp)=sigbchs(3*i-2,ms,md,icomp)
                  sigbchs(3*i,ms,md,icomp)=sigbchs(3*i-1,ms,md,icomp)

                  sigbchs(3*i-2+9:3*i+9,ms,md,icomp)=sigbchs(3*i-2:3*i,ms,md,icomp)
                  sigbchs(3*i-2+18:3*i+18,ms,md,icomp)=sigbchs(3*i-2:3*i,ms,md,icomp)
                enddo
              else                
                read(ifile,*) sigbchs(:,ms,md,icomp)
              endif
              enddo
            enddo
            iblock=iblock+1
          case(6) ! adf
            do m=1,ng
              do while(true)
                read(ifile,'(a512)',end=1000) oneline            
                if(probe.ne.BANG .and. probe.ne.AST) exit
              enddo  
              backspace(ifile)
              if(ncomppnt.eq.3) then
                read(ifile,*) (sigbchadf(3*i-2,m,icomp), i=1,3)
                do i=1,3
                  sigbchadf(3*i-1,m,icomp)=sigbchadf(3*i-2,m,icomp)
                  sigbchadf(3*i,m,icomp)=sigbchadf(3*i-1,m,icomp)

                  sigbchadf(3*i-2+9:3*i+9,m,icomp)=sigbchadf(3*i-2:3*i,m,icomp)
                  sigbchadf(3*i-2+18:3*i+18,m,icomp)=sigbchadf(3*i-2:3*i,m,icomp)
                enddo
              else
                read(ifile,*) sigbchadf(:,m,icomp)
              endif
            enddo

            iblock=iblock+1
          case(7) ! fission spectrum
            do while(true)
              read(ifile,'(a512)',end=1000) oneline            
              if(probe.ne.BANG .and. probe.ne.AST) exit
            enddo  
            read(oneline,*) sigchi(:,icomp)
            sigchid(:,icomp)=sigchi(:,icomp)
            iblock=iblock+1
          case(8) ! inverse velocity
            do while(true)
              read(ifile,'(a512)',end=1000) oneline            
              if(probe.ne.BANG .and. probe.ne.AST) exit
            enddo  
            read(oneline,*) velof(:,icomp)
            velof(:,icomp)=1/velof(:,icomp)
            iblock=iblock+1            
          case(9) ! delay neutron decay constant
            do while(true)
              read(ifile,'(a512)',end=1000) oneline            
              if(probe.ne.BANG .and. probe.ne.AST) exit
            enddo  
            read(oneline,*) lmbdk(:,icomp)
            iblock=iblock+1
          case(10) ! delay neutron fraction
            do while(true)
              read(ifile,'(a512)',end=1000) oneline            
              if(probe.ne.BANG .and. probe.ne.AST) exit
            enddo  
            read(oneline,*) betak(:,icomp)
            iblock=iblock+1
          case(11)
            icomp=icomp+1
            if(icomp.le.ncomp) then
              npnt(:,icomp)=npnt(:,icomp-1)
              pnts(:,:,icomp)=pnts(:,:,icomp-1)
              iffuelc(icomp) = iffuelc(icomp-1)
            endif
            iblock=1
        end select
      enddo
!
1000  continue

      return
    end subroutine