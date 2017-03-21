!#define DLL_test  ! 2013_10_22 . scb
    subroutine drivethss
! calculates t/h in steady-state 
      use param
      use timer

      include 'global.h'
      include 'times.h'
      include 'thgeom.inc'
      include 'thcntl.inc'
      include 'thfuel.inc'
      include 'thcool.inc'
      include 'thfdbk.inc'
      include 'frzfdbk.inc'
      include 'thop.inc'      
      include 'pow.h'
      include 'xsec.h'
! added in ARTOS ver. 0.2 ( for Thermal Expansion, Core Height ). 2012_07_03 by SCB      
      include 'thexpan.inc' ! Thermal Expansion
      include 'geom.h' ! coreheight
! added end	
      
      integer, save :: iii=0
      character*79 :: filetmp
      
      real :: tlap,tfmaxt,toutavgt,fnchan
! 2013_10_22 . scb      
#ifdef DLL_test
      logical,save :: flagfinal=.false.
      
      common /dlltest / flagfinal
#endif      
! added end
      
      call timeron()
      ntth=ntth+1
      
      tfmax=0.d0
      tmmax=0.d0   ! 2014_10_06 . scb
      toutavg=0.d0
      dcoolavg=0.d0
      tcoolavg=0.d0
      tdoplmax=0.d0
      
      fac=1./(nchan*nzth)      
      if(freezedm) then
        do lchan=1,nchan
          do kth=1,nzth
            dcool(kth,lchan)=frozendm
            dcoolavg=frozendm
           enddo
         enddo
      else
! coolant temperature calculation
        do l=1,nchan
          k=0
          rhouin=rhou(k,l)
          rhohuin=rhohu(k,l)
          do k=1,nzth
            qprime=plevel*powlin*relp(k,l)
            qf=fracdf*qprime/afp
            qvol(k,l)=qf
            qc=fracdc*qprime/acf
            qflux=qf*afp/zeta
            qeffnew=qflux*zetap+qc
            qeff(k,l)=qeffnew
            
            rhohuout=rhohuin+qeffnew*hzth(k)            
            hcool(k,l)=0.5*(rhohuout+rhohuin)/rhouin
            tcool(k,l)=ftemp(hcool(k,l))
            dcool(k,l)=fdensh(hcool(k,l))
            rhohuin=rhohuout
            rhohu(k,l)=rhohuout
            dcoolavg=dcoolavg+dcool(k,l)
            tcoolavg=tcoolavg+tcool(k,l)
          enddo
          hout=rhohu(nzth,l)/rhou(nzth,l)
          tout=ftemp(hout)
          toutavg=toutavg+tout
        enddo      
        dcoolavg=dcoolavg*fac
        tcoolavg=tcoolavg*fac
        toutavg=toutavg/nchan
      endif
      
! fuel temperature calculation     
!      tfuelavg=0   ! commented in ARTOS ver. 0.2 . 2012_07_04 by SCB    
      if(freezetf) then
        do lchan=1,nchan
          do kth=1,nzth
            tfuel(:,kth,lchan)=frozentf
            tdopl(kth,lchan)=sqrt(frozentf+CKELVIN)
            tfuelavg=frozentf
           enddo
         enddo        
      else
        do l=1,nchan
          k=0
          do k=1,nzth
            qprime=plevel*powlin*relp(k,l)
            qf=fracdf*qprime/afp
! heat transfer coefficient
            htcoef(k,l)=fhtcoef(tcool(k,l),deq,rhou(k,l))
            call caltfss(k,l,tcool(k,l),htcoef(k,l),qf)
            tfmax=max(tfmax,tfuel(1,k,l))
            tmmax=max(tmmax,tcool(k,l))   ! 2014_10_06 . scb
!            tfuelavg=tfuelavg+tdopl(k,l)**2   ! commented in ARTOS ver. 0.2 . 2012_07_04 by SCB  
          enddo
        enddo
!        tfuelavg=tfuelavg*fac   ! commented in ARTOS ver. 0.2 . 2012_07_04 by SCB  
      endif
      
! added in ARTOS ver. 0.2 . 2012_07_04 by SCB    
! average temperature
      tfuelavg=0.d0
      tcoolavg=0.d0   
      dcoolavg=0.d0   ! 2014_09_05 . scb
      do l=1,nchan
        do k=kfsth,kfeth
          !tfuelavg=tfuelavg+tfuel(nrp5,k,l)*hzth(k)
		  tfuelavg=tfuelavg+tdopl(k,l)**2*hzth(k)
		  tcoolavg=tcoolavg+tcool(k,l)*hzth(k)
      dcoolavg=dcoolavg+dcool(k,l)*hzth(k)   ! 2014_09_05 . scb
        enddo
      enddo
      tfuelavg=tfuelavg/(nchan*0.01*coreheight) !+CKELVIN
      !tcoolavg=tcoolavg/(nchan*0.01*coreheight)+CKELVIN	
      tcoolavg=tcoolavg/(nchan*0.01*coreheight)	  ! 2014_10_06 . scb
      dcoolavg=dcoolavg/(nchan*0.01*coreheight)   ! 2014_09_05 . scb
! added end
      
      if(.not.freezetf .and. tdoplmax.lt.epstf) then
         flagth=true
      else
         flagth=FALSE
      endif
      
      write(mesg,601) tfmax,toutavg,tdoplmax,flagth
      call message(true,true,mesg)

! added in ARTOS ver. 0.2 . 2012_07_04 by SCB    
!      write(mesg,602) tfuelavg,tcoolavg+CKELVIN,dcoolavg
      write(mesg,602) tfuelavg,tcoolavg,dcoolavg
! added end
      call message(true,true,mesg)

  601 format(" Max. Tf=",f7.2,"C,  Avg. Outlet Temp.=",f7.2,"C", ", Max. Doppler Change=",1p,e9.2,l2)  
  !602 format(" Avg. Fuel. Temp.=",f7.2,"K,  Avg. Cool. Temp.=",f7.2,"K", ", Avg. Cool. Density=",f9.2)          
  602 format(" Avg. Fuel. Temp.=",f7.2,"K,  Avg. Cool. Temp.=",f7.2,"C", ", Avg. Cool. Density=",f9.2)          
     
      call timeroff(tth)    
      
! 2013_10_22 . scb      
#ifdef DLL_test
      open(unit=1022,file='DLL_TH',status='unknown')
      
      write(unit,*) nxy, nz
      
      do ixy=1,nxy
        la=ltola(l)                  
        iat=iassytyp(la)           
        if(iffuela(iat)) then      
          ichan=ichan+1           
        endif                     
      enddo                        
      nchan=ichan                  
      
#endif      
! added end      
      
    end subroutine
      