    subroutine drivethtr(deltm)
! calculates t/h in transient state
      use param
      use timer

      include 'global.h'
      include 'geom.h'
      include 'times.h'
      include 'thgeom.inc'
      include 'thcntl.inc'
      include 'thfuel.inc'
      include 'thcool.inc'
      include 'thfdbk.inc'
      include 'frzfdbk.inc'
      include 'thop.inc'
      include 'pow.h'
! added in ARTOS ver. 0.2 ( for Thermal Expansion ). 2012_07_03 by SCB
      include 'thexpan.inc'
! added end
      
      call timeron()
      
      tfmax=0.
      tmmax=0.   ! 2014_10_06 . scb
      
      dr2odt=delr2/deltm
      tw2odt=tw2/deltm

      if(relp(1,1).eq.0) then
         kbeg=2
      else
         kbeg=1
      endif
      if(relp(nzth,1).eq.0) then
         kend=nzth-1
      else
         kend=nzth
      endif
      
      toutavg=0  
      do lchan=1,nchan
        kth=0
        rhouin=rhou(kth,lchan)
        rhohuin=rhohu(kth,lchan)
        rhouinn=rhouin
        rhohuinn=rhohuin
        do kth=1,nzth
          dz=hzth(kth)
          qprime=plevel*powlin*relp(kth,lchan)
          qfn=fracdf*qprime/afp
          qf=thetafb*qvol(kth,lchan)+thetaf*qfn
          qvol(kth,lchan)=qfn
          htcoef(kth,lchan)=fhtcoef(tcool(kth,lchan),deq,rhou(kth,lchan))
          call caltftr(kth,lchan,tcool(kth,lchan),htcoef(kth,lchan),qf,dr2odt,tw2odt,TRUE)
          tfmax=max(tfmax,tfuel(1,kth,lchan))
          tmmax=max(tmmax,tcool(kth,lchan))   ! 2014_10_06 . scb
          qflux=htcoef(kth,lchan)*(tfuel(nrp4,kth,lchan)-tcool(kth,lchan))
          dzcetadt=dz/(thetac*deltm)
          qc=fracdc*qprime/acf
          qeffn=qflux*zetap+qc
          qc=qeffn+rthetac*qeff(kth,lchan)
          qeff(kth,lchan)=qeffn
          call caltctr(kth,lchan,rhouin,rhohuin,rhouinn,rhohuinn,qc,dzcetadt)
        enddo
        hout=rhohu(kend,lchan)/rhou(kend,lchan)
        tout=ftemp(hout)
        toutavg=toutavg+tout
      enddo
      toutavg=toutavg/nchan

! added in ARTOS ver. 0.2 ( for Thermal Expansion ). 2012_07_03 by SCB
! average temperature
      tfuelavg=0.d0
      tcoolavg=0.d0
      dcoolavg=0.d0    ! 2014_09_05 . scb added dcoolavg
      do l=1,nchan
        do k=kfsth,kfeth
          tfuelavg=tfuelavg+tdopl(k,l)**2*hzth(k)
          tcoolavg=tcoolavg+tcool(k,l)*hzth(k)
          dcoolavg=dcoolavg+dcool(k,l)*hzth(k)    ! 2014_09_05 . scb added dcoolavg
        enddo
      enddo
      tfuelavg=tfuelavg/(nchan*0.01*coreheight) 
      !tcoolavg=tcoolavg/(nchan*0.01*coreheight)+CKELVIN 
      tcoolavg=tcoolavg/(nchan*0.01*coreheight)   ! 2014_10_06 . scb
      dcoolavg=dcoolavg/(nchan*0.01*coreheight)    ! 2014_09_05 . scb added dcoolavg
! added end      
	
      call timeroff(tth)

      write(mesg,'(" Max. Tf=",f9.4," C,  Avg. Outlet Temp.=",f9.4," C")') tfmax,toutavg
      call message(true,true,mesg)
      
      return
      
    end subroutine