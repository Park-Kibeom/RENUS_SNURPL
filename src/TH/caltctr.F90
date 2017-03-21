    subroutine caltctr(kth,lchan,rhouin,rhohuin,rhouinn,rhohuinn,qc,dzcetadt)
! transient coolant temperature calculation
      use param
      use allocs

      include 'global.h'
      include 'thgeom.inc'
      include 'thfuel.inc'
      include 'thcntl.inc'
      include 'thfdbk.inc'
      include 'thcool.inc'
      
      data eps/0.001/     
!
      h=hcool(kth,lchan)
      rho=dcool(kth,lchan)
      rhouout=rhou(kth,lchan)
      rhohuout=rhohu(kth,lchan)
      uoutd=ud(kth,lchan)
      uout=u(kth,lchan)
      dz=hzth(kth)

! calculate derivatives
      hr=(1+eps)*h
      hl=(1-eps)*h
      dh=hr-hl
      rhor=fdensh(hr)
      rhol=fdensh(hl)
      drhodh=(rhor-rhol)/dh
      drhohdh=(rhor*hr-rhol*hl)/dh

! calculate coefficients
      delrhou=rhouout-rhouin
      delrhohu=rhohuout-rhohuin
      hinn=rhohuinn/rhouinn
      alpha=-dzcetadt*drhodh
      beta=rhouinn-rthetac*delrhou
      tmpterm=2*h-hinn
      gamma=alpha*tmpterm+2*beta
      delta=beta*tmpterm
      a=2*alpha
      b=gamma+dzcetadt*drhohdh
      c=rthetac*delrhohu+delta-rhohuinn-dz*qc
      sqterm=sqrt(b*b-4*a*c)
      delh1=0.5*(-b+sqterm)/a
      delh2=0.5*(-b-sqterm)/a

      hout=rhohuout/rhouout
      biga=dzcetadt*drhohdh+alpha*hout
      bigb=dz*qc-rhouinn*hout+rhohuinn+rthetac*(delrhou*hout-delrhohu)
      delh3=bigb/biga
      if(abs(delh1-delh3).lt.abs(delh2-delh3)) then
         delh=delh1
      else
         delh=delh2
      endif
      hn=h+delh
     
      rhon=fdensh(hn)
      tn=ftemp(hn)
      rhoh=rho*h
      rhohn=rhon*hn
      rhououtn=rhouinn-dzcetadt*(rhon-rho)-rthetac*delrhou
      rhohuoutn=rhohuinn-dzcetadt*(rhohn-rhoh)-rthetac*delrhohu+dz*qc
      houtn=rhohuoutn/rhououtn
      toutn=ftemp(houtn)
      rhooutn=fdensh(houtn)
      uoutn=rhououtn/rhooutn
      uavg=thetac*(uoutn+rthetac*uout)
      if(((uoutd-uout)*(uoutn-uout).lt.-1.e-8*uavg*uavg) .and. (thetac.lt.1.0)) then
         uoutn=uavg
         rhououtn=rhooutn*uoutn
         rhohuoutn=rhououtn*houtn
      endif

! update variables
      tcool(kth,lchan)=tn
      dcool(kth,lchan)=rhon
      hcool(kth,lchan)=hn
      ud(kth,lchan)=u(kth,lchan)
      u(kth,lchan)=uoutn
      rhouin=rhou(kth,lchan)
      rhohuin=rhohu(kth,lchan)
      rhou(kth,lchan)=rhououtn
      rhohu(kth,lchan)=rhohuoutn
      rhouinn=rhououtn
      rhohuinn=rhohuoutn

      return
      
    end subroutine
