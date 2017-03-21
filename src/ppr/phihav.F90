    function phihav(xlower,xupper,ylower,yupper,l,k,m)

      use param
      use ppr, only : kappa, pc2d,hc2d, sqrt2
      
      include 'global.h'
      include 'geom.h'
      include 'xsec.h'
      include 'nodal.h'
      include 'pinpr.h'

      real::xupper,xlower,yupper,ylower
      real::tpc(0:14),thc(1:8),temp(1:8)
      real::vl(1:2),vr(1:2),dvl(1:2),dvr(1:2)
      real::kl,rkl,ml
      real::coshm(1:2),sinhm(1:2),dif(1:2)
      real::sinhkvl2,sinhkvr2,sinhkvl1,sinhkvr1
      logical::flag


      kl=kappa(l,k,m)
      rkl=1./kl
      ml=kl/sqrt2
      oe=1./8.
      ot=0.5
      of=0.25
            
      vl(1)=ylower
      vl(2)=xlower
      vr(1)=yupper
      vr(2)=xupper
      dvl(1)=ylower*ylower
      dvl(2)=xlower*xlower
      dvr(1)=yupper*yupper
      dvr(2)=xupper*xupper

      tpc(0)=pc2d(0,l,k,m)

      do i=1,2 !1:y, 2:x
        temp(1)=ot*(vl(i)+vr(i))
        temp(2)=ot*(-1.+dvl(i)+vl(i)*vr(i)+dvr(i))
        temp(3)=oe*(vl(i)+vr(i))*(-6.+5.*dvl(i)+5.*dvr(i))
        temp(4)=oe*(3.+7.*dvl(i)*dvr(i)-10.*vl(i)*vr(i)  &
               -10.*dvl(i)+7.*dvl(i)*dvl(i)+7.*dvl(i)*vl(i)*vr(i)  &
               -10.*dvr(i)+7.*dvr(i)*dvr(i)+7.*dvr(i)*vr(i)*vl(i))
        if(i==1) then
          nadd=0
        elseif(i==2) then
          nadd=4
        endif
        do j=1,4
          tpc(j+nadd)=temp(j)*pc2d(j+nadd,l,k,m)
        enddo
      enddo  

      do i=1,2 ! 1:y 2:x
        temp(i)=of*(vl(i)+vr(i))
        temp(i+2)=-1.+dvl(i)+vl(i)*vr(i)+dvr(i)
      enddo
      tpc(9)=4.*temp(1)*temp(2)*pc2d(9,l,k,m)
      tpc(10)=temp(2)*temp(3)*pc2d(10,l,k,m)
      tpc(11)=temp(1)*temp(4)*pc2d(11,l,k,m)
      tpc(12)=of*temp(3)*temp(4)*pc2d(12,l,k,m)  
      tpc(13)=temp(1)*temp(2)*(-6+5*(dvl(1)+dvr(1)))*pc2d(13,l,k,m)
      tpc(14)=temp(1)*temp(2)*(-6+5*(dvl(2)+dvr(2)))*pc2d(14,l,k,m)
      do i=1,2 ! 1:yl-yr 2:xl-xr
        dif(i)=1./(vl(i)-vr(i))
      enddo

      coshkvl2=cosh(kl*vl(2))
      tsinh=sqrt(coshkvl2*coshkvl2-1)
      if(vl(2).ge.zero) then
        sinhkvl2=tsinh
      else
        sinhkvl2=-tsinh
      endif    
  
      coshkvr2=cosh(kl*vr(2))
      tsinh=sqrt(coshkvr2*coshkvr2-1)
      if(vr(2).ge.zero) then
        sinhkvr2=tsinh
      else
        sinhkvr2=-tsinh
      endif 

      coshkvl1=cosh(kl*vl(1))
      tsinh=sqrt(coshkvl1*coshkvl1-1)
      if(vl(1).ge.zero) then
        sinhkvl1=tsinh
      else
        sinhkvl1=-tsinh
      endif 

      coshkvr1=cosh(kl*vr(1))
      tsinh=sqrt(coshkvr1*coshkvr1-1)
      if(vr(1).ge.zero) then
        sinhkvr1=tsinh
      else
        sinhkvr1=-tsinh
      endif     
  
      thc(1)=(coshkvl2-coshkvr2)*dif(2)*rkl*hc2d(1,l,k,m)     
      thc(2)=(sinhkvl2-sinhkvr2)*dif(2)*rkl*hc2d(2,l,k,m)     
      thc(3)=(coshkvl1-coshkvr1)*dif(1)*rkl*hc2d(3,l,k,m)     
      thc(4)=(sinhkvl1-sinhkvr1)*dif(1)*rkl*hc2d(4,l,k,m)     

      bunmo=2.*dif(1)*dif(2)*rkl*rkl
      do i=1,2 ! 1:y 2:x
        coshm(i)=cosh(ml*vl(i))-cosh(ml*vr(i))
        sinhm(i)=sinh(ml*vl(i))-sinh(ml*vr(i))
      enddo
      thc(5)=bunmo*coshm(2)*sinhm(1)*hc2d(5,l,k,m)
      thc(6)=bunmo*coshm(2)*coshm(1)*hc2d(6,l,k,m)
      thc(7)=bunmo*sinhm(2)*coshm(1)*hc2d(7,l,k,m)
      thc(8)=bunmo*sinhm(2)*sinhm(1)*hc2d(8,l,k,m)
      
      tsum=0
      tsum=tsum+tpc(0)
      do i=1,14
        if(i>=9) goto 120
        tsum=tsum+thc(i)      
  120   tsum=tsum+tpc(i)
      enddo
      if(tsum.lt.zero) then
!FIXME uncomment
!        print *, "NAGATIVE PIN FLUX at ", l,k,m
        tsum=zero
      endif
      phihav=tsum
!     
      return     
      
    end function