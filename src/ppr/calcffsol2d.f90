subroutine calcffsol2d(m,l,k,kappa1)
    use ppr
    use geom,   only : hmesh
    use xsec,   only : xstf, xsdf
    implicit none

    integer,intent(in)  :: m,l,k
    real                :: kappa1

    real                :: rh2,sigr,rsigr,          &
                           rsigr2,rh2rsigr2,        &
                           sigdc,rsigdc,rsigdc2,    &
                           cko,sko,rsko,            &
                           ckt,skt,ckt2,skt2,       &
                           temp,rbk,rxt,rbkckosko,  &
                           r2rsigdc,bkh

! for coefficients of particular solutions
! square cell only
! rh
      rh2=1/(hmesh(XDIR,l,k)*hmesh(XDIR,l,k))

! rxstf
      sigr=xstf(m,l,k)
      rsigr=1/sigr
      rsigr2=rsigr*rsigr

! rh*rh*rxstf*rxstf
      rh2rsigr2=rh2*rsigr2  

! xsdf*xsdf
      sigdc=xsdf(m,l,k)
      rsigdc=1/sigdc
      rsigdc2=rsigdc*rsigdc
      
! constant variable for coefficients of particular solutions
      cpc02(l,k,m)=12*sigdc*rh2rsigr2
      cpc04(l,k,m)=4*sigdc*rh2rsigr2*(420*sigdc*rh2*rsigr+10)
      cpc022(l,k,m)=288*sigdc*sigdc*rh2rsigr2*rh2*rsigr
      temp=12*sigdc*rh2rsigr2
      cpc11(l,k,m)=5*temp
      cpc12(l,k,m)=temp
      temp=4*sigdc*rh2rsigr2
      cpc21(l,k,m)=35*temp
      cpc22(l,k,m)=3*temp
!
! for coefficients of homogeneous solutions
      rbk=1/kappa1
      cko=cosh(kappa1)          ! kappa-Original
      sko=sqrt(cko*cko-1)
      rsko=1/sko 
      ckt=cosh(kappa1*rsqrt2)   ! kappa-Tilda
      skt=sqrt(ckt*ckt-1)
  
      ckt2=ckt*ckt
      skt2=skt*skt
      rxt=1/(skt*ckt)      
      rbkckosko=1/(kappa1*cko-sko)
      bkh=kappa1*hmesh(XDIR,l,k)
      r2rsigdc=0.5*rsigdc

      chc6(l,k,m)=1/skt2
      chc13j(l,k,m)=bkh*r2rsigdc*(-rbkckosko)
      chc13p(l,k,m)=-rbkckosko
      chc57j(l,k,m)=bkh*sko*r2rsigdc*rbkckosko*rxt
      chc57p(l,k,m)=kappa1*cko*rbkckosko*rxt
      temp=kappa1*1/(kappa1*sko*ckt2-2*cko*skt2)
      chc8j(l,k,m)=temp*cko*hmesh(XDIR,l,k)*0.5*rsigdc
      chc8p(l,k,m)=temp*sko
      chc24j(l,k,m)=-hmesh(XDIR,l,k)*r2rsigdc*rsko
      chc24a(l,k,m)=-skt2*rsko*rbk
!
      return
end subroutine