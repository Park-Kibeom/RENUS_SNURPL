! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
	subroutine onenemz                      &
          (reigv,cntz0i,cntz1i,dtl,dtlu,dtll, &
           dc,dcu,dcl,                        &
           hz,hzu,hzl,hs,sigr,signf,sigs,     &
           aflx,                              &
           cntz0o,cntz1o,srcz)


	implicit double precision (a-h,o-z)

! nem solver for z-direction of hexagonal node  
! axial direction source calculation to solve radial direction 
! input : cntz0i - incoming current from lower node in z-direction
!         cntz1i - incoming current from upper node in z-direction
!         aleaku - average trans-leak of upper hexagon
!         dcu - dc(diffusion coeff) of upper node
!         aleakl - average trans-leak of lower hexagon
!         dcl - dc(diffusion coeff) of lower node
! Working : aflx - node average flux
! output : cntzo - outgoing current in z-direction
!          lkz - delta cross section corresponding to axial leakage
!          srcz - source from axial direction
      dimension cntz0i(2),cntz1i(2)
      dimension dtl(2),dtlu(2),dtll(2),tl0(2),tl1(2)
      dimension dc(2),dcu(2),dcl(2)  
      dimension sigr(2),signf(2)  
      dimension aflx(2)
      dimension cntz0o(2),cntz1o(2),srcz(2)
      dimension alp1(2),alp2(2)
      dimension sbal(2),smm1(2),smm2(2),sj0(2),sj1(2)
! Transverse Leakage and coeff's at interface
! 2013_07_10 . scb
#define scb_correction_tl
#ifdef scb_correction_tl
      real :: ratl, ratr, difsl(2), difsr(2)
      
      ratl=hzl/hz
      ratr=hzu/hz
      
      difsl=dtll-dtl
      difsr=dtlu-dtl
      
      if(hzl.lt.1.e-10) then
        alp1(:)=(dtl(:)+difsr(:)/((1+ratr)*(1+2*ratr))) / (1+1/(1+2*ratr))
        alp2(:)=(-dtl(:)+difsr(:)/(1+ratr)) / (2*(1+ratr))
        !print *, 'left : vacuum'
      elseif(hzu.lt.1.e-10) then
        alp1(:)=-(dtl(:)+difsl(:)/((1+ratl)*(1+2*ratl))) / (1+1/(1+2*ratl))
        alp2(:)=(-dtl(:)+difsl(:)/(1+ratl)) / (2*(1+ratl))       
        !print *, 'right : vacuum'
      else
        alp1(:)=(-difsl(:)/((1+ratl)*(1+2*ratl))+difsr(:)/((1+ratr)*(1+2*ratr))) / (1/(1+2*ratl)+1/(1+2*ratr))
        alp2(:)=(difsl(:)/(1+ratl)+difsr(:)/(1+ratr)) / (2*(1+ratl+ratr))        
      endif
 
      do ig=1,2
        sbal(ig)=-dtl(ig)
        smm1(ig)=-alp1(ig)/3.
        smm2(ig)=-alp2(ig)/5.
      enddo
! added end      
#else
      do ig=1,2
!        tl0(ig)=(dtll(ig)/hzl+dtl(ig)/hz)
!     +             /(dcl(ig)/hzl+dc(ig)/hz)
!        tl1(ig)=(dtlu(ig)/hzu+dtl(ig)/hz)
!     +             /(dcu(ig)/hzu+dc(ig)/hz)
        tl0(ig)=(dtll(ig)/hzl+dtl(ig)/hz) &
                  /(dc(ig)/hzl+dc(ig)/hz)
        tl1(ig)=(dtlu(ig)/hzu+dtl(ig)/hz) &
                  /(dc(ig)/hzu+dc(ig)/hz)
      enddo
      do ig=1,2
        alp1(ig)=(tl1(ig)-tl0(ig))/2.
        alp2(ig)=dtl(ig)/dc(ig)-(tl1(ig)+tl0(ig))/2.
      enddo

      do ig=1,2
        sbal(ig)=-dtl(ig)
        smm1(ig)=-alp1(ig)*dc(ig)/3.
        smm2(ig)=-alp2(ig)*dc(ig)/5.
      enddo
#endif      
      sbal(1)=sbal(1)+(signf(1)*aflx(1)+signf(2)*aflx(2))*reigv
      sl1=sigr(1)-signf(1)*reigv
      sr1=sigr(1)
      sr2=sigr(2)
      sf2=signf(2)*reigv
      ss=sigs
      d1=dc(1)
      d2=dc(2)
      h=hz

      rdt=-1./sr1/sr2
      h2=h*h
      h3=h*h*h
      h4=h*h*h*h
      rh=1./h
      a11=10.*d1+h2*sl1
      a12=10.*d2+h2*sr2
      a21=42.*d1+h2*sl1
      a22=42.*d2+h2*sr2
      a3=h4*sf2*ss-a11*a12
      a4=h4*sf2*ss-a21*a22
      ra3=1./a3
      ra4=1./a4
      c1=-rh/sr1+h/60.*(-4./d1+125.*a12*ra3+147.*a22*ra4)
      c2=-rh/sr1+h/60.*(1./d1-125.*a12*ra3+147.*a22*ra4)
      c3=sf2*h3*(25./12.*ra3+49./20.*ra4)
      c4=sf2*h3*(-25./12.*ra3+49./20.*ra4)
      c5=ss*(rdt*rh+h3*(25./12.*ra3+49./20.*ra4))
      c6=ss*(rdt*rh+h3*(-25./12.*ra3+49./20.*ra4))
      c7=-rh/sr2+h/60.*(-4./d2+125.*a11*ra3+147*a21*ra4)
      c8=-rh/sr2+h/60.*(1./d2-125.*a11*ra3+147*a21*ra4)

      sj1(1)=-(2.+c1)*cntz1i(1)-c2*cntz0i(1)-c3*cntz1i(2)-c4*cntz0i(2) &
              +1./sr1*sbal(1)                                          &
              -2.5*(a12*smm1(1)+h2*sf2*smm1(2))*h2*ra3                 &
              +3.5*(a22*smm2(1)+h2*sf2*smm2(2))*h2*ra4
      sj0(1)=-(2.+c1)*cntz0i(1)-c2*cntz1i(1)-c3*cntz0i(2)-c4*cntz1i(2) &
              +1./sr1*sbal(1)                                          &
              +2.5*(a12*smm1(1)+h2*sf2*smm1(2))*h2*ra3                 &
              +3.5*(a22*smm2(1)+h2*sf2*smm2(2))*h2*ra4                 
      sj1(2)=-c5*cntz1i(1)-c6*cntz0i(1)-(2.+c7)*cntz1i(2)-c8*cntz0i(2) &
              -(ss*sbal(1)+sr1*sbal(2))*rdt                            &
              -2.5*(h2*ss*smm1(1)+a11*smm1(2))*h2*ra3                  &
              +3.5*(h2*ss*smm2(1)+a21*smm2(2))*h2*ra4              
      sj0(2)=-c5*cntz0i(1)-c6*cntz1i(1)-(2.+c7)*cntz0i(2)-c8*cntz1i(2) &
              -(ss*sbal(1)+sl1*sbal(2))*rdt                            &
              +2.5*(h2*ss*smm1(1)+a11*smm1(2))*h2*ra3                  &
              +3.5*(h2*ss*smm2(1)+a21*smm2(2))*h2*ra4
!      tm1=4.-2.*c1+2.*c2-c3*c5+c4*c5+c3*c6-c4*c6-2.*c7+c1*c7
!     +     -c2*c7+2.*c8-c1*c8+c2*c8
!      tm2=4.-2.*c1-2.*c2-c3*c5-c4*c5-c3*c6-c4*c6-2.*c7+c1*c7+c2*c7 
!     +     -2.*c8+c1*c8+c2*c8
!      rdt2=1./tm1/tm2
      c9=-2.+c1
      c0=-2.+c7
      tm1=(c0-c8)*(c2-c9)+(c3-c4)*(c5-c6)
      tm2=(c0+c8)*(c2+c9)-(c3+c4)*(c5+c6)
      rdt2=-1./tm1/tm2
!      cio=rdt2*((-2.+c7)*(-4.+c3*c5+c4*c6-c1*(-2.+c7)+2.*c7)
!     +         -(c4*c5+c3*c6)*c8+(-2.+c1)*c8*c8)
!      ctr=rdt2*(-(c4*c5+c3*c6)*(-2.+c7)+(c3*c5+c4*c6)*c8
!     +          +c2*((-2.+c7)*(-2.+c7)-c8*c8))
!      cgo=rdt2*(-c3*c3*c5+(-2.+c1)*c3*(-2.+c7)+c2*c3*c8+
!     +          +c4*(c4*c5-c2*(-2.+c7)+2.*c8-c1*c8))
!      cgt=rdt2*(-c4*c4*c6+(-2.+c1)*c4*(-2.+c7)+c3*(c3*c6+2.*c8-c1*c8)
!     +          +c2*(-c3*(-2.+c7)+c4*c8))
      t1=c3*c5+c4*c6
      t2=c4*c5+c3*c6
      t3=c0*c2+c8*c9
      t4=c2*c8+c0*c9
      cio=rdt2*(c0*(t1-c0*c9)+c8*(c8*c9-t2))
      ctr=rdt2*(c0*(c0*c2-t2)+c8*(t1-c2*c8))
      cgo=rdt2*(c4*(c4*c5-t3)+c3*(t4-c3*c5))
      cgt=rdt2*(c3*(c3*c6-t3)+c4*(t4-c4*c6))
      cntz1o(1)=cio*sj1(1)+ctr*sj0(1)+cgo*sj1(2)+cgt*sj0(2)
      cntz0o(1)=cio*sj0(1)+ctr*sj1(1)+cgo*sj0(2)+cgt*sj1(2)
      cgo=rdt2*(c6*(c3*c6-t3)+c5*(t4-c3*c5))
      cgt=rdt2*(c5*(c4*c5-t3)+c6*(t4-c4*c6))
      cio=rdt2*(c2*(c0*c2-t2)+c9*(t1-c0*c9))
      ctr=rdt2*(c2*(t1-c2*c8)+c9*(c8*c9-t2))
!      cgo=rdt2*(-c3*c5*c5+(-2.+c1)*c5*(-2.+c7)+c2*c5*c8
!     +          +c6*(c3*c6-c2*(-2.+c7)+2.*c8-c1*c8))
!      cgt=rdt2*(c4*(c5-c6)*(c5+c6)+(-2.+c1)*(c6*(-2.+c7)-c5*c8)
!     +          +c2*(-c5*(-2.+c7)+c6*c8))
!      cio=rdt2*(-c2*(c4*c5+c3*c6)+(c2*c2-c1*c1)*(-2.+c7)
!     +          +(c1-2.)*(-4.+c3*c5+c4*c6+2.*c7)+c1*(-4.+2.*c7))
!      ctr=rdt2*(c2*(c3*c5+c4*c6)-c2*c2*c8
!     +          -(-2.+c1)*(c4*c5+c3*c6+2.*c8-c1*c8))
      cntz1o(2)=cio*sj1(2)+ctr*sj0(2)+cgo*sj1(1)+cgt*sj0(1)
      cntz0o(2)=cio*sj0(2)+ctr*sj1(2)+cgo*sj0(1)+cgt*sj1(1)
! outgoing partial current calculation
      do ig=1,2
        srcz(ig)=(cntz0i(ig)+cntz1i(ig)-cntz0o(ig)-cntz1o(ig))*rh
      enddo
 
      return
      end 
