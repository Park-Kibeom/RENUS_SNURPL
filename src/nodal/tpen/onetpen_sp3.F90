! added in ARTOS ver. 0.2 . 2012_07_24 by SCB        
   subroutine onetpen_sp3(       &
         xsr,xst,   &
         xsd,xsd2,     &
         hside,pflx,cnti,     &
         aflx,xmom,ymom,       &
         cnto,hflx,            &
         asrc,xsrc,ysrc,  &
         sflb)

      implicit none

! input arguements
      real :: hside
      real :: xsr, xst, xsd, xsd2
      real :: pflx(2,6), cnti(2,6)
      real :: aflx(2,6), xmom(2,6), ymom(2,6)
      real :: cnto(2,6), hflx(2)

      real :: sfli(2,6), sflb(2,6)

! local variables
      real :: caxy(4), ca2(2), ca3(2), raxy(4)
      real :: cx2(2), cx3(2), cx4(2)
      real :: cy2(2), cy3(2)
      real :: cj1(2), cj2(2), cj3(4), cj4(2)
      real :: bsrc(2,6), ssrc(2,6), psrc(2)
      real :: asrc(2,6), xsrc(2,6), ysrc(2,6)

! condensed matrix(13x13)
      real :: csb1(4), csb2(4), csb3(4)
      real :: csi1(4), csi2(4), csi3(4), csi4(4)
      real :: cpf1(4), cpf2(4), cpf3(4)

      real :: j1ra2(4), ra2(4), j1r(4)
      real :: det, sfa(2), sfx(2), sfy(2), sfax(2)

! condensed matrix(7x7)
      real :: rsb(4)
      real :: t(4), q1(4), q2(4)
      real :: c(2), fac1(4), fac2(4), fac(4), rcsi(4)
      real :: cpfl(2)
      real :: addsfl1, addsfl2
      real :: addsfl, subsfl, ocur


      integer :: im, im2, it
      real :: r80h2, r32h2, h2, r2o3, r4o3, r5o3, r1o6
      real :: r20r3oh, r8r3oh, r1o12, r1o18, r1o3

	   integer :: mp1(6),mp2(6),mp3(6),mp4(6),mp5(6)
	   integer :: mg1(4),mg2(4),mg3(4),mg4(4),mg5(4) 
	   data mp1/2,3,4,5,6,1/ 
	   data mp2/3,4,5,6,1,2/ 
	   data mp3/4,5,6,1,2,3/ 
	   data mp4/5,6,1,2,3,4/ 
	   data mp5/6,1,2,3,4,5/
	   data mg1/1,1,2,2/
	   data mg2/1,2,1,2/
	   data mg3/1,1,3,3/
	   data mg4/2,2,4,4/
	   data mg5/3,4,3,4/

      h2=hside*hside
      r80h2=80./h2; r32h2=32./h2
      r2o3=2./3.; r4o3=4./3.; r5o3=5./3.; r1o6=1./6.
      r1o12=1./12.; r1o18=1./18.; r1o3=1./3.
      r20r3oh=20./sqrt(3.)/hside
      r8r3oh=8.*sqrt(3.)/hside

! FIX_ME ::
#ifdef DEBUG
      pflx(1,:)=20.
      pflx(2,1)=-5.
      pflx(2,2)=-6.
      pflx(2,3)=-7.
      pflx(2,4)=-8.
      pflx(2,5)=-9.
      pflx(2,6)=-10.

      cnti(1,:)=10.
      cnti(2,:)=-5.

      reigv=1.

      asrc=0.
      xsrc=0.
      ysrc=0.

      caxy(1)=r80h2*xsd(1)+xsr(1)-reigv*xsnf(1)
      caxy(2)=-reigv*xsnf(2)
      caxy(3)=-xssf(1,2)
      caxy(4)=r80h2*xsd(2)+xsr(2)

      cj3(1)=0.5+r8r3oh*xsd(1)
      cj3(2)=0.
      cj3(3)=0.
      cj3(4)=0.5+r8r3oh*xsd(2)

      ca2(1)=-r32h2*xsd(1)
      ca2(2)=-r32h2*xsd(2)

      cj1(1)=-r20r3oh*xsd(1)
      cj1(2)=-r20r3oh*xsd(2)
#endif

! initialization
      caxy(1)=r80h2*xsd+xsr
      caxy(2)=-2.*xsr
      caxy(3)=-r2o3*xsr
      caxy(4)=r80h2*xsd2+r4o3*xsr+r5o3*xst

      cj3(1)=0.5+r8r3oh*xsd
      cj3(2)=-0.375
      cj3(3)=-0.125
      cj3(4)=0.875+r8r3oh*xsd2

      ca2(1)=-r32h2*xsd
      ca2(2)=-r32h2*xsd2

      cj1(1)=-r20r3oh*xsd
      cj1(2)=-r20r3oh*xsd2

      do im=1,2
         ca3(im)=-r1o6*ca2(im)
         cx2(im)=r1o6*ca2(im)
         cx3(im)=-r1o12*ca2(im)
         cx4(im)=-r1o18*ca2(im)
         cy2(im)=0.25*ca2(im)
         cy3(im)=-cy2(im)
         cj4(im)=0.1*cj1(im)
      enddo

      do it=1,6
         do im=1,2
            asrc(im,it)=asrc(im,it)-ca3(im)*(pflx(im,it)+pflx(im,mp1(it)))
            xsrc(im,it)=xsrc(im,it)+0.5*cx4(im)*(pflx(im,it)+pflx(im,mp1(it)))
            ysrc(im,it)=ysrc(im,it)-r1o3*cy2(im)*(pflx(im,it)-pflx(im,mp1(it)))
                    
            bsrc(im,it)=-0.5*cj4(im)*(pflx(im,it)+pflx(im,mp1(it)))+2.*cnti(im,it)
            ssrc(im,it)=pflx(im,it)+pflx(im,mp1(it))+pflx(im,mp5(it))
         enddo
      enddo

! condensed matrix(13x13)
      det=1./(caxy(1)*caxy(4)-caxy(2)*caxy(3))
      raxy(1)=det*caxy(4)
      raxy(2)=-det*caxy(2)
      raxy(3)=-det*caxy(3)
      raxy(4)=det*caxy(1)

      j1ra2(1)=cj1(1)*raxy(1)*ca2(1)
      j1ra2(2)=cj1(1)*raxy(2)*ca2(2)
      j1ra2(3)=cj1(2)*raxy(3)*ca2(1)
      j1ra2(4)=cj1(2)*raxy(4)*ca2(2)

      ra2(1)=ca2(1)*raxy(1)
      ra2(2)=ca2(2)*raxy(2)
      ra2(3)=ca2(1)*raxy(3)
      ra2(4)=ca2(2)*raxy(4)

      j1r(1)=cj1(1)*raxy(1)
      j1r(2)=cj1(1)*raxy(2)
      j1r(3)=cj1(2)*raxy(3)
      j1r(4)=cj1(2)*raxy(4)

! left had side
! line 1
      do im=1,4
         csb1(im)=cj3(im)-2.*j1ra2(im)
         csb2(im)=-0.5*j1ra2(im)
         csb3(im)=0.5*j1ra2(im)
      enddo
      csb3(1)=csb3(1)+0.1*cj1(1)
      csb3(4)=csb3(4)+0.1*cj1(2)

! line 2
      do im=1,4
         csi1(im)=5.*ra2(im)
         csi2(im)=40.*ra2(im)
         csi3(im)=5.*ra2(im)
      enddo
      csi2(1)=csi2(1)+24.
      csi2(4)=csi2(4)+24.
      csi4(1)=-1.
      csi4(2)=0.
      csi4(3)=0.
      csi4(4)=-1.

! line 3
      do im=1,4
         cpf1(im)=-2.5*ra2(im)
         cpf2(im)=2.5*ra2(im)
         cpf3(im)=5.*ra2(im)
      enddo
      cpf1(1)=cpf1(1)-1.   
      cpf1(4)=cpf1(4)-1.  
      cpf3(1)=cpf3(1)+6.  
      cpf3(4)=cpf3(4)+6. 

! right hand side
      psrc=0.
      do it=1,6
         do im=1,2
            sfax(im)=asrc(im,it)+6.*xsrc(im,it)
            sfa(im)=asrc(im,it)+asrc(im,mp5(it))
            sfx(im)=-xsrc(im,it)-xsrc(im,mp5(it))
            sfy(im)=ysrc(im,it)-ysrc(im,mp5(it))
         enddo

         bsrc(1,it)=bsrc(1,it)-j1r(1)*sfax(1) &
                                -j1r(2)*sfax(2)
         bsrc(2,it)=bsrc(2,it)-j1r(3)*sfax(1) &
                                -j1r(4)*sfax(2)

         ssrc(1,it)=ssrc(1,it)+raxy(1)*(10.*sfa(1)+30.*sfx(1)+30.*sfy(1)) &
                              +raxy(2)*(10.*sfa(2)+30.*sfx(2)+30.*sfy(2))
         ssrc(2,it)=ssrc(2,it)+raxy(3)*(10.*sfa(1)+30.*sfx(1)+30.*sfy(1)) &
                              +raxy(4)*(10.*sfa(2)+30.*sfx(2)+30.*sfy(2))
         psrc(1)=psrc(1)-15.*(raxy(1)*xsrc(1,it)+raxy(2)*xsrc(2,it))
         psrc(2)=psrc(2)-15.*(raxy(3)*xsrc(1,it)+raxy(4)*xsrc(2,it))
      enddo

! condensed matrix(7x7)
      det=1./(csb1(1)*csb1(4)-csb1(2)*csb1(3))
      rsb(1)=det*csb1(4)
      rsb(2)=-det*csb1(2)
      rsb(3)=-det*csb1(3)
      rsb(4)=det*csb1(1) 

      ! first row
      call axb(csi1,rsb,t)
      call bxc(t,csb2,q1)
      call bxc(t,csb3,q2)
      do im=1,4
         csi2(im)=csi2(im)-q1(im)*2.
         csi3(im)=csi3(im)-q1(im)
         csi4(im)=csi4(im)-q2(im)*2.
      enddo 
      do it=1,6
         c(1)=bsrc(1,it)+bsrc(1,mp5(it))
         c(2)=bsrc(2,it)+bsrc(2,mp5(it))
         ssrc(1,it)=ssrc(1,it)-t(1)*c(1)-t(2)*c(2)
         ssrc(2,it)=ssrc(2,it)-t(3)*c(1)-t(4)*c(2)
      enddo

      ! second row
      call axb(cpf1,rsb,t)
      call bxc(t,csb2,q1)
      call bxc(t,csb3,q2)
      do im=1,4
         cpf2(im)=cpf2(im)-2.*q1(im)
         cpf3(im)=cpf3(im)-6.*q2(im)
      enddo

      do im=1,2
         c(im)=0.
         do it=1,6
            c(im)=c(im)+bsrc(im,it)
         enddo
      enddo
      psrc(1)=psrc(1)-t(1)*c(1)-t(2)*c(2)
      psrc(2)=psrc(2)-t(3)*c(1)-t(4)*c(2)

! calculate the center point flx
      do im=1,4
         fac1(im)=csi2(im)+2.*csi3(im)
         fac2(im)=6.*csi4(im)
      enddo
      det=1./(fac1(1)*fac1(4)-fac1(2)*fac1(3))
      rcsi(1)=det*fac1(4)
      rcsi(2)=-det*fac1(2)
      rcsi(3)=-det*fac1(3)
      rcsi(4)=det*fac1(1)
      do im=1,4
         fac(im)=-cpf2(mg3(im))*rcsi(mg2(im)) &
                 -cpf2(mg4(im))*rcsi(mg5(im))
      enddo
      do im=1,4
         cpf3(im)=cpf3(im)+fac(mg3(im))*fac2(mg2(im)) &
                          +fac(mg4(im))*fac2(mg5(im))
      enddo
      do im=1,2
         c(im)=0.
         do it=1,6
            c(im)=c(im)+ssrc(im,it)
         enddo
      enddo
      psrc(1)=psrc(1)+fac(1)*c(1)+fac(2)*c(2)
      psrc(2)=psrc(2)+fac(3)*c(1)+fac(4)*c(2)

      det=1./(cpf3(1)*cpf3(4)-cpf3(2)*cpf3(3))
      cpfl(1)=(cpf3(4)*psrc(1)-cpf3(2)*psrc(2))*det
      cpfl(2)=(-cpf3(3)*psrc(1)+cpf3(1)*psrc(2))*det

! calculate the inner surface fluxes
      do it=1,6
         ssrc(1,it)=ssrc(1,it)-csi4(1)*cpfl(1)-csi4(2)*cpfl(2)
         ssrc(2,it)=ssrc(2,it)-csi4(3)*cpfl(1)-csi4(4)*cpfl(2)
      enddo
      call sol_mat(csi2,csi3,ssrc,sfli)

! calculate the boundary surface fluxes  
      do it=1,6
         addsfl1=sfli(1,it)+sfli(1,mp1(it))
         addsfl2=sfli(2,it)+sfli(2,mp1(it))
         bsrc(1,it)=bsrc(1,it)-csb3(1)*cpfl(1)-csb3(2)*cpfl(2) &
                              -csb2(1)*addsfl1-csb2(2)*addsfl2
         bsrc(2,it)=bsrc(2,it)-csb3(3)*cpfl(1)-csb3(4)*cpfl(2) &
                              -csb2(3)*addsfl1-csb2(4)*addsfl2
         sflb(1,it)=rsb(1)*bsrc(1,it)+rsb(2)*bsrc(2,it)
         sflb(2,it)=rsb(3)*bsrc(1,it)+rsb(4)*bsrc(2,it)
      enddo

! calculate average fluxes, x-moments and y-moments
      do it=1,6
        do im=1,2
          addsfl=sfli(im,it)+sfli(im,mp1(it))
          subsfl=sfli(im,it)-sfli(im,mp1(it))
          ocur=sflb(im,it)

          ysrc(im,it)=ysrc(im,it)-cy2(im)*subsfl
          xsrc(im,it)=xsrc(im,it)-cx2(im)*ocur &
                      -cx3(im)*addsfl-cx4(im)*cpfl(im)
          asrc(im,it)=asrc(im,it)-ca2(im)*ocur &
                      -ca2(im)*addsfl-ca3(im)*cpfl(im)
        enddo
        ymom(1,it)=raxy(1)*ysrc(1,it)+raxy(2)*ysrc(2,it)
        ymom(2,it)=raxy(3)*ysrc(1,it)+raxy(4)*ysrc(2,it)

        xmom(1,it)=raxy(1)*xsrc(1,it)+raxy(2)*xsrc(2,it)
        xmom(2,it)=raxy(3)*xsrc(1,it)+raxy(4)*xsrc(2,it)

        aflx(1,it)=raxy(1)*asrc(1,it)+raxy(2)*asrc(2,it)
        aflx(2,it)=raxy(3)*asrc(1,it)+raxy(4)*asrc(2,it)
      enddo 

      do it=1,6
         cnto(1,it)=-cnti(1,it)+1./2.*sflb(1,it)-3./8.*sflb(2,it) 
         cnto(2,it)=-cnti(2,it)-1./8.*sflb(1,it)+7./8.*sflb(2,it)  
      enddo

      do im=1,2
         hflx(im)=0.
         do it=1,6 
            hflx(im)=hflx(im)+aflx(im,it)
         enddo
         hflx(im)=hflx(im)/6.
      enddo 
                                          
      return
   end subroutine 

   subroutine axb(a,b,d)
      implicit none
      real :: a(4), b(4), d(4)

      d(1)=a(1)*b(1)+a(2)*b(3)
      d(2)=a(1)*b(2)+a(2)*b(4)
      d(3)=a(3)*b(1)+a(4)*b(3)
      d(4)=a(3)*b(2)+a(4)*b(4)
         
      return
   end subroutine

   subroutine bxc(b,c,d)
      implicit none
      real :: b(4), c(4), d(4)

      d(1)=b(1)*c(1)+b(2)*c(3)
      d(2)=b(1)*c(2)+b(2)*c(4)
      d(3)=b(3)*c(1)+b(4)*c(3)
      d(4)=b(3)*c(2)+b(4)*c(4)
         
      return
   end subroutine

   subroutine sol_mat(a1,a2,b,x)
      implicit none
      real :: a1(4), a2(4), b(2,6), x(2,6)
      real :: ra(4), ttt(2)
      real :: det, a5, a6, a7, a8, d1, d2, d3, d4
      real :: c1, c2, c3, c4, c5, c6, c7, c8
      real :: e1, e2, e3, e4, e5, e6, e7, e8
      real :: addsrc1, addsrc2, addsrc3, addsrc4
      real :: dat1, dat2, dat3, dat4
      integer :: it

      integer :: mp1(6),mp2(6),mp3(6),mp4(6),mp5(6)
      data mp1/2,3,4,5,6,1/ 
      data mp2/3,4,5,6,1,2/ 
      data mp3/4,5,6,1,2,3/ 
      data mp4/5,6,1,2,3,4/ 
      data mp5/6,1,2,3,4,5/

      det=1./(a1(1)*a1(4)-a1(2)*a1(3))
      ra(1)=det*a1(4)
      ra(2)=-det*a1(2)
      ra(3)=-det*a1(3)
      ra(4)=det*a1(1)
      a5=ra(1)*a2(1)+ra(2)*a2(3)
      a6=ra(1)*a2(2)+ra(2)*a2(4)
      a7=ra(3)*a2(1)+ra(4)*a2(3)
      a8=ra(3)*a2(2)+ra(4)*a2(4)
      do it=1,6
        ttt(1)=b(1,it)
        ttt(2)=b(2,it)
        b(1,it)=ra(1)*ttt(1)+ra(2)*ttt(2)
        b(2,it)=ra(3)*ttt(1)+ra(4)*ttt(2)
      enddo
      d1=1./6./((1.-a8)*(1.-a5)-a6*a7)
      d2=1./6./((1.+a8)*(1.+a5)-a6*a7)
      d3=1./6./((1.-2.*a8)*(1.-2.*a5)-4.*a6*a7)
      d4=1./6./((1.+2.*a8)*(1.+2.*a5)-4.*a6*a7)
      c1=2.*(d1+d2)+d3+d4+2.*a8*(-d1+d2-d3+d4)
      c2=2.*a6*(d1-d2+d3-d4)
      c3=-d1+d2-d3+d4+a8*(d1+d2+2.*(d3+d4))
      c4=-a6*(d1+d2+2.*(d3+d4))
      c5=-d1-d2+d3+d4+a8*(d1-d2-2.*(d3-d4))
      c6=a6*(-d1+d2+2.*(d3-d4))
      c7=2.*(d1-d2)-d3+d4+2.*a8*(-d1-d2+d3+d4)
      c8=2.*a6*(d1+d2-d3-d4)
      e1=2.*a7*(d1-d2+d3-d4)
      e2=2.*(d1+d2)+d3+d4+2.*a5*(-d1+d2-d3+d4)
      e3=-a7*(d1+d2+2.*(d3+d4))
      e4=-d1+d2-d3+d4+a5*(d1+d2+2.*(d3+d4))
      e5=a7*(-d1+d2+2.*(d3-d4))
      e6=-d1-d2+d3+d4+a5*(d1-d2-2.*(d3-d4))
      e7=2.*a7*(d1+d2-d3-d4)
      e8=2.*(d1-d2)-d3+d4+2.*a5*(-d1-d2+d3+d4)
      do it=1,6
        addsrc1=b(1,mp1(it))+b(1,mp5(it))
        addsrc2=b(2,mp1(it))+b(2,mp5(it))
        addsrc3=b(1,mp2(it))+b(1,mp4(it))
        addsrc4=b(2,mp2(it))+b(2,mp4(it))
        dat1=b(1,it)
        dat2=b(2,it)
        dat3=b(1,mp3(it))
        dat4=b(2,mp3(it))
        x(1,it)=c1*dat1+c2*dat2+c3*addsrc1+c4*addsrc2 &
                  +c5*addsrc3+c6*addsrc4+c7*dat3+c8*dat4
        x(2,it)=e1*dat1+e2*dat2+e3*addsrc1+e4*addsrc2 &
                  +e5*addsrc3+e6*addsrc4+e7*dat3+e8*dat4
      enddo

      return
   end subroutine




