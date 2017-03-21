! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
!#define DEBUG  
      subroutine onetpen(                      &
          reigv,sigr,signf,sigs,dc             &
         ,pflx,cnti,sz,srcbal,srcmomx,srcmomy  &
         ,cnto,hflx,errh)

	   implicit double precision (a-h,o-z)

! HOPEN(Higher Order Polynomial Expansion Nodal method) solver for hexagon 
!    using HOPEN solution based on Triangle 
!    input : boundary incoming current, boundary point flux
!    output : boundary outgoing current  
!  
!    pflx : boundary point flux(input) 
!    cnti : boundary incoming current(input) 
!    aflx : node average flux at previous step 
!    xmom : x-moment at previous step 
!    ymom : y-moment at previous step 
!    sfli : inner surface flux
!    cpfl : center point flux
!    cnto : boundary outgoing current(output) 
!
!                        2     
!                   2---------3 
!                   / \  2  / \ 
!                 1/   2   3   \3 
!                 /  1  \ /   3 \  
!                1--1--------4---4 
!                 \   6 / \   4 / 
!                 6\   6   5   /4  
!                   \ /  5  \ /   
!                   6---------5 
!                        5             
!                3
!               / \                     
!           2 /     \ 3                   
!        2  /         \  4                 
!          l           l                      
!        1 l           l 4                     
!          l           l                         
!        1  \         /  5                  
!           6 \     / 5                    
!               \ /                          
!                6                      
!                                      
!
!
      parameter (myn=6,mxs=6,mxo=6,mxp=1)

      dimension pflx(2,mxo),cnti(2,mxo),cnto(2,mxo),aflx(2,myn),hflx(2)
      dimension xmom(2,myn),ymom(2,myn),sfli(2,mxs),cpfl(2)  
      dimension asrc(2,myn),xsrc(2,myn),ysrc(2,myn),ssrc(2,mxs)
      dimension ojsrc(2,mxo),psrc(2)
      dimension sigr(2),signf(2),dc(2)
      dimension ca(2,3),cx(2,3),cy(2),cs(2,5),cj(2,4),cp(2,3),caxy(4)
      dimension csp(4,4),cjp(4,3),cpp(4,3)
      dimension rcaxy(4),rcjp(4),rcsp(4), bsflx(2,6)
      dimension oflx(2,myn),fac(4),fac1(4),fac2(4),fac3(4),ttt(2)
      dimension srcbal(2,myn),srcmomx(2,myn),srcmomy(2,myn)
      dimension mp1(6),mp2(6),mp3(6),mp4(6),mp5(6)
      dimension mg1(4),mg2(4),mg3(4),mg4(4),mg5(4) 
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
      
! 2014_01_17 . scb      
      real,pointer :: aflxvtk(:,:,:,:)
      common / hexvtk / aflxvtk
! added end      

#ifdef DEBUG
      real :: rmat(62,62)=0.d0, rvec(62,1)=0.d0
      
      !pflx(1,:)=20.
      !pflx(2,1)=5.
      !pflx(2,2)=6.
      !pflx(2,3)=7.
      !pflx(2,4)=8.
      !pflx(2,5)=9.
      !pflx(2,6)=10.
      !
      !cnti(1,:)=10.
      !cnti(2,:)=-5.
      !
      !reigv=1.d0/1.02
      
918   format(a15,1p,e15.5)
919   format(a15,1p,2e15.5)      
      
      write(140918, '(a)')  "=========================  Input ========================="
      write(140918, 918) 'h(cm)', sz
      write(140918, 919) 'D(cm)', dc(1:2)
      write(140918, 919) 'XS_total(/cm)', sigr(1:2)
      write(140918, 919) 'XS_nfis(/cm)', signf(1:2)
      write(140918, 918) 'XS_s12(/cm)', sigs
      write(140918, 918) 'keff', 1.d0/reigv
      
      write(140918, 919) 'avg_src_1', srcbal(:,1)
      write(140918, 919) 'avg_src_2', srcbal(:,2)
      write(140918, 919) 'avg_src_3', srcbal(:,3)
      write(140918, 919) 'avg_src_4', srcbal(:,4)
      write(140918, 919) 'avg_src_5', srcbal(:,5)
      write(140918, 919) 'avg_src_6', srcbal(:,6)
      
      write(140918, 919) 'x_src_1', srcmomx(:,1)
      write(140918, 919) 'x_src_2', srcmomx(:,2)
      write(140918, 919) 'x_src_3', srcmomx(:,3)
      write(140918, 919) 'x_src_4', srcmomx(:,4)
      write(140918, 919) 'x_src_5', srcmomx(:,5)
      write(140918, 919) 'x_src_6', srcmomx(:,6)
      
      write(140918, 919) 'y_src_1', srcmomy(:,1)
      write(140918, 919) 'y_src_2', srcmomy(:,2)
      write(140918, 919) 'y_src_3', srcmomy(:,3)
      write(140918, 919) 'y_src_4', srcmomy(:,4)
      write(140918, 919) 'y_src_5', srcmomy(:,5)
      write(140918, 919) 'y_src_6', srcmomy(:,6)
      
      write(140918, 919) 'pflx_1', pflx(:,1)
      write(140918, 919) 'pflx_2', pflx(:,2)
      write(140918, 919) 'pflx_3', pflx(:,3)
      write(140918, 919) 'pflx_4', pflx(:,4)
      write(140918, 919) 'pflx_5', pflx(:,5)
      write(140918, 919) 'pflx_6', pflx(:,6)
      
      write(140918, 919) 'J_in_1', cnti(:,1)
      write(140918, 919) 'J_in_2', cnti(:,2)
      write(140918, 919) 'J_in_3', cnti(:,3)
      write(140918, 919) 'J_in_4', cnti(:,4)
      write(140918, 919) 'J_in_5', cnti(:,5)
      write(140918, 919) 'J_in_6', cnti(:,6)
      
      !stop
#endif


      rt3=sqrt(3.d0)
      
!#ifndef DEBUG
      do it=1,myn 
        asrc(1,it)=srcbal(1,it)
        xsrc(1,it)=srcmomx(1,it)
        ysrc(1,it)=srcmomy(1,it)
        asrc(2,it)=srcbal(2,it)
        xsrc(2,it)=srcmomx(2,it)
        ysrc(2,it)=srcmomy(2,it)
      enddo

      caxy(1)=sigr(1)-signf(1)*reigv
      caxy(2)=-signf(2)*reigv
      caxy(3)=-sigs
      caxy(4)=sigr(2)
      do ig=1,2
        bt=dc(ig)/sz/sz  
        caxy(3*ig-2)=caxy(3*ig-2)+80.*bt    ! caxy  : a1=x1=x2
        ca(ig,1)=-32.*bt                    ! ca(1) : a2
        !ca(ig,2)=-64.*bt                   
        ca(ig,2)=-32.*bt                    ! ca(2) : a3
        ca(ig,3)=16.*bt/3.                  ! ca(3) : a5=-1/6*a2
        cx(ig,1)=8.*bt/3.                   ! cx(1) : x2=x3
        !cx(ig,2)=-32.*bt/3.
        cx(ig,2)=-32.*bt/6.                 ! cx(2) : x4?
        cx(ig,3)=16.*bt/9.                  ! cx(3) : x5
        cy(ig)=-8.*bt                       ! cy    : y2

        cs(ig,1)=20. 
        cs(ig,2)=-60.
        cs(ig,3)=60.
        cs(ig,4)=-48.
        cs(ig,5)=2.

        tt=dc(ig)/rt3/sz
        cj(ig,1)=20.*tt                     ! cj(1) : alpha1
        cj(ig,2)=120.*tt                    ! cj(2) : alpha2
        !cj(ig,3)=-1.-48.*tt  
        cj(ig,3)=-0.5-24.*tt                ! cj(3) : 1/2 * alpha3

        cj(ig,4)=2.*tt                      ! cj(4) : alpha4
        cp(ig,1)=-15.
        !cp(ig,2)=2.
        cp(ig,2)=1.
        cp(ig,3)=-6.

        psrc(ig)=0.
        do it=1,mxo
!#ifdef DEBUG
!          asrc(ig,it)=asrc(ig,it)-ca(ig,2)*cnti(ig,it) &
!                      -ca(ig,3)*(pflx(ig,it)+pflx(ig,mp1(it)))
!          xsrc(ig,it)=xsrc(ig,it)-cx(ig,2)*cnti(ig,it) &
!                       +0.5*cx(ig,3)*(pflx(ig,it)+pflx(ig,mp1(it)))
!          ysrc(ig,it)=ysrc(ig,it) &
!                       +8.*bt/3.*(pflx(ig,it)-pflx(ig,mp1(it)))
!
!          ssrc(ig,it)=-2.*(pflx(ig,it) &
!                       +pflx(ig,mp1(it))+pflx(ig,mp5(it)))
!
!          ojsrc(ig,it)=(-1.+48.*tt)*cnti(ig,it) &
!                      -tt*(pflx(ig,it)+pflx(ig,mp1(it)))
!          psrc(ig)=psrc(ig)-2.*cnti(ig,it)
!#endif

          asrc(ig,it)=asrc(ig,it) &
                      -ca(ig,3)*(pflx(ig,it)+pflx(ig,mp1(it)))
          xsrc(ig,it)=xsrc(ig,it) &
                       +0.5*cx(ig,3)*(pflx(ig,it)+pflx(ig,mp1(it)))
          ysrc(ig,it)=ysrc(ig,it) &
                       +8.*bt/3.*(pflx(ig,it)-pflx(ig,mp1(it)))

          ssrc(ig,it)=-2.*(pflx(ig,it) &
                       +pflx(ig,mp1(it))+pflx(ig,mp5(it)))

          ojsrc(ig,it)=-2.*cnti(ig,it) &
                      -tt*(pflx(ig,it)+pflx(ig,mp1(it)))
        enddo
      enddo

      errh=0.
      do it=1,myn
        do ig=1,2
          oflx(ig,it)=aflx(ig,it)
        enddo
      enddo

      det=1./(caxy(1)*caxy(4)-caxy(2)*caxy(3))
      rcaxy(1)=det*caxy(4)
      rcaxy(2)=-det*caxy(2)
      rcaxy(3)=-det*caxy(3)
      rcaxy(4)=det*caxy(1)

      do it=1,3
        do ig=1,4
          csp(ig,it)=0
          cjp(ig,it)=0
          cpp(ig,it)=0
        enddo
      enddo
      do ig=1,4
        csp(ig,4)=0
      enddo
      do ig=1,4,3
        csp(ig,1)=cs(mg1(ig),4)
        csp(ig,4)=cs(mg1(ig),5)
        cjp(ig,2)=cj(mg1(ig),3)
        cjp(ig,3)=cj(mg1(ig),4)
        cpp(ig,2)=cp(mg1(ig),2)
        cpp(ig,3)=cp(mg1(ig),3)
      enddo

      do ig=1,4
        fac1(ig)=-rcaxy(ig)*cs(mg1(ig),1)
        fac2(ig)=-rcaxy(ig)*cs(mg1(ig),2)
        fac3(ig)=-rcaxy(ig)*cs(mg1(ig),3)
      enddo
      do ig=1,4
        sfax=ca(mg2(ig),1)*fac1(ig)+cx(mg2(ig),1)*fac2(ig)
        sfy=cy(mg2(ig))*fac3(ig)

        csp(ig,1)=csp(ig,1)+2.*(sfax+sfy)
        csp(ig,2)=csp(ig,2)+sfax-sfy

        csp(ig,3)=csp(ig,3)+ca(mg2(ig),2)*fac1(ig) &
                 +cx(mg2(ig),2)*fac2(ig)
        csp(ig,4)=csp(ig,4)+2.*(ca(mg2(ig),3)*fac1(ig) &
                 +cx(mg2(ig),3)*fac2(ig))
      enddo
      do it=1,myn
        sfa1=asrc(1,it)+asrc(1,mp5(it))
        sfa2=asrc(2,it)+asrc(2,mp5(it))
        sfx1=xsrc(1,it)+xsrc(1,mp5(it))
        sfx2=xsrc(2,it)+xsrc(2,mp5(it))
        sfy1=ysrc(1,it)-ysrc(1,mp5(it))
        sfy2=ysrc(2,it)-ysrc(2,mp5(it))
        ssrc(1,it)=ssrc(1,it)+fac1(1)*sfa1+fac1(2)*sfa2+fac2(1)*sfx1 &
                         +fac2(2)*sfx2+fac3(1)*sfy1+fac3(2)*sfy2
        ssrc(2,it)=ssrc(2,it)+fac1(3)*sfa1+fac1(4)*sfa2+fac2(3)*sfx1 &
                         +fac2(4)*sfx2+fac3(3)*sfy1+fac3(4)*sfy2
      enddo

      do ig=1,4
        fac1(ig)=-rcaxy(ig)*cj(mg1(ig),1)
        fac2(ig)=-rcaxy(ig)*cj(mg1(ig),2)
      enddo
      do ig=1,4
         cjp(ig,1)=cjp(ig,1)+ca(mg2(ig),1)*fac1(ig)+cx(mg2(ig),1)*fac2(ig)
         cjp(ig,2)=cjp(ig,2)+ca(mg2(ig),2)*fac1(ig)+cx(mg2(ig),2)*fac2(ig)
         cjp(ig,3)=cjp(ig,3)+ca(mg2(ig),3)*fac1(ig)+cx(mg2(ig),3)*fac2(ig)
      enddo
      do it=1,myn
        ojsrc(1,it)=ojsrc(1,it)+fac1(1)*asrc(1,it)+fac2(1)*xsrc(1,it) &
                               +fac1(2)*asrc(2,it)+fac2(2)*xsrc(2,it)
        ojsrc(2,it)=ojsrc(2,it)+fac1(3)*asrc(1,it)+fac2(3)*xsrc(1,it) &
                               +fac1(4)*asrc(2,it)+fac2(4)*xsrc(2,it)
      enddo
      do ig=1,4
        fac(ig)=-rcaxy(ig)*cp(mg1(ig),1)
      enddo
      do ig=1,4
        cpp(ig,1)=cpp(ig,1)+2.*cx(mg2(ig),1)*fac(ig)
        cpp(ig,2)=cpp(ig,2)+cx(mg2(ig),2)*fac(ig)
        cpp(ig,3)=cpp(ig,3)+6.*cx(mg2(ig),3)*fac(ig)
      enddo
      do ig=1,2
        ttt(ig)=0.
        do it=1,myn
          ttt(ig)=ttt(ig)+xsrc(ig,it)
        enddo
      enddo
      psrc(1)=psrc(1)+fac(1)*ttt(1)+fac(2)*ttt(2)
      psrc(2)=psrc(2)+fac(3)*ttt(1)+fac(4)*ttt(2)

      det=1./(cjp(1,2)*cjp(4,2)-cjp(2,2)*cjp(3,2))
      rcjp(1)=det*cjp(4,2)
      rcjp(2)=-det*cjp(2,2)
      rcjp(3)=-det*cjp(3,2)
      rcjp(4)=det*cjp(1,2)

      do ig=1,4
        fac(ig)=-csp(mg3(ig),3)*rcjp(mg2(ig)) &
                -csp(mg4(ig),3)*rcjp(mg5(ig))
      enddo
      do ig=1,4
        dat1=fac(mg3(ig))*cjp(mg2(ig),1)+fac(mg4(ig))*cjp(mg5(ig),1)
        dat2=fac(mg3(ig))*cjp(mg2(ig),3)+fac(mg4(ig))*cjp(mg5(ig),3)
        csp(ig,1)=csp(ig,1)+2.*dat1
        csp(ig,2)=csp(ig,2)+dat1
        csp(ig,4)=csp(ig,4)+2.*dat2
      enddo
      do it=1,myn
        dat1=ojsrc(1,it)+ojsrc(1,mp5(it))
        dat2=ojsrc(2,it)+ojsrc(2,mp5(it))
        ssrc(1,it)=ssrc(1,it)+fac(1)*dat1+fac(2)*dat2
        ssrc(2,it)=ssrc(2,it)+fac(3)*dat1+fac(4)*dat2
      enddo
      do ig=1,4
        fac(ig)=-cpp(mg3(ig),2)*rcjp(mg2(ig)) &
                -cpp(mg4(ig),2)*rcjp(mg5(ig))
      enddo
      do ig=1,4
        cpp(ig,1)=cpp(ig,1)+2.*fac(mg3(ig))*cjp(mg2(ig),1) &
                           +2.*fac(mg4(ig))*cjp(mg5(ig),1)
        cpp(ig,3)=cpp(ig,3)+6.*fac(mg3(ig))*cjp(mg2(ig),3) &
                           +6.*fac(mg4(ig))*cjp(mg5(ig),3)
      enddo
      do ig=1,2
        ttt(ig)=0.
        do it=1,myn
          ttt(ig)=ttt(ig)+ojsrc(ig,it)
        enddo
      enddo
      psrc(1)=psrc(1)+fac(1)*ttt(1)+fac(2)*ttt(2)
      psrc(2)=psrc(2)+fac(3)*ttt(1)+fac(4)*ttt(2)

      do ig=1,4
        fac1(ig)=csp(ig,1)+2.*csp(ig,2)
        fac2(ig)=6.*csp(ig,4)
      enddo
      rdet=1./(fac1(1)*fac1(4)-fac1(2)*fac1(3))
      rcsp(1)=rdet*fac1(4)
      rcsp(2)=-rdet*fac1(2)
      rcsp(3)=-rdet*fac1(3)
      rcsp(4)=rdet*fac1(1)
      do ig=1,4
        fac(ig)=-cpp(mg3(ig),1)*rcsp(mg2(ig)) &
                -cpp(mg4(ig),1)*rcsp(mg5(ig))
      enddo
      do ig=1,4
        cpp(ig,3)=cpp(ig,3)+fac(mg3(ig))*fac2(mg2(ig)) &
                           +fac(mg4(ig))*fac2(mg5(ig))
      enddo
      do ig=1,2
        ttt(ig)=0.
        do it=1,myn
          ttt(ig)=ttt(ig)+ssrc(ig,it)
        enddo
      enddo
      psrc(1)=psrc(1)+fac(1)*ttt(1)+fac(2)*ttt(2)
      psrc(2)=psrc(2)+fac(3)*ttt(1)+fac(4)*ttt(2)
      rdet=1./(cpp(1,3)*cpp(4,3)-cpp(2,3)*cpp(3,3))
      cpfl(1)=(cpp(4,3)*psrc(1)-cpp(2,3)*psrc(2))*rdet
      cpfl(2)=(-cpp(3,3)*psrc(1)+cpp(1,3)*psrc(2))*rdet

      do it=1,myn
        ssrc(1,it)=ssrc(1,it)-csp(1,4)*cpfl(1)-csp(2,4)*cpfl(2)
        ssrc(2,it)=ssrc(2,it)-csp(3,4)*cpfl(1)-csp(4,4)*cpfl(2)
      enddo

      det=1./(csp(1,1)*csp(4,1)-csp(2,1)*csp(3,1))
      rcsp(1)=det*csp(4,1)
      rcsp(2)=-det*csp(2,1)
      rcsp(3)=-det*csp(3,1)
      rcsp(4)=det*csp(1,1)
      a5=rcsp(1)*csp(1,2)+rcsp(2)*csp(3,2)
      a6=rcsp(1)*csp(2,2)+rcsp(2)*csp(4,2)
      a7=rcsp(3)*csp(1,2)+rcsp(4)*csp(3,2)
      a8=rcsp(3)*csp(2,2)+rcsp(4)*csp(4,2)
      do it=1,myn
        ttt(1)=ssrc(1,it)
        ttt(2)=ssrc(2,it)
        ssrc(1,it)=rcsp(1)*ttt(1)+rcsp(2)*ttt(2)
        ssrc(2,it)=rcsp(3)*ttt(1)+rcsp(4)*ttt(2)
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
      do it=1,myn
        addsrc1=ssrc(1,mp1(it))+ssrc(1,mp5(it))
        addsrc2=ssrc(2,mp1(it))+ssrc(2,mp5(it))
        addsrc3=ssrc(1,mp2(it))+ssrc(1,mp4(it))
        addsrc4=ssrc(2,mp2(it))+ssrc(2,mp4(it))
        dat1=ssrc(1,it)
        dat2=ssrc(2,it)
        dat3=ssrc(1,mp3(it))
        dat4=ssrc(2,mp3(it))
        sfli(1,it)=c1*dat1+c2*dat2+c3*addsrc1+c4*addsrc2 &
                  +c5*addsrc3+c6*addsrc4+c7*dat3+c8*dat4
        sfli(2,it)=e1*dat1+e2*dat2+e3*addsrc1+e4*addsrc2 &
                  +e5*addsrc3+e6*addsrc4+e7*dat3+e8*dat4
      enddo

      do it=1,myn
        addsfl1=sfli(1,it)+sfli(1,mp1(it))
        addsfl2=sfli(2,it)+sfli(2,mp1(it))
        ojsrc(1,it)=ojsrc(1,it)-cjp(1,3)*cpfl(1)-cjp(2,3)*cpfl(2) &
                   -cjp(1,1)*addsfl1-cjp(2,1)*addsfl2
        ojsrc(2,it)=ojsrc(2,it)-cjp(3,3)*cpfl(1)-cjp(4,3)*cpfl(2) &
                   -cjp(3,1)*addsfl1-cjp(4,1)*addsfl2
        cnto(1,it)=rcjp(1)*ojsrc(1,it)+rcjp(2)*ojsrc(2,it)
        cnto(2,it)=rcjp(3)*ojsrc(1,it)+rcjp(4)*ojsrc(2,it)
      enddo

      do it=1,myn
        do ig=1,2
          addsfl=sfli(ig,it)+sfli(ig,mp1(it))
          subsfl=sfli(ig,it)-sfli(ig,mp1(it))
          ocur=cnto(ig,it)
!#ifndef DEBUG
          ysrc(ig,it)=ysrc(ig,it)-cy(ig)*subsfl
          xsrc(ig,it)=xsrc(ig,it)-cx(ig,1)*addsfl &
                      -cx(ig,2)*ocur-cx(ig,3)*cpfl(ig)
!#endif
          asrc(ig,it)=asrc(ig,it)-ca(ig,1)*addsfl &
                      -ca(ig,2)*ocur-ca(ig,3)*cpfl(ig)
        enddo
!#ifndef DEBUG
        ymom(1,it)=rcaxy(1)*ysrc(1,it)+rcaxy(2)*ysrc(2,it)
        ymom(2,it)=rcaxy(3)*ysrc(1,it)+rcaxy(4)*ysrc(2,it)
        xmom(1,it)=rcaxy(1)*xsrc(1,it)+rcaxy(2)*xsrc(2,it)
        xmom(2,it)=rcaxy(3)*xsrc(1,it)+rcaxy(4)*xsrc(2,it)
!#endif
        aflx(1,it)=rcaxy(1)*asrc(1,it)+rcaxy(2)*asrc(2,it)
        aflx(2,it)=rcaxy(3)*asrc(1,it)+rcaxy(4)*asrc(2,it)
      enddo

      do ig=1,2
        do it=1,myn
! 2013_06_13 . scb          
          !err=(oflx(ig,it)-aflx(ig,it))/aflx(ig,it)
          err=(oflx(ig,it)-aflx(ig,it))/aflx(ig,it)
          if(err.lt.0.) err=-err
          if(err.gt.errh) errh=err 
        enddo
      enddo

      do ig=1,2
        hflx(ig)=0
        do it=1,6 
          hflx(ig)=hflx(ig)+aflx(ig,it)
        enddo
        hflx(ig)=hflx(ig)/6.
      enddo

      do it=1,6
         cnto(1,it)=-cnti(1,it)+1./2.*cnto(1,it)
         cnto(2,it)=-cnti(2,it)+1./2.*cnto(2,it)
      enddo
!#endif

!#ifdef DEBUG      
!      do it=1,myn 
!        asrc(1,it)=srcbal(1,it)
!        xsrc(1,it)=srcmomx(1,it)
!        ysrc(1,it)=srcmomy(1,it)
!        asrc(2,it)=srcbal(2,it)
!        xsrc(2,it)=srcmomx(2,it)
!        ysrc(2,it)=srcmomy(2,it)
!      enddo
!
!      caxy(1)=sigr(1)-signf(1)*reigv
!      caxy(2)=-signf(2)*reigv
!      caxy(3)=-sigs
!      caxy(4)=sigr(2)
!      do ig=1,2
!        bt=dc(ig)/sz/sz  
!        caxy(3*ig-2)=caxy(3*ig-2)+80.*bt
!        ca(ig,1)=-32.*bt
!        ca(ig,2)=-64.*bt
!        ca(ig,3)=16.*bt/3.
!        cx(ig,1)=8.*bt/3.
!        cx(ig,2)=-32.*bt/3.
!        cx(ig,3)=16.*bt/9.
!        cy(ig)=-8.*bt
!        cs(ig,1)=20. 
!        cs(ig,2)=-60.
!        cs(ig,3)=60.
!        cs(ig,4)=-48.
!        cs(ig,5)=2.
!        tt=dc(ig)/rt3/sz
!        cj(ig,1)=20.*tt
!        cj(ig,2)=120.*tt
!        cj(ig,3)=-1.-48.*tt
!        cj(ig,4)=2.*tt
!        cp(ig,1)=-15.
!        cp(ig,2)=2.
!        cp(ig,3)=-6.
!
!        psrc(ig)=0.
!        do it=1,mxo
!          asrc(ig,it)=asrc(ig,it)-ca(ig,2)*cnti(ig,it) &
!                      -ca(ig,3)*(pflx(ig,it)+pflx(ig,mp1(it)))
!          xsrc(ig,it)=xsrc(ig,it)-cx(ig,2)*cnti(ig,it) &
!                      +0.5*cx(ig,3)*(pflx(ig,it)+pflx(ig,mp1(it)))
!          ysrc(ig,it)=ysrc(ig,it) &
!                      +8.*bt/3.*(pflx(ig,it)-pflx(ig,mp1(it)))
!          ssrc(ig,it)=-2.*(pflx(ig,it) &
!                      +pflx(ig,mp1(it))+pflx(ig,mp5(it)))
!          ojsrc(ig,it)=(-1.+48.*tt)*cnti(ig,it) &
!                      -tt*(pflx(ig,it)+pflx(ig,mp1(it)))
!          psrc(ig)=psrc(ig)-2.*cnti(ig,it)
!        enddo
!      enddo
!!      
!!#ifdef DEBUG
!!918   format(a15,1p,e15.5)
!!919   format(a15,1p,2e15.5)      
!!1002  format(a15,1p,4e15.5)
!!1003  format(a15,1p,2e15.5)
!!      
!!      write(141002, 1002) 'a1 matrix', caxy(1:4)
!!      write(141002, 1003) 'a2 matrix', ca(1:2,1)
!!      write(141002, 1003) 'a3 matrix', ca(1:2,1)
!!      write(141002, 1003) 'a4 matrix', ca(1:2,2)
!!      write(141002, 1003) 'a5 matrix', ca(1:2,3)
!!      write(14
!!      
!!      
!!#endif      
!!
!      errh=0.
!      do it=1,myn
!        do ig=1,2
!          oflx(ig,it)=aflx(ig,it)
!        enddo
!      enddo
!
!      det=1./(caxy(1)*caxy(4)-caxy(2)*caxy(3))
!      rcaxy(1)=det*caxy(4)
!      rcaxy(2)=-det*caxy(2)
!      rcaxy(3)=-det*caxy(3)
!      rcaxy(4)=det*caxy(1)
!
!      do it=1,3
!        do ig=1,4
!          csp(ig,it)=0
!          cjp(ig,it)=0
!          cpp(ig,it)=0
!        enddo
!      enddo
!      do ig=1,4
!        csp(ig,4)=0
!      enddo
!      do ig=1,4,3
!        csp(ig,1)=cs(mg1(ig),4)
!        csp(ig,4)=cs(mg1(ig),5)
!        cjp(ig,2)=cj(mg1(ig),3)
!        cjp(ig,3)=cj(mg1(ig),4)
!        cpp(ig,2)=cp(mg1(ig),2)
!        cpp(ig,3)=cp(mg1(ig),3)
!      enddo
!
!      do ig=1,4
!        fac1(ig)=-rcaxy(ig)*cs(mg1(ig),1)
!        fac2(ig)=-rcaxy(ig)*cs(mg1(ig),2)
!        fac3(ig)=-rcaxy(ig)*cs(mg1(ig),3)
!      enddo
!      do ig=1,4
!        sfax=ca(mg2(ig),1)*fac1(ig)+cx(mg2(ig),1)*fac2(ig)
!        sfy=cy(mg2(ig))*fac3(ig)
!        csp(ig,1)=csp(ig,1)+2.*(sfax+sfy)
!        csp(ig,2)=csp(ig,2)+sfax-sfy
!        csp(ig,3)=csp(ig,3)+ca(mg2(ig),2)*fac1(ig) &
!                      +cx(mg2(ig),2)*fac2(ig)
!        csp(ig,4)=csp(ig,4)+2.*(ca(mg2(ig),3)*fac1(ig) &
!                      +cx(mg2(ig),3)*fac2(ig))
!      enddo
!      do it=1,myn
!        sfa1=asrc(1,it)+asrc(1,mp5(it))
!        sfa2=asrc(2,it)+asrc(2,mp5(it))
!        sfx1=xsrc(1,it)+xsrc(1,mp5(it))
!        sfx2=xsrc(2,it)+xsrc(2,mp5(it))
!        sfy1=ysrc(1,it)-ysrc(1,mp5(it))
!        sfy2=ysrc(2,it)-ysrc(2,mp5(it))
!        ssrc(1,it)=ssrc(1,it)+fac1(1)*sfa1+fac1(2)*sfa2+fac2(1)*sfx1 &
!                      +fac2(2)*sfx2+fac3(1)*sfy1+fac3(2)*sfy2
!        ssrc(2,it)=ssrc(2,it)+fac1(3)*sfa1+fac1(4)*sfa2+fac2(3)*sfx1 &
!                      +fac2(4)*sfx2+fac3(3)*sfy1+fac3(4)*sfy2
!      enddo
!
!      do ig=1,4
!        fac1(ig)=-rcaxy(ig)*cj(mg1(ig),1)
!        fac2(ig)=-rcaxy(ig)*cj(mg1(ig),2)
!      enddo
!      do ig=1,4
!       cjp(ig,1)=cjp(ig,1)+ca(mg2(ig),1)*fac1(ig)+cx(mg2(ig),1)*fac2(ig)
!       cjp(ig,2)=cjp(ig,2)+ca(mg2(ig),2)*fac1(ig)+cx(mg2(ig),2)*fac2(ig)
!       cjp(ig,3)=cjp(ig,3)+ca(mg2(ig),3)*fac1(ig)+cx(mg2(ig),3)*fac2(ig)
!      enddo
!      do it=1,myn
!        ojsrc(1,it)=ojsrc(1,it)+fac1(1)*asrc(1,it)+fac2(1)*xsrc(1,it) &
!                      +fac1(2)*asrc(2,it)+fac2(2)*xsrc(2,it)
!        ojsrc(2,it)=ojsrc(2,it)+fac1(3)*asrc(1,it)+fac2(3)*xsrc(1,it) &
!                      +fac1(4)*asrc(2,it)+fac2(4)*xsrc(2,it)
!      enddo
!      do ig=1,4
!        fac(ig)=-rcaxy(ig)*cp(mg1(ig),1)
!      enddo
!      do ig=1,4
!        cpp(ig,1)=cpp(ig,1)+2.*cx(mg2(ig),1)*fac(ig)
!        cpp(ig,2)=cpp(ig,2)+cx(mg2(ig),2)*fac(ig)
!        cpp(ig,3)=cpp(ig,3)+6.*cx(mg2(ig),3)*fac(ig)
!      enddo
!      do ig=1,2
!        ttt(ig)=0.
!        do it=1,myn
!          ttt(ig)=ttt(ig)+xsrc(ig,it)
!        enddo
!      enddo
!      psrc(1)=psrc(1)+fac(1)*ttt(1)+fac(2)*ttt(2)
!      psrc(2)=psrc(2)+fac(3)*ttt(1)+fac(4)*ttt(2)
!
!      det=1./(cjp(1,2)*cjp(4,2)-cjp(2,2)*cjp(3,2))
!      rcjp(1)=det*cjp(4,2)
!      rcjp(2)=-det*cjp(2,2)
!      rcjp(3)=-det*cjp(3,2)
!      rcjp(4)=det*cjp(1,2)
!
!      do ig=1,4
!        fac(ig)=-csp(mg3(ig),3)*rcjp(mg2(ig)) &
!                      -csp(mg4(ig),3)*rcjp(mg5(ig))
!      enddo
!      do ig=1,4
!        dat1=fac(mg3(ig))*cjp(mg2(ig),1)+fac(mg4(ig))*cjp(mg5(ig),1)
!        dat2=fac(mg3(ig))*cjp(mg2(ig),3)+fac(mg4(ig))*cjp(mg5(ig),3)
!        csp(ig,1)=csp(ig,1)+2.*dat1
!        csp(ig,2)=csp(ig,2)+dat1
!        csp(ig,4)=csp(ig,4)+2.*dat2
!      enddo
!      do it=1,myn
!        dat1=ojsrc(1,it)+ojsrc(1,mp5(it))
!        dat2=ojsrc(2,it)+ojsrc(2,mp5(it))
!        ssrc(1,it)=ssrc(1,it)+fac(1)*dat1+fac(2)*dat2
!        ssrc(2,it)=ssrc(2,it)+fac(3)*dat1+fac(4)*dat2
!      enddo
!      do ig=1,4
!        fac(ig)=-cpp(mg3(ig),2)*rcjp(mg2(ig)) &
!                      -cpp(mg4(ig),2)*rcjp(mg5(ig))
!      enddo
!      do ig=1,4
!        cpp(ig,1)=cpp(ig,1)+2.*fac(mg3(ig))*cjp(mg2(ig),1) &
!                      +2.*fac(mg4(ig))*cjp(mg5(ig),1)
!        cpp(ig,3)=cpp(ig,3)+6.*fac(mg3(ig))*cjp(mg2(ig),3) &
!                      +6.*fac(mg4(ig))*cjp(mg5(ig),3)
!      enddo
!      do ig=1,2
!        ttt(ig)=0.
!        do it=1,myn
!          ttt(ig)=ttt(ig)+ojsrc(ig,it)
!        enddo
!      enddo
!      psrc(1)=psrc(1)+fac(1)*ttt(1)+fac(2)*ttt(2)
!      psrc(2)=psrc(2)+fac(3)*ttt(1)+fac(4)*ttt(2)
!
!      do ig=1,4
!        fac1(ig)=csp(ig,1)+2.*csp(ig,2)
!        fac2(ig)=6.*csp(ig,4)
!      enddo
!      rdet=1./(fac1(1)*fac1(4)-fac1(2)*fac1(3))
!      rcsp(1)=rdet*fac1(4)
!      rcsp(2)=-rdet*fac1(2)
!      rcsp(3)=-rdet*fac1(3)
!      rcsp(4)=rdet*fac1(1)
!      do ig=1,4
!        fac(ig)=-cpp(mg3(ig),1)*rcsp(mg2(ig)) &
!                      -cpp(mg4(ig),1)*rcsp(mg5(ig))
!      enddo
!      do ig=1,4
!        cpp(ig,3)=cpp(ig,3)+fac(mg3(ig))*fac2(mg2(ig)) &
!                      +fac(mg4(ig))*fac2(mg5(ig))
!      enddo
!      do ig=1,2
!        ttt(ig)=0.
!        do it=1,myn
!          ttt(ig)=ttt(ig)+ssrc(ig,it)
!        enddo
!      enddo
!      psrc(1)=psrc(1)+fac(1)*ttt(1)+fac(2)*ttt(2)
!      psrc(2)=psrc(2)+fac(3)*ttt(1)+fac(4)*ttt(2)
!      rdet=1./(cpp(1,3)*cpp(4,3)-cpp(2,3)*cpp(3,3))
!      cpfl(1)=(cpp(4,3)*psrc(1)-cpp(2,3)*psrc(2))*rdet
!      cpfl(2)=(-cpp(3,3)*psrc(1)+cpp(1,3)*psrc(2))*rdet
!
!      do it=1,myn
!        ssrc(1,it)=ssrc(1,it)-csp(1,4)*cpfl(1)-csp(2,4)*cpfl(2)
!        ssrc(2,it)=ssrc(2,it)-csp(3,4)*cpfl(1)-csp(4,4)*cpfl(2)
!      enddo
!
!      det=1./(csp(1,1)*csp(4,1)-csp(2,1)*csp(3,1))
!      rcsp(1)=det*csp(4,1)
!      rcsp(2)=-det*csp(2,1)
!      rcsp(3)=-det*csp(3,1)
!      rcsp(4)=det*csp(1,1)
!      a5=rcsp(1)*csp(1,2)+rcsp(2)*csp(3,2)
!      a6=rcsp(1)*csp(2,2)+rcsp(2)*csp(4,2)
!      a7=rcsp(3)*csp(1,2)+rcsp(4)*csp(3,2)
!      a8=rcsp(3)*csp(2,2)+rcsp(4)*csp(4,2)
!      do it=1,myn
!        ttt(1)=ssrc(1,it)
!        ttt(2)=ssrc(2,it)
!        ssrc(1,it)=rcsp(1)*ttt(1)+rcsp(2)*ttt(2)
!        ssrc(2,it)=rcsp(3)*ttt(1)+rcsp(4)*ttt(2)
!      enddo
!      d1=1./6./((1.-a8)*(1.-a5)-a6*a7)
!      d2=1./6./((1.+a8)*(1.+a5)-a6*a7)
!      d3=1./6./((1.-2.*a8)*(1.-2.*a5)-4.*a6*a7)
!      d4=1./6./((1.+2.*a8)*(1.+2.*a5)-4.*a6*a7)
!      c1=2.*(d1+d2)+d3+d4+2.*a8*(-d1+d2-d3+d4)
!      c2=2.*a6*(d1-d2+d3-d4)
!      c3=-d1+d2-d3+d4+a8*(d1+d2+2.*(d3+d4))
!      c4=-a6*(d1+d2+2.*(d3+d4))
!      c5=-d1-d2+d3+d4+a8*(d1-d2-2.*(d3-d4))
!      c6=a6*(-d1+d2+2.*(d3-d4))
!      c7=2.*(d1-d2)-d3+d4+2.*a8*(-d1-d2+d3+d4)
!      c8=2.*a6*(d1+d2-d3-d4)
!      e1=2.*a7*(d1-d2+d3-d4)
!      e2=2.*(d1+d2)+d3+d4+2.*a5*(-d1+d2-d3+d4)
!      e3=-a7*(d1+d2+2.*(d3+d4))
!      e4=-d1+d2-d3+d4+a5*(d1+d2+2.*(d3+d4))
!      e5=a7*(-d1+d2+2.*(d3-d4))
!      e6=-d1-d2+d3+d4+a5*(d1-d2-2.*(d3-d4))
!      e7=2.*a7*(d1+d2-d3-d4)
!      e8=2.*(d1-d2)-d3+d4+2.*a5*(-d1-d2+d3+d4)
!      do it=1,myn
!        addsrc1=ssrc(1,mp1(it))+ssrc(1,mp5(it))
!        addsrc2=ssrc(2,mp1(it))+ssrc(2,mp5(it))
!        addsrc3=ssrc(1,mp2(it))+ssrc(1,mp4(it))
!        addsrc4=ssrc(2,mp2(it))+ssrc(2,mp4(it))
!        dat1=ssrc(1,it)
!        dat2=ssrc(2,it)
!        dat3=ssrc(1,mp3(it))
!        dat4=ssrc(2,mp3(it))
!        sfli(1,it)=c1*dat1+c2*dat2+c3*addsrc1+c4*addsrc2 &
!                      +c5*addsrc3+c6*addsrc4+c7*dat3+c8*dat4
!        sfli(2,it)=e1*dat1+e2*dat2+e3*addsrc1+e4*addsrc2 &
!                      +e5*addsrc3+e6*addsrc4+e7*dat3+e8*dat4
!      enddo
!
!      do it=1,myn
!        addsfl1=sfli(1,it)+sfli(1,mp1(it))
!        addsfl2=sfli(2,it)+sfli(2,mp1(it))
!        ojsrc(1,it)=ojsrc(1,it)-cjp(1,3)*cpfl(1)-cjp(2,3)*cpfl(2) &
!                      -cjp(1,1)*addsfl1-cjp(2,1)*addsfl2
!        ojsrc(2,it)=ojsrc(2,it)-cjp(3,3)*cpfl(1)-cjp(4,3)*cpfl(2) &
!                      -cjp(3,1)*addsfl1-cjp(4,1)*addsfl2
!        cnto(1,it)=rcjp(1)*ojsrc(1,it)+rcjp(2)*ojsrc(2,it)
!        cnto(2,it)=rcjp(3)*ojsrc(1,it)+rcjp(4)*ojsrc(2,it)
!      enddo
!
!      do it=1,myn
!        do ig=1,2
!          addsfl=sfli(ig,it)+sfli(ig,mp1(it))
!          subsfl=sfli(ig,it)-sfli(ig,mp1(it))
!          ocur=cnto(ig,it)
!          ysrc(ig,it)=ysrc(ig,it)-cy(ig)*subsfl
!          xsrc(ig,it)=xsrc(ig,it)-cx(ig,1)*addsfl &
!                      -cx(ig,2)*ocur-cx(ig,3)*cpfl(ig)
!          asrc(ig,it)=asrc(ig,it)-ca(ig,1)*addsfl &
!                      -ca(ig,2)*ocur-ca(ig,3)*cpfl(ig)
!        enddo
!        ymom(1,it)=rcaxy(1)*ysrc(1,it)+rcaxy(2)*ysrc(2,it)
!        ymom(2,it)=rcaxy(3)*ysrc(1,it)+rcaxy(4)*ysrc(2,it)
!        xmom(1,it)=rcaxy(1)*xsrc(1,it)+rcaxy(2)*xsrc(2,it)
!        xmom(2,it)=rcaxy(3)*xsrc(1,it)+rcaxy(4)*xsrc(2,it)
!        aflx(1,it)=rcaxy(1)*asrc(1,it)+rcaxy(2)*asrc(2,it)
!        aflx(2,it)=rcaxy(3)*asrc(1,it)+rcaxy(4)*asrc(2,it)
!      enddo
!
!      do ig=1,2
!        do it=1,myn
!          err=(oflx(ig,it)-aflx(ig,it))/aflx(ig,it)
!          if(err.lt.0.) err=-err
!          if(err.gt.errh) errh=err 
!        enddo
!      enddo
!
!      do ig=1,2
!        hflx(ig)=0
!        do it=1,6 
!          hflx(ig)=hflx(ig)+aflx(ig,it)
!        enddo
!        hflx(ig)=hflx(ig)/6.
!      enddo
!#endif      
      
#ifdef DEBUG
      write(140918, '(a)')  "========================= Output ========================="
      write(140918, 919) 'avg_phi_1', aflx(:,1)
      write(140918, 919) 'avg_phi_2', aflx(:,2)
      write(140918, 919) 'avg_phi_3', aflx(:,3)
      write(140918, 919) 'avg_phi_4', aflx(:,4)
      write(140918, 919) 'avg_phi_5', aflx(:,5)
      write(140918, 919) 'avg_phi_6', aflx(:,6)
      
      write(140918, 919) 'x_mom_1', xmom(:,1)
      write(140918, 919) 'x_mom_2', xmom(:,2)
      write(140918, 919) 'x_mom_3', xmom(:,3)
      write(140918, 919) 'x_mom_4', xmom(:,4)
      write(140918, 919) 'x_mom_5', xmom(:,5)
      write(140918, 919) 'x_mom_6', xmom(:,6)
      
      write(140918, 919) 'y_mom_1', ymom(:,1)
      write(140918, 919) 'y_mom_2', ymom(:,2)
      write(140918, 919) 'y_mom_3', ymom(:,3)
      write(140918, 919) 'y_mom_4', ymom(:,4)
      write(140918, 919) 'y_mom_5', ymom(:,5)
      write(140918, 919) 'y_mom_6', ymom(:,6)
      
      write(140918, 919) 'surf_phi_1', sfli(:,1)
      write(140918, 919) 'surf_phi_2', sfli(:,2)
      write(140918, 919) 'surf_phi_3', sfli(:,3)
      write(140918, 919) 'surf_phi_4', sfli(:,4)
      write(140918, 919) 'surf_phi_5', sfli(:,5)
      write(140918, 919) 'surf_phi_6', sfli(:,6)
      
      write(140918, 919) 'J_out_1', cnto(:,1)
      write(140918, 919) 'J_out_2', cnto(:,2)
      write(140918, 919) 'J_out_3', cnto(:,3)
      write(140918, 919) 'J_out_4', cnto(:,4)
      write(140918, 919) 'J_out_5', cnto(:,5)
      write(140918, 919) 'J_out_6', cnto(:,6)
            
      write(140918, 919) 'cent_phi', cpfl(:)

      close(140918)
      !stop
#endif      
      return
      end

