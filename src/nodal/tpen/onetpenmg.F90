! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
      subroutine onetpenmg(ng,nassy,nz,l,k,jupstt, &
           iscatib,iscatie,                        &
           xst,xsnf,xss,xsd,chi,sz,hz,hzu,hzl,     &
           pflbt,cntit,cntz0i,cntz1i,seff,srczn,   &
           dtlu,dtll,rxkeff,                       &
           aflxt,xmomt,ymomt,zmom1,zmom2,          &
           cntot,cntz0o,cntz1o,hflx) ! output

      use sfam_cntl, only : epsl2      ! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
      
			implicit double precision (a-h,o-z)

      dimension xst(ng),xsnf(ng),xsd(ng),chi(ng),xss(ng,ng)
      dimension pflbt(ng,6),cntit(ng,6),cntz0i(ng)
      dimension cntz1i(ng),seff(ng),srczn(ng,6),dtlu(ng),dtll(ng)
      dimension cntot(ng,6),cntz0o(ng),cntz1o(ng)
      dimension hflx(ng),aflxt(ng,6),xmomt(ng,6),ymomt(ng,6)
      dimension zmom1(ng),zmom2(ng)

      dimension asrc(6),xsrc(6),ysrc(6),ssrcb(6),ssrce(6),ssrcold(6), &
                sflb(6),sfle(6), &
                tdt1(6),tdt2(6),tdt3(6)

      real,allocatable,save,dimension(:) :: tem1,tem2,hflxold, &
                raxy,ca1,csfi2,rptj,cbd1,cbd2,cbd3,rbd,rsfsm,  &
                aax,czj1,czj2,ctl1i,cmom1,cmom2,ccpt1,ccpt2,   &
                pflcbd,sj1bd,sj0bd,csj,smm1bd,smm2bd,csmm1,csmm2, &
                csmm2hfl,csmm2sjs,csmm1sj,xax,scgr,scgrb, &
                cntzisum,cntisum,pflbsum
      real,allocatable,save,dimension(:,:) :: asrcbd,xsrcbd,ysrcbd, &
                aflxbd,xmombd,ymombd,cntim,pflbm,ccsf,sflebd,ssrcbbd, &
                ccgr,rcgr
      dimension iscatib(ng),iscatie(ng)
      dimension pflb(6),cnti(6),cnto(6),aflx(6),xmom(6),ymom(6)

      logical if2d,if1d,if3d

      save onenine,onesix,one12,oneten,onethree,rsz,rsz2,areavol, &
           rfac1,rfac2,rfac3,if2d,if1d,if3d,rt3

      dimension mp1(6),mp2(6),mp3(6),mp4(6),mp5(6)
      data mp1/2,3,4,5,6,1/ 
      data mp2/3,4,5,6,1,2/ 
      data mp3/4,5,6,1,2,3/ 
      data mp4/5,6,1,2,3,4/ 
      data mp5/6,1,2,3,4,5/

	  logical,save :: first=.TRUE.

      if(first) then
        onethree=1/3.0D0
        onenine=onethree*onethree
        onesix=0.5*onethree
        one12=0.5*onesix
        oneten=0.1D0
        rsz=1./sz
        rsz2=rsz*rsz
        rfac1=1/540.0D0
        rfac2=1/3240.0D0
        rfac3=1/360.0D0
	    rt3=sqrt(3.)
        areavol=2.*rt3*onenine*rsz

        if1d=.false.
        if2d=.false.
        if3d=.false.
        if(nassy.eq.1) if1d=.true.
        if(nz.eq.1) if2d=.true.
        if(.not.if1d.and..not.if2d) if3d=.true.
        allocate(tem1(ng),tem2(ng),hflxold(ng), &
                 asrcbd(ng,6),xsrcbd(ng,6),ysrcbd(ng,6), &
                 raxy(ng),ca1(ng),csfi2(ng),rptj(ng), &
                 cbd1(ng),cbd2(ng),cbd3(ng),rbd(ng),rsfsm(ng), &
                 aax(ng),czj1(ng),czj2(ng), &
                 ctl1i(ng),cmom1(ng),cmom2(ng), &
                 ccpt1(ng),ccpt2(ng),ccsf(9,ng), &
                 pflcbd(ng),sflebd(6,ng),ssrcbbd(6,ng), &
                 aflxbd(6,ng),xmombd(6,ng),ymombd(6,ng), &
                 sj1bd(ng),sj0bd(ng),csj(ng), &
                 smm1bd(ng),smm2bd(ng),csmm1(ng),csmm2(ng), &
                 csmm2hfl(ng),csmm2sjs(ng),csmm1sj(ng),xax(ng), &
                 cntim(6,ng),pflbm(6,ng), &
                 scgr(ng-jupstt+1),ccgr(ng-jupstt+1,ng-jupstt+1), &
                 rcgr(ng-jupstt+1,ng-jupstt+1),scgrb(ng-jupstt+1), &
                 cntzisum(ng),cntisum(ng),pflbsum(ng))
        tem1=0
        tem2=0
        hflxold=0
 
        raxy=0
        ca1=0
        csfi2=0
        rptj=0
        cbd1=0
        cbd2=0
        cbd3=0
        rbd=0
        rsfsm=0
        aax=0
        czj1=0
        czj2=0
        ctl1i=0
        cmom1=0
        cmom2=0
        ccpt1=0
        ccpt2=0
        pflcbd=0
        sj1bd=0
        sj0bd=0
        csj=0
        smm1bd=0
        smm2bd=0
        csmm1=0
        csmm2=0
        csmm2hfl=0
        csmm2sjs=0
        csmm1sj=0
        xax=0
        scgr=0
        scgrb=0
        cntzisum=0
        cntisum=0
        pflbsum=0
        asrcbd=0
        xsrcbd=0
        ysrcbd=0         
        aflxbd=0
        xmombd=0
        ymombd=0
        cntim=0
        pflbm=0
        ccsf=0
        sflebd=0
        ssrcbbd=0
        ccgr=0
        rcgr=0
        first=.FALSE.
      endif

      rhz=1./hz
      rhz2=rhz*rhz
      rhzl=1./hzl
      rhzu=1./hzu
      faccsf=-190.*onenine*rhz
      faccpt=-5.*onethree*rhz
      faccbd=-80.*onenine*rhz
      facmx=rhz/54.

      do ii=1,6
        do j=1,ng
          cntim(ii,j)=cntit(j,ii)
          pflbm(ii,j)=pflbt(j,ii)
          asrcbd(j,ii)=0
          xsrcbd(j,ii)=0
          ysrcbd(j,ii)=0
        enddo
      enddo
      do j=1,ng
        hflxold(j)=hflx(j)
        tem1(j)=seff(j)
        tem2(j)=seff(j)
        cntzisum(j)=cntz0i(j)+cntz1i(j)
        cntisum(j)=cntim(1,j)+cntim(2,j)+cntim(3,j)+cntim(4,j) &
                  +cntim(5,j)+cntim(6,j)
        pflbsum(j)=pflbm(1,j)+pflbm(2,j)+pflbm(3,j)+pflbm(4,j) &
                  +pflbm(5,j)+pflbm(6,j)
      enddo

      do j=1,ng
        if(.not.if1d) then
          bt=xsd(j)*rsz2
          gam=rt3*sz/xsd(j)
          raxy(j)=1./(xst(j)+80.*bt)
          ca1(j)=32.*bt*raxy(j)
          rsf=1./(40.*ca1(j)-24.)
          csf1=5.*ca1(j)*rsf
          rpt=1./(5.*ca1(j)-6.)
          cpt1=-5.*onesix*ca1(j)*rpt
          cpt3=-3.*cpt1
          cpt2=-cpt3+rpt
          rbd(j)=4./(80.*ca1(j)-48.-gam)
          cbd1(j)=5.*ca1(j)*rbd(j)
          cbd2(j)=rbd(j)-cbd1(j)
          do ii=1,6
            ssrcbbd(ii,j)=rbd(j)*(-0.5*(pflbm(ii,j)+pflbm(mp1(ii),j)) &
                         -gam*cntim(ii,j))
          enddo
          ttt=onesix*ca1(j)
          do ii=1,6
             tdt3(ii)=ttt*pflbm(ii,j)
          enddo
          do ii=1,6
            i1=mp1(ii)
            aflxbd(ii,j)=-tdt3(ii)-tdt3(i1)
            xmombd(ii,j)=onesix*(tdt3(ii)+tdt3(i1))
            ymombd(ii,j)=0.5*(tdt3(ii)-tdt3(i1))
          enddo
!   delete the boundary surfaces
          fac1=csf1*cbd1(j)
          rsfi=1./(1.-2.*fac1)
          csfi1=rsfi*(csf1-fac1)
          csfi2(j)=rsfi*(rsf-2.*csf1*cbd2(j))
          rpti=1./(1.-6.*cpt2*cbd2(j))
          cpti1=rpti*(cpt3-2.*cpt2*cbd1(j))
!   delele the inner surfaces
          fac1=onesix/(-1.+csfi1)
          fac2=onesix/(1.+csfi1)
          fac3=onesix/(-1.+2.*csfi1)
          fac4=onesix/(1.+2.*csfi1)
          rsi1=-2.*(fac1-fac2)-fac3+fac4
          rsi2=fac1+fac2+fac3+fac4
          rsi3=fac1-fac2-fac3+fac4
          rsi4=-2.*(fac1+fac2)+fac3+fac4
          rsfsm(j)=rsi1+2.*(rsi2+rsi3)+rsi4
          rptj(j)=1./(1.-6.*cpti1*rsfsm(j)*csfi2(j))

          t1=cpti1*rsfsm(j)*rsfi
          t2=rbd(j)*(rpti*cpt2-2.*t1*csf1)
          t3=t1*rsf
          ccpt1(j)=10.*(t2+2.*t3)
          ccpt2(j)=15.*(rpti*rpt+4.*(t2-t3))
          ccpt3=-rpti*cpt1-t1*0.25-3.*t3+t2
          ccpt4=gam*t2
          pflcbd(j)=ccpt3*pflbsum(j)+ccpt4*cntisum(j)

          t1=rsfi*rsi1
          t2=rsfi*rsi2
          t3=rsfi*rsi3
          t4=rsfi*rsi4
          tt=csf1*rbd(j)
          q4=-rsf+tt
          q1=10.*q4
          ccsf(1,j)=(t1+t2)*q1
          ccsf(2,j)=(t2+t3)*q1
          ccsf(3,j)=(t3+t4)*q1
          q2=60.*tt
          q3=-30.*rsf
          q7=q2-q3
          ccsf(4,j)=(t1+t2)*q7
          ccsf(5,j)=(t2+t3)*q7
          ccsf(6,j)=(t3+t4)*q7
          ccsf(7,j)=(t1-t2)*q3
          ccsf(8,j)=(t2-t3)*q3
          ccsf(9,j)=(t3-t4)*q3
          q5=2.*csf1-rsf+q4
          ccsf11=q4*t1+q5*t2
          ccsf12=0.5*q5*(t1+t3)+q4*t2
          ccsf13=0.5*q5*(t2+t4)+q4*t3
          ccsf14=q5*t3+q4*t4
          tt=tt*gam
          ccsf15=tt*(t1+t2)
          ccsf16=tt*(t2+t3)
          ccsf17=tt*(t3+t4)

          do ii=1,3
            datj1=cntim(ii,j)+cntim(mp5(ii),j)
            datj2=ccsf16*(cntim(mp1(ii),j)+cntim(mp4(ii),j))
            datj3=cntim(mp2(ii),j)+cntim(mp3(ii),j)
            datc1=pflbm(mp1(ii),j)+pflbm(mp5(ii),j)
            datc2=pflbm(mp2(ii),j)+pflbm(mp4(ii),j)
            sflebd(ii,j)=ccsf11*pflbm(ii,j)+ccsf12*datc1 &
                        +ccsf13*datc2+ccsf14*pflbm(mp3(ii),j) &
                        +ccsf15*datj1+datj2+ccsf17*datj3
            sflebd(mp3(ii),j)=ccsf11*pflbm(mp3(ii),j)+ccsf12*datc2 &
                        +ccsf13*datc1+ccsf14*pflbm(ii,j) &
                        +ccsf15*datj3+datj2+ccsf17*datj1
          enddo
        endif

        if(if3d) then
          csf4=faccsf*raxy(j)*rsf
          cbd3(j)=faccbd*raxy(j)*rbd(j)
          csfi3=rsfi*(csf4-2.*csf1*cbd3(j))
          cpt4=rpti*(faccpt*raxy(j)*rpt-6.*cpt2*cbd3(j))
          cpt4=rptj(j)*(cpt4-6.*cpti1*rsfsm(j)*csfi3)
          csf4=rsfsm(j)*(csfi3-csfi2(j)*cpt4)
          cbd3(j)=cbd3(j)-cbd2(j)*cpt4-2.*cbd1(j)*csf4
          aax(j)=rhz*raxy(j)+ca1(j)*(2.*csf4+cbd3(j)-onesix*cpt4)
          aax(j)=0.5*aax(j)
          xax(j)=-onesix*ca1(j)*(cbd3(j)-csf4-onethree*cpt4) &
                 +facmx*raxy(j)
        endif

        if(if1d) then
          aax(j)=0.5*rhz/xst(j)
        endif
        if(.not.if2d) then
          gamz=hz/xsd(j)
          fac1=10.*aax(j)+8.+0.25*gamz
          fac2=10.*aax(j)+2.
          rdet=1./(fac1*fac1-fac2*fac2)
          rzj1=rdet*fac1
          rzj2=rdet*fac2
          sj1bd(j)=gamz*(rzj1*cntz1i(j)-rzj2*cntz0i(j))
          sj0bd(j)=gamz*(rzj1*cntz0i(j)-rzj2*cntz1i(j))
          csj(j)=10.*(rzj1-rzj2)
          czj1(j)=rzj1+rzj2
          czj2(j)=rzj1-rzj2
        endif
        if(if1d) then
          ctl1i(j)=0
          ctl2i=0
          smm1bd(j)=0
          smm2bd(j)=0
        else
          fac2=105.*areavol*cbd3(j)*czj2(j)
          rbtmb=1./(rhzl+rhz)
          rbtmu=1./(rhzu+rhz)
          tl0bd=dtll(j)*rhzl*rbtmb
          tl1bd=dtlu(j)*rhzu*rbtmu
          smm1bd(j)=-onesix*(tl1bd-tl0bd)
          smm2bd(j)=oneten*(tl1bd+tl0bd)
          csmm1(j)=-onesix*rhz*(rbtmu-rbtmb)
          csmm2(j)=oneten*(-2.+rhz*(rbtmu+rbtmb))
          ctl1i(j)=-csmm1(j)*fac2
          ctl2i=-csmm2(j)*fac2
        endif
        if(.not.if2d) then
          btz=xsd(j)*rhz2
          dat=14.*btz*(1.+2.*aax(j))
          cmom2(j)=xst(j)+140.*btz-70.*dat*czj2(j)+ctl2i
          cmom1(j)=xst(j)+60.*btz*(1.-5.*czj1(j))
          csmm2hfl(j)=28.*btz
          csmm2sjs(j)=-14.*btz*(1.+2.*aax(j))
          csmm1sj(j)=10.*btz
        endif
      enddo

      sumfs=0.
      sumcnti=0
      sumcntzi=0
      do j=1,ng
        sumfs=sumfs+xsnf(j)*hflx(j)
      enddo
      sumfs=rxkeff*sumfs
      do j=1,ng-jupstt+1
        scgrb(j)=0
      enddo
      if(.not.if1d) then
        do j=1,ng
          tem1(j)=tem1(j)+rhz*cntzisum(j)
          sumcnti=sumcnti+cntisum(j)
        enddo
        do j1=jupstt,ng
          j=j1-jupstt+1
          scgrb(j)=scgrb(j)+areavol*cntisum(j1)
        enddo
        do ii=1,6
          tfsrc=0.
          tfsrcx=0.
          tfsrcy=0.
          do j=1,ng
            tfsrc=tfsrc+xsnf(j)*aflxt(j,ii)
            tfsrcx=tfsrcx+xsnf(j)*xmomt(j,ii)
            tfsrcy=tfsrcy+xsnf(j)*ymomt(j,ii)
          enddo
          tfsrc=rxkeff*tfsrc
          tfsrcx=rxkeff*tfsrcx
          tfsrcy=rxkeff*tfsrcy
          do j=1,ng
            ttt1=srczn(j,mp1(ii))+srczn(j,mp5(ii))
            ttt2=srczn(j,mp2(ii))+srczn(j,mp4(ii))
            ttt3=srczn(j,mp1(ii))-srczn(j,mp5(ii))
            ttt4=srczn(j,mp2(ii))-srczn(j,mp4(ii))
            asrcbd(j,ii)=chi(j)*tfsrc+tem1(j)+(83.*srczn(j,ii)+17.*ttt1 &
                        -37.*ttt2-43.*srczn(j,mp3(ii)))*rfac1 
            xsrcbd(j,ii)=chi(j)*tfsrcx+(-60.*tem1(j)+59.*srczn(j,ii) &
                        +14.*ttt1-10.*ttt2-7.*srczn(j,mp3(ii)))*rfac2 
            ysrcbd(j,ii)=chi(j)*tfsrcy-(9.*ttt3+ttt4)*rfac3
          enddo
        enddo
      endif

      if(.not.if2d) then
        do j=1,ng
          tem2(j)=tem2(j)+areavol*cntisum(j)
          sumcntzi=sumcntzi+cntzisum(j)
        enddo
        do j1=jupstt,ng
          j=j1-jupstt+1
          scgrb(j)=scgrb(j)+rhz*cntzisum(j1)
        enddo
        tfsrcz1=0.
        tfsrcz2=0.
        do j=1,ng
          tfsrcz1=tfsrcz1+xsnf(j)*zmom1(j)
          tfsrcz2=tfsrcz2+xsnf(j)*zmom2(j)
        enddo
        tfsrcz1=rxkeff*tfsrcz1
        tfsrcz2=rxkeff*tfsrcz2
        tfsrcz=sumfs
        do j=1,ng
          smm1bd(j)=smm1bd(j)+chi(j)*tfsrcz1
          smm2bd(j)=smm2bd(j)+chi(j)*tfsrcz2
        enddo
      endif

      nitrsrc=1
      do 1000 iter=1,nitrsrc
        if(iter.eq.1) then
          jstt=1
        else
          jstt=jupstt
        endif
        do 1100 j=jstt,ng
          do ii=1,6
            asrc(ii)=asrcbd(j,ii)
            xsrc(ii)=xsrcbd(j,ii)
            ysrc(ii)=ysrcbd(j,ii)
            do j2=iscatib(j),iscatie(j)
              asrc(ii)=asrc(ii)+xss(j2,j)*aflxt(j2,ii)
              xsrc(ii)=xsrc(ii)+xss(j2,j)*xmomt(j2,ii)
              ysrc(ii)=ysrc(ii)+xss(j2,j)*ymomt(j2,ii)
            enddo
            asrc(ii)=raxy(j)*asrc(ii)
            xsrc(ii)=raxy(j)*xsrc(ii)
            ysrc(ii)=raxy(j)*ysrc(ii)
          enddo

          pflc=rptj(j)*(pflcbd(j) &
            +ccpt1(j)*(asrc(1)+asrc(2)+asrc(3)+asrc(4)+asrc(5)+asrc(6)) &
            +ccpt2(j)*(xsrc(1)+xsrc(2)+xsrc(3)+xsrc(4)+xsrc(5)+xsrc(6)))
! calculate surface flux
          dat0=rsfsm(j)*csfi2(j)*pflc
          do ii=1,3
            i3=mp3(ii)
            ad1=asrc(ii)+asrc(mp5(ii))
            ad2=ccsf(2,j)*(asrc(mp1(ii))+asrc(mp4(ii)))
            ad3=asrc(mp2(ii))+asrc(i3)
            xd1=xsrc(ii)+xsrc(mp5(ii))
            xd2=ccsf(5,j)*(xsrc(mp1(ii))+xsrc(mp4(ii)))
            xd3=xsrc(mp2(ii))+xsrc(i3)
            yd1=ysrc(ii)-ysrc(mp5(ii))
            yd2=ccsf(8,j)*(ysrc(mp1(ii))-ysrc(mp4(ii)))
            yd3=ysrc(mp2(ii))-ysrc(i3)
            sfle(ii)=ccsf(1,j)*ad1+ad2+ccsf(3,j)*ad3 &
                    +ccsf(4,j)*xd1+xd2+ccsf(6,j)*xd3 &
                    +ccsf(7,j)*yd1+yd2+ccsf(9,j)*yd3-dat0+sflebd(ii,j)
            sfle(i3)=ccsf(1,j)*ad3+ad2+ccsf(3,j)*ad1 &
                    +ccsf(4,j)*xd3+xd2+ccsf(6,j)*xd1 &
                    -ccsf(7,j)*yd3-yd2-ccsf(9,j)*yd1-dat0+sflebd(i3,j)
          enddo
! calculate boundary surface flux
          dat0=cbd2(j)*pflc
          do ii=1,6
            sflb(ii)=ssrcbbd(ii,j)-10.*rbd(j)*(asrc(ii)+6.*xsrc(ii)) &
                    -cbd1(j)*(sfle(ii)+sfle(mp1(ii)))-dat0
            cnto(ii)=0.5*sflb(ii)-cntim(ii,j)
          enddo
!  calculate node average flux for each small triangle 
          hflxr=0
          do ii=1,6
            tdt1(ii)=ca1(j)*sfle(ii)
            tdt2(ii)=ca1(j)*sflb(ii)
          enddo
          dat1=onesix*ca1(j)*pflc
          do ii=1,6
            i1=mp1(ii)
            aflx(ii)=aflxbd(ii,j) &
                    +asrc(ii)+tdt1(ii)+tdt1(i1)+tdt2(ii)-dat1 
            xmom(ii)=xmombd(ii,j)+xsrc(ii) &
                    +one12*(2.*tdt2(ii)-tdt1(ii)-tdt1(i1)-4.*dat1) 
            ymom(ii)=ymombd(ii,j)+ysrc(ii)-0.25*(tdt1(i1)-tdt1(ii))
            hflxr=hflxr+aflx(ii)
          enddo
          hflxr=onesix*hflxr
          if(if2d) then
            hflx(j)=hflxr
            do ii=1,6
              aflxt(j,ii)=aflx(ii)
              xmomt(j,ii)=xmom(ii)
              ymomt(j,ii)=ymom(ii)
              cntot(j,ii)=cnto(ii)
            enddo
            goto 1100
          endif

!         calculate coefficients for axial nem                 
!            phi_h=aax*(J_b+J_t)+hsrcax                        
!            sum of phi_bd=csf4*(J_b+J_t)+bsrcax                

          dat=2.*aax(j)*cntzisum(j)
          hflxr=hflxr+dat
          do ii=1,6
            aflx(ii)=aflx(ii)+dat
          enddo

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  100     continue
          smm1=smm1bd(j)
          smm2=smm2bd(j)
          do j2=iscatib(j),iscatie(j)
            smm1=smm1+xss(j2,j)*zmom1(j2)
            smm2=smm2+xss(j2,j)*zmom2(j2)
          enddo
          if(if1d) then
            sbal=chi(j)*tfsrcz+seff(j)
            do j2=iscatib(j),iscatie(j)
              sbal=sbal+xss(j2,j)*hflx(j2)
            enddo
            hflxr=sbal/xst(j)+4.*aax(j)*cntzisum(j)
          endif
! calculation for boundary surface flux
          ttt=csj(j)*hflxr
          sj1=sj1bd(j)+ttt
          sj0=sj0bd(j)+ttt
          sjsum=sj0+sj1
! Transverse Leakage and coeff's at interface
          if(.not.if1d) then
            dtlm=areavol*(cnto(1)+cnto(2)+cnto(3) &
                         +cnto(4)+cnto(5)+cnto(6))-tem2(j) &
                -1.5*areavol*cbd3(j)*(sjsum-2.*cntzisum(j))
            smm1=smm1+csmm1(j)*dtlm
            smm2=smm2+csmm2(j)*dtlm
          endif
!  z-moments
          zmom2(j)=(smm2+csmm2sjs(j)*sjsum+csmm2hfl(j)*hflxr)/cmom2(j)
          zmom1(j)=(smm1+csmm1sj(j)*(sj1-sj0)-ctl1i(j)*zmom2(j))/cmom1(j)
!  calculate outgoing axial currents
          fac1=15.*czj1(j)*zmom1(j)
          fac2=35.*czj2(j)*zmom2(j)
          phit=fac1-fac2+sj1
          phib=-fac1-fac2+sj0
          cntz1o(j)=0.5*phit-cntz1i(j)
          cntz0o(j)=0.5*phib-cntz0i(j)
          flxzsum=phit+phib
          cntzosum=cntz1o(j)+cntz0o(j)
!  calculate hex-avg. flux
          hflx(j)=hflxr-aax(j)*flxzsum
          if(if1d) goto 1100
!  calculate outgoing radial currents
          fac1=0.5*cbd3(j)*cntzosum
          dat=aax(j)*flxzsum
          dat2=xax(j)*cntzosum
          do ii=1,6
            cntot(j,ii)=cnto(ii)-fac1
            aflxt(j,ii)=aflx(ii)-dat
            xmomt(j,ii)=xmom(ii)+dat2
            ymomt(j,ii)=ymom(ii)
          enddo
 1100   enddo
!   check residual
      sumcnto=0
      if(.not.if1d) then
        do it=1,6
          do j=1,ng
            sumcnto=sumcnto+cntot(j,it)
         enddo
        enddo
      endif
      sumcntzo=0
      sumscat=0
      sumsrc=0
      sumtot=0
      do j=1,ng
        sumcntzo=sumcntzo+cntz0o(j)+cntz1o(j)
        sumsrc=sumsrc+seff(j)
        sumtot=sumtot+xst(j)*hflx(j)
        do j2=iscatib(j),iscatie(j)
          sumscat=sumscat+xss(j2,j)*hflx(j2)
        enddo
      enddo
      sumrhs=sumsrc+sumfs+areavol*sumcnti+rhz*sumcntzi
      sumlhs=sumtot-sumscat+areavol*sumcnto+rhz*sumcntzo

      errimbal=abs((sumrhs-sumlhs)/sumtot)
      if(errimbal.lt.0.05*epsl2) goto 999

!  coarse group acceleration
      do j1=jupstt,ng
        j=j1-jupstt+1
        scgr(j)=scgrb(j)+seff(j1)+chi(j1)*sumfs
        do j2=iscatib(j1),jupstt-1
          scgr(j)=scgr(j)+xss(j2,j1)*hflx(j2)
        enddo
        do j2=jupstt,ng
          ccgr(j2-jupstt+1,j)=-xss(j2,j1)*hflx(j2)
        enddo
        ccgr(j,j)=xst(j1)*hflx(j1)+ccgr(j,j)
      enddo
      if(.not.if1d) then
        do j1=jupstt,ng
          j=j1-jupstt+1
          ccgr(j,j)=ccgr(j,j) &
                   +areavol*(cntot(j1,1)+cntot(j1,2)+cntot(j1,3) &
                            +cntot(j1,4)+cntot(j1,5)+cntot(j1,6))
        enddo
      endif
      if(.not.if2d) then
        do j1=jupstt,ng
          j=j1-jupstt+1
          ccgr(j,j)=ccgr(j,j)+rhz*(cntz0o(j1)+cntz1o(j1))
        enddo
      endif
      call invag(ng-jupstt+1,ccgr,rcgr)
      do j=1,ng-jupstt+1
        fcgr=0
        do j2=1,ng-jupstt+1
          fcgr=fcgr+rcgr(j2,j)*scgr(j2)
        enddo
        j3=j+jupstt-1
        hflx(j3)=fcgr*hflx(j3)
        do ii=1,6
          cntot(j3,ii)=fcgr*cntot(j3,ii)
          aflxt(j3,ii)=fcgr*aflxt(j3,ii)
          xmomt(j3,ii)=fcgr*xmomt(j3,ii)
          ymomt(j3,ii)=fcgr*ymomt(j3,ii)
        enddo
      enddo
 1000 enddo
      iter=iter-1
  999 continue
#ifdef HPPA
c
c fixes to prevent underflow on HPs
      do ii=1,6
        do j=1,ng
          if(abs(ymomt(j,ii)).lt.1.e-30) ymomt(j,ii)=0
        enddo 
      enddo
#endif

      return
      end

      subroutine invag(n,a,r)

!  inverse all group nxn matrix 

			use param
			implicit double precision (a-h,o-z)

!      dimension a(n,n),r(n,n),t(12,12)
      dimension a(n,n),r(n,n),t(n,n)
      do j=1,n
        do j2=1,n
          t(j2,j)=a(j2,j)
          r(j2,j)=0.
        enddo
        r(j,j)=1.
      enddo

      do j=1,n
! change by yunlin xu for DEBUG: overflow
!        rdet=1./t(j,j)
        IF (t(j,j) .EQ. 0.0) THEN
          rdet = 1.0e30
        ELSE
          rdet = 1./t(j,j)
        ENDIF
!end of change
        do j2=j+1,n
          t(j2,j)=rdet*t(j2,j)
        enddo
        do j2=1,j
          r(j2,j)=rdet*r(j2,j)
        enddo
        do j2=1,j-1
          tem=-t(j,j2)
          do j3=j+1,n
            t(j3,j2)=t(j3,j2)+tem*t(j3,j)
          enddo
          do j3=1,j
            r(j3,j2)=r(j3,j2)+tem*r(j3,j)
          enddo
        enddo
        do j2=j+1,n
          tem=-t(j,j2)
          do j3=j+1,n
            t(j3,j2)=t(j3,j2)+tem*t(j3,j)
          enddo
          do j3=1,j
            r(j3,j2)=r(j3,j2)+tem*r(j3,j)
          enddo
        enddo
      enddo

      return
      end
