Subroutine UpdDen(iffdbk,m,l,k,lchan,kth,icomp,icond,del)

      use MASTERXSL  ! 2014_05_22 . pkb
      use param
      use timer
      use decusping1n,  only : xsset ! 2012_09_28 . scb
      use xsec,         only : xsrsp3, xstfsp3   ! 2014_10_07 . scb
      use readN2A,       only : assm    ! 2016.10.12   jjh
!
      include 'global.h'
      include 'times.h'
      include 'xsec.h'
      include 'geom.h'
      include 'srchppm.h'
      include 'thfdbk.inc'
      include 'thfuel.inc'
      include 'thcntl.inc'
      include 'thgeom.inc'
      include 'thop.inc'
! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
      include 'files.h'
!      include 'mslb.inc' ! MSLB
! added end     
      include 'xesm.h'  ! 2013_09_27 . scb
      include 'trancntl.inc'   ! 2014_05_22 . scb
      INCLUDE 'GEOMH.H'     ! 2014_08_05 . PKB
      include 'defhex.h'    ! 2014_08_23 . scb

 ! Update Base and Gradient XS using moderator density
      
      Logical :: iffdbk
      Integer :: m,l,k,lchan,kth,icomp,icond
      INTEGER :: ID                          ! 2016.12.22. JJH
      Real :: nH2O,nB10,fact,keff4, keff5, del(iver)

      !fact = 2.38725e-8
      fact = 0.998921
      nB10 = ppm*B10_ppm
      nH2O = NUCND_T(ICOMP,IH2O)
      If(iffdbk) then
        if (test==.false.) then
          nB10 = dcool(kth,lchan)/basecond(DDM)*nB10
          nH2O = dcool(kth,lchan)/basecond(DDM)*nH2O
        else if (test==.true.) then 
          nB10 = (basecond(ddm)+del(ddm)*1000)/basecond(DDM)*nB10
          nH2O = (basecond(ddm)+del(ddm)*1000)/basecond(DDM)*nH2O
          if (NDcorrect==.true.) then   ! correction of the difference due to ND calculation way between nTRACER and ARTOS
            nB10 = nB10 * fact
            nH2O = nH2O * fact
          end if
        end if
      Endif
       keff4 = (XSnff(1,l,k) * (xsaf(2,l,k) + xssf(2,1,l,k)) + XSnff(2,l,k) * xssf(1,2,l,k) ) &          ! 2016.9.12.jjh
                / (xsaf(1,l,k) * (xsaf(2,l,k) + xssf(2,1,l,k)) + xsaf(2,l,k) * xssf(1,2,l,k))
        
      If(icond .eq. dppm) then
!------------------------------BASE---------------------------------------      
        xstrf(m,l,k)=xstrf(m,l,k)+ OriginXS(icomp,tran,ih2o,m)*nH2O &
                                 + OriginXS(icomp,tran,ib10,m)*nB10
        xsaf(m,l,k) =xsaf(m,l,k) + OriginXS(icomp,abso,ih2o,m)*nH2O &
                                 + OriginXS(icomp,abso,ib10,m)*nB10
        xsnff(m,l,k)=xsnff(m,l,k)+ OriginXS(icomp,nufs,ih2o,m)*nH2O &
                                 + OriginXS(icomp,nufs,ib10,m)*nB10
        xsff(m,l,k) =xsff(m,l,k) + OriginXS(icomp,fiss,ih2o,m)*nH2O &
                                 + OriginXS(icomp,fiss,ib10,m)*nB10
        xskpf(m,l,k)=xskpf(m,l,k)+ OriginXS(icomp,kapa,ih2o,m)*nH2O &
                                 + OriginXS(icomp,kapa,ib10,m)*nB10
        Do ms=1,ng
          If(ms .eq. m) then
          Else
            xssf(ms,m,l,k)=xssf(ms,m,l,k)+ OriginXS(icomp,sca,ih2o,m)*nH2O &
                                         + OriginXS(icomp,sca,ib10,m)*nB10
          Endif
        Enddo
        
       keff5 = (XSnff(1,l,k) * (xsaf(2,l,k) + xssf(2,1,l,k)) + XSnff(2,l,k) * xssf(1,2,l,k) ) &          ! 2016.9.12.jjh
                / (xsaf(1,l,k) * (xsaf(2,l,k) + xssf(2,1,l,k)) + xsaf(2,l,k) * xssf(1,2,l,k))
        
!----------------------------PPM GRADIENT---------------------------------                                ! 2016.12.22. JJH
        do id=1, nboron(icomp)
          dsigd_tr1(icond,ID,m,icomp)=dsigd_tr10(icond,ID,m,icomp)+ ppmdxs(icomp,tran,ih2o,m,ID)*nH2O &
                                                          + ppmdxs(icomp,tran,ib10,m,ID)*nB10
          dsigt_a1(icond,ID,m,icomp) =dsigt_a10(icond,ID,m,icomp) + ppmdxs(icomp,abso,ih2o,m,ID)*nH2O &
                                                          + ppmdxs(icomp,abso,ib10,m,ID)*nB10
          dsignf1(icond,ID,m,icomp)  =dsignf10(icond,ID,m,icomp)  + ppmdxs(icomp,nufs,ih2o,m,ID)*nH2O &
                                                          + ppmdxs(icomp,nufs,ib10,m,ID)*nB10
          dsigf1(icond,ID,m,icomp)   =dsigf10(icond,ID,m,icomp)   + ppmdxs(icomp,fiss,ih2o,m,ID)*nH2O &
                                                          + ppmdxs(icomp,fiss,ib10,m,ID)*nB10
          dsigkf1(icond,ID,m,icomp)  =dsigkf10(icond,ID,m,icomp)  + ppmdxs(icomp,kapa,ih2o,m,ID)*nH2O &
                                                          + ppmdxs(icomp,kapa,ib10,m,ID)*nB10
          Do ms=1,ng
            If(ms .eq. m) then
            Else
              dsigsm1(ms,icond,ID,m,icomp)=dsigsm10(ms,icond,ID,m,icomp)+ ppmdxs(icomp,sca,ih2o,m,ID)*nH2O &
                                                                + ppmdxs(icomp,sca,ib10,m,ID)*nB10
            Endif
          Enddo
        end do
      
        call Decide_delta(icond,boronstep(icomp,:),nboron(icomp),del(icond), dsigd_tr1(icond,:,m,icomp), dsigd_tr(icond,m,icomp))   
        call Decide_delta(icond,boronstep(icomp,:),nboron(icomp),del(icond), dsigt_a1(icond,:,m,icomp), dsigt_a(icond,m,icomp))           
        call Decide_delta(icond,boronstep(icomp,:),nboron(icomp),del(icond), dsignf1(icond,:,m,icomp), dsignf(icond,m,icomp))     
        call Decide_delta(icond,boronstep(icomp,:),nboron(icomp),del(icond), dsigf1(icond,:,m,icomp), dsigf(icond,m,icomp))     
        call Decide_delta(icond,boronstep(icomp,:),nboron(icomp),del(icond), dsigkf1(icond,:,m,icomp), dsigkf(icond,m,icomp))   
        call Decide_delta(icond,boronstep(icomp,:),nboron(icomp),del(icond), dadfactor0(icond,:,m,icomp), dadfactor(icond,m,icomp))  
        Do ms=1,ng
          If(ms .eq. m) then
          Else        
            call Decide_delta(icond,boronstep(icomp,:),nboron(icomp),del(icond), dsigsm1(ms,icond,:,m,icomp), dsigsm(ms,icond,m,icomp))     
          Endif
        Enddo            
        
!----------------------------TMOD GRADIENT------------------------------
      Elseif(icond .eq. dtm) then
        do id=1, ntmod(icomp)
          dsigd_tr1(icond,ID,m,icomp)=dsigd_tr10(icond,ID,m,icomp)+ tmdxs(icomp,tran,ih2o,m,ID)*nH2O &
                                                          + tmdxs(icomp,tran,ib10,m,ID)*nB10
          dsigt_a1(icond,ID,m,icomp) =dsigt_a10(icond,ID,m,icomp) + tmdxs(icomp,abso,ih2o,m,ID)*nH2O &
                                                          + tmdxs(icomp,abso,ib10,m,ID)*nB10
          dsignf1(icond,ID,m,icomp)  =dsignf10(icond,ID,m,icomp)  + tmdxs(icomp,nufs,ih2o,m,ID)*nH2O &
                                                          + tmdxs(icomp,nufs,ib10,m,ID)*nB10
          dsigf1(icond,ID,m,icomp)   =dsigf10(icond,ID,m,icomp)   + tmdxs(icomp,fiss,ih2o,m,ID)*nH2O &
                                                          + tmdxs(icomp,fiss,ib10,m,ID)*nB10
          dsigkf1(icond,ID,m,icomp)  =dsigkf10(icond,ID,m,icomp)  + tmdxs(icomp,kapa,ih2o,m,ID)*nH2O &
                                                          + tmdxs(icomp,kapa,ib10,m,ID)*nB10
          Do ms=1,ng
            If(ms .eq. m) then
            Else
              dsigsm1(ms,icond,ID,m,icomp)=dsigsm10(ms,icond,ID,m,icomp)+ tmdxs(icomp,sca,ih2o,m,ID)*nH2O &
                                                                + tmdxs(icomp,sca,ib10,m,ID)*nB10
            Endif
          Enddo
        end do
        
        call Decide_delta(icond,tmodstep(icomp,:),ntmod(icomp),del(icond), dsigd_tr1(icond,:,m,icomp), dsigd_tr(icond,m,icomp))   
        call Decide_delta(icond,tmodstep(icomp,:),ntmod(icomp),del(icond), dsigt_a1(icond,:,m,icomp), dsigt_a(icond,m,icomp))           
        call Decide_delta(icond,tmodstep(icomp,:),ntmod(icomp),del(icond), dsignf1(icond,:,m,icomp), dsignf(icond,m,icomp))     
        call Decide_delta(icond,tmodstep(icomp,:),ntmod(icomp),del(icond), dsigf1(icond,:,m,icomp), dsigf(icond,m,icomp))     
        call Decide_delta(icond,tmodstep(icomp,:),ntmod(icomp),del(icond), dsigkf1(icond,:,m,icomp), dsigkf(icond,m,icomp))
        call Decide_delta(icond,tmodstep(icomp,:),ntmod(icomp),del(icond), dadfactor0(icond,:,m,icomp), dadfactor(icond,m,icomp))  
        Do ms=1,ng
          If(ms .eq. m) then
          Else        
            call Decide_delta(icond,tmodstep(icomp,:),ntmod(icomp),del(icond), dsigsm1(ms,icond,:,m,icomp), dsigsm(ms,icond,m,icomp))     
          Endif
        Enddo     
        
        
!----------------------------DMOD GRADIENT------------------------------
      Elseif(icond .eq. ddm) then        
        do id=1, ndmod(icomp)
          dsigd_tr1(icond,ID,m,icomp)=dsigd_tr10(icond,ID,m,icomp)+ dmdxs(icomp,tran,ih2o,m,ID)*nH2O &
                                                          + dmdxs(icomp,tran,ib10,m,ID)*nB10
          dsigt_a1(icond,ID,m,icomp) =dsigt_a10(icond,ID,m,icomp) + dmdxs(icomp,abso,ih2o,m,ID)*nH2O &
                                                          + dmdxs(icomp,abso,ib10,m,ID)*nB10
          dsignf1(icond,ID,m,icomp)  =dsignf10(icond,ID,m,icomp)  + dmdxs(icomp,nufs,ih2o,m,ID)*nH2O &
                                                          + dmdxs(icomp,nufs,ib10,m,ID)*nB10
          dsigf1(icond,ID,m,icomp)   =dsigf10(icond,ID,m,icomp)   + dmdxs(icomp,fiss,ih2o,m,ID)*nH2O &
                                                          + dmdxs(icomp,fiss,ib10,m,ID)*nB10
          dsigkf1(icond,ID,m,icomp)  =dsigkf10(icond,ID,m,icomp)  + dmdxs(icomp,kapa,ih2o,m,ID)*nH2O &
                                                          + dmdxs(icomp,kapa,ib10,m,ID)*nB10
          Do ms=1,ng
            If(ms .eq. m) then
            Else
              dsigsm1(ms,icond,ID,m,icomp)=dsigsm10(ms,icond,ID,m,icomp)+ dmdxs(icomp,sca,ih2o,m,ID)*nH2O &
                                                                + dmdxs(icomp,sca,ib10,m,ID)*nB10
              Endif
          Enddo
        end do
        
        call Decide_delta(icond,modstep(icomp,:),ndmod(icomp),del(icond), dsigd_tr1(icond,:,m,icomp), dsigd_tr(icond,m,icomp))   
        call Decide_delta(icond,modstep(icomp,:),ndmod(icomp),del(icond), dsigt_a1(icond,:,m,icomp), dsigt_a(icond,m,icomp))           
        call Decide_delta(icond,modstep(icomp,:),ndmod(icomp),del(icond), dsignf1(icond,:,m,icomp), dsignf(icond,m,icomp))     
        call Decide_delta(icond,modstep(icomp,:),ndmod(icomp),del(icond), dsigf1(icond,:,m,icomp), dsigf(icond,m,icomp))     
        call Decide_delta(icond,modstep(icomp,:),ndmod(icomp),del(icond), dsigkf1(icond,:,m,icomp), dsigkf(icond,m,icomp)) 
        call Decide_delta(icond,modstep(icomp,:),ndmod(icomp),del(icond), dadfactor0(icond,:,m,icomp), dadfactor(icond,m,icomp))  
        Do ms=1,ng
          If(ms .eq. m) then
          Else        
            call Decide_delta(icond,modstep(icomp,:),ndmod(icomp),del(icond), dsigsm1(ms,icond,:,m,icomp), dsigsm(ms,icond,m,icomp))     
          Endif
        Enddo      
         
!----------------------------TFUEL GRADIENT------------------------------
      Elseif(icond .eq. dtf) then 
        do id=1, ntfuel(icomp)
          dsigd_tr1(icond,ID,m,icomp)=dsigd_tr10(icond,ID,m,icomp)+ tfdxs(icomp,tran,ih2o,m,ID)*nH2O &
                                                          + tfdxs(icomp,tran,ib10,m,ID)*nB10
          dsigt_a1(icond,ID,m,icomp) =dsigt_a10(icond,ID,m,icomp) + tfdxs(icomp,abso,ih2o,m,ID)*nH2O &
                                                          + tfdxs(icomp,abso,ib10,m,ID)*nB10
          dsignf1(icond,ID,m,icomp)  =dsignf10(icond,ID,m,icomp)  + tfdxs(icomp,nufs,ih2o,m,ID)*nH2O &
                                                          + tfdxs(icomp,nufs,ib10,m,ID)*nB10
          dsigf1(icond,ID,m,icomp)   =dsigf10(icond,ID,m,icomp)   + tfdxs(icomp,fiss,ih2o,m,ID)*nH2O &
                                                          + tfdxs(icomp,fiss,ib10,m,ID)*nB10
          dsigkf1(icond,ID,m,icomp)  =dsigkf10(icond,ID,m,icomp)  + tfdxs(icomp,kapa,ih2o,m,ID)*nH2O &
                                                          + tfdxs(icomp,kapa,ib10,m,ID)*nB10 
          Do ms=1,ng
            If(ms .eq. m) then
            Else
              dsigsm1(ms,icond,ID,m,icomp)=dsigsm10(ms,icond,ID,m,icomp)+ tfdxs(icomp,sca,ih2o,m,ID)*nH2O &
                                                                + tfdxs(icomp,sca,ib10,m,ID)*nB10
            Endif
          Enddo
        end do
        
        call Decide_delta(icond,tfuelstep(icomp,:),ntfuel(icomp),del(icond), dsigd_tr1(icond,:,m,icomp), dsigd_tr(icond,m,icomp))   
        call Decide_delta(icond,tfuelstep(icomp,:),ntfuel(icomp),del(icond), dsigt_a1(icond,:,m,icomp), dsigt_a(icond,m,icomp))           
        call Decide_delta(icond,tfuelstep(icomp,:),ntfuel(icomp),del(icond), dsignf1(icond,:,m,icomp), dsignf(icond,m,icomp))     
        call Decide_delta(icond,tfuelstep(icomp,:),ntfuel(icomp),del(icond), dsigf1(icond,:,m,icomp), dsigf(icond,m,icomp))     
        call Decide_delta(icond,tfuelstep(icomp,:),ntfuel(icomp),del(icond), dsigkf1(icond,:,m,icomp), dsigkf(icond,m,icomp)) 
        call Decide_delta(icond,tfuelstep(icomp,:),ntfuel(icomp),del(icond), dadfactor0(icond,:,m,icomp), dadfactor(icond,m,icomp)) 
        Do ms=1,ng
          If(ms .eq. m) then
          Else        
            call Decide_delta(icond,tfuelstep(icomp,:),ntfuel(icomp),del(icond), dsigsm1(ms,icond,:,m,icomp), dsigsm(ms,icond,m,icomp))     
          Endif
        Enddo     
        
        
      Endif

End Subroutine