SUBROUTINE UpdateXsec_Depletion(L,K,IFFDBK,IFTRAN)

! calculated nodal cross sections at a given t/h condition
      use MASTERXSL  ! 2014_05_22 . pkb
      use param
      use deplmod
      use timer
      use decusping1n,  only : xsset ! 2012_09_28 . scb
      use xsec,         only : xsrsp3, xstfsp3   ! 2014_10_07 . scb
      use readN2A,       only : assm, maxnb
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
!
      logical, intent(in) :: iffdbk
      logical, intent(in) :: iftran  ! 2012_09_28 . scb
      integer, intent(in) :: l
      integer, intent(in) :: k
!
      logical, save :: first=TRUE
!      
      real :: del(NUMOFDXS-1)
      equivalence(del(1), delppm)
      equivalence(del(2), deltm)
      equivalence(del(3), deldm)
      equivalence(del(4), deltf)

      real, external  :: calbchxs1st, calbchxs2nd      

! 2012_09_28 . scb
      type(xsset)      :: xsurod(ng), xsdel(ng), xsgen(ng)
      real             :: xssurod(ng,ng), xssdel(ng,ng), xssgen(ng,ng)
      logical          :: flagdecusp, dcspOK
! added end

      real :: delxsa ! 2013_10_01 . scb
! 2013_10_18 . scb for dbg      
      real :: delxea, delsma, delxeat, delsmat, delxsat  ! 2013_10_18 . scb for dbg
      
      real :: fnorm

      real :: branch(maxnb), nbranch
      integer :: ms, id=1          ! 2016.10.11.  jjh
      
      common / xedbg / fnorm
! added end

      integer,save :: nrdirdf=4  ! 2014_05_22 . scb
      LOGICAL :: FLAGUPDCRXS=.FALSE.    ! 2014_08_21 . SCB
      
      nrdirdf=4
      if(hex)  nrdirdf=6      
      
      if(flagmasxsl) then
        usesigtr=.true.  ! 2014_05_22 . pkb
        usesiga=.true.   ! 2014_05_22 . pkb
      endif
      
      sumdkp=0
! update xsec
        ka=ktoka(k)
          la=ltola(l)
          iat=iassytyp(la)
          if(rect) then  ! 2014_05_22 . scb
            ia=latoia(la)
            ja=latoja(la)
          endif          
          icomp=icompnumz(ka,iat)

          do m=1,ng         
            if(usesigtr) then
              xstrf(m,l,k)=NodeTR(m,l,k)   ! 2013_07_17 . scb
              xsdf(m,l,k)=1.d0/(3.d0*xstrf(m,l,k))
              xsdf2(m,l,k)=9.d0/7.d0*xsdf(m,l,k)
            else
              xsdf(m,l,k)=NodeDI(m,l,k)
              xstrf(m,l,k)=1.d0/(3.d0*xsdf(m,l,k))              
              xsdf2(m,l,k)=9.d0/7.d0*sigd(m,icomp)   ! 2013_07_15 . scb
            endif
            
            if(usesiga) then
              xsaf(m,l,k)=NodeAB(m,l,k)
            else
!              xstf(m,l,k)=sigt(m,icomp)
            endif
            
            xsnff(m,l,k)=NodeNF(m,l,k)
            xsff(m,l,k)=NodeFI(m,l,k)    ! 2014_07_31 . scb
            xsnuf(m,l,k)=xsnff(m,l,k)/(xsff(m,l,k)+1.e-20)   ! 2014_10_24 . scb            
            if(ifsp3 .and. usesigt_p3)  xstfsp3(m,l,k)=sigt_p3(m,icomp)   ! 2014_11_28 . scb
            xskpf(m,l,k)=NodeKF(m,l,k)
!            xssf(sigsms(m,icomp):sigsme(m,icomp),m,l,k) = sigsm(sigsms(m,icomp):sigsme(m,icomp),m,icomp)
            do ms=1,ng
              if(ms.eq.m) then
              else
                xssf(ms,m,l,k) = NodeSC(m,l,k)
              endif
            enddo
           
            if(flagadfmas) then
              do irdir=1,nrdirdf
                xsadf(irdir,m,l,k) = NodeADF(l,k,m)
                xscdf(irdir,m,l,k) = NodeADF(l,k,m)
              enddo
            endif                
          enddo
          !xssf(1,2,l,k)=macorigin(5,5,2)   ! 2014_08_05 . scb for dbg       
          
! 2014_05_22 . pkb          
          if(flagmasxsl) then
            if(icm2ica(icomp).gt.0) then
              BASECOND(DPPM) = RPPM(icm2ica(ICOMP))
              BASECOND(DTF)  = RTF(icm2ica(ICOMP))
              BASECOND(DTM)  = RTM(icm2ica(ICOMP))
              BASECOND(DDM)  = RDM(icm2ica(ICOMP))
            else
              !basecond=0.d0  ! 2014_07_31 . scb   
            endif   
! added end
          END IF
! ADDED END
          del(:) = 0.
          delppm=ppm-basecond(DPPM)
          if(iffdbk) then
            kth=ktokth(k)
            lchan=ltochan(l)              
            deltm=tcool(kth,lchan)-basecond(DTM)   ! 2014_08_01 . scb removed version dependency
            deldm=dcool(kth,lchan)-basecond(DDM)   ! 2014_05_22 . pkb
            if(flagmasxsl)  deldm=deldm*1.e-3   ! 2014_08_01 . scb
            deltf=tdopl(kth,lchan)-basecond(DTF)               
            deldm = 0._8        ! deldm means pressure depression, and deltm includes delta of H2O density

! 2014_08_07 . SCB       
            if(flagmasxsl) then     
              IF(ICM2ICA(ICOMP).GT.MASCOM .AND. ICM2ICA(ICOMP).LE.MASCOM+3) THEN
                IXS=IXSDERV(ICM2ICA(ICOMP)-MASCOM)
                IF(IXS.EQ.0) THEN
                  DELTM = DELPPM
                ELSEIF(IXS.NE.1) THEN
                  STOP 'IXS should be 0 or 1 !'
                ELSEIF(ixsecver.eq.4) THEN
                  DELTM=DELDM
                ENDIF
              ENDIF
            ENDIF
! ADDED END
          endif
            
          do m=1,ng
! 2014_05_22 . scb
            IF(N2AXS) call Updden_Depl(iffdbk,m,l,k,lchan,kth,icomp,del)  ! 2016_11_23 . pkb
            if(flagadfmas) then
              do irdir=1,nrdirdf                  
                xsadf(irdir,m,l,k) = xsadf(irdir,m,l,k) + NppmdADF0(l,k,m)*del(dppm)&
                                                        + NtfdADF0(l,k,m)*del(dtf)&
                                                        + NdmdADF0(l,k,m)*del(ddm)&
                                                        + NtmdADF0(l,k,m)*del(dtm)
                xscdf(irdir,m,l,k) = xscdf(irdir,m,l,k) + NppmdADF0(l,k,m)*del(dppm)&
                                                        + NtfdADF0(l,k,m)*del(dtf)&
                                                        + NdmdADF0(l,k,m)*del(ddm)&
                                                        + NtmdADF0(l,k,m)*del(dtm)
              enddo
            endif              
! added end            
            if(usesigtr) then
              xstrf(m,l,k)=xstrf(m,l,k)+Nodedppm_tr0(m,l,k)*del(dppm)&
                                       +Nodedtf_tr0(m,l,k)*del(dtf)&
                                       +Nodedtm_tr0(m,l,k)*del(dtm)&
                                       +Nodeddm_tr0(m,l,k)*del(ddm)
              xsdf(m,l,k)=1.d0/(3.d0*xstrf(m,l,k))
            else
!              xsdf(m,l,k)=xsdf(m,l,k)+dsigd_tr(icond,m,icomp)*del(icond)
!              xstrf(m,l,k)=1.d0/(3.d0*xsdf(m,l,k))        
            endif
            xsdf2(m,l,k)=9.d0/7.d0*xsdf(m,l,k)   ! 2013_07_15 . scb

            if(usesiga) then
              xsaf(m,l,k)=xsaf(m,l,k)+Nodedppm_ab0(m,l,k)*del(dppm)&
                                     +Nodedtf_ab0(m,l,k)*del(dtf)&
                                     +Nodedtm_ab0(m,l,k)*del(dtm)&
                                     +Nodeddm_ab0(m,l,k)*del(ddm)
            else
!              xstf(m,l,k)=xstf(m,l,k)+dsigt_a(icond,m,icomp)*del(icond)
            endif

            xsnff(m,l,k)=xsnff(m,l,k)+Nodedppm_nf0(m,l,k)*del(dppm)&
                                     +Nodedtf_nf0(m,l,k)*del(dtf)&
                                     +Nodedtm_nf0(m,l,k)*del(dtm)&
                                     +Nodeddm_nf0(m,l,k)*del(ddm)
            xsff(m,l,k)=xsff(m,l,k)+Nodedppm_fi0(m,l,k)*del(dppm)&
                                   +Nodedtf_fi0(m,l,k)*del(dtf)&
                                   +Nodedtm_fi0(m,l,k)*del(dtm)&
                                   +Nodeddm_fi0(m,l,k)*del(ddm)
            xsnuf(m,l,k)=xsnff(m,l,k)/(xsff(m,l,k)+1.e-20)
            xskpf(m,l,k)=xskpf(m,l,k)+Nodedppm_kf0(m,l,k)*del(dppm)&
                                     +Nodedtf_kf0(m,l,k)*del(dtf)&
                                     +Nodedtm_kf0(m,l,k)*del(dtm)&
                                     +Nodeddm_kf0(m,l,k)*del(ddm)
            do ms=1,ng
              if(ms.eq.m) then
              else
                xssf(ms,m,l,k) = xssf(ms,m,l,k)+Nodedppm_sc0(m,l,k)*del(dppm)&
                                               +Nodedtf_sc0(m,l,k)*del(dtf)&
                                               +Nodedtm_sc0(m,l,k)*del(dtm)&
                                               +Nodeddm_sc0(m,l,k)*del(ddm)
              endif
            enddo
            
            if(ifsp3 .and. usesigt_p3)  xstfsp3(m,l,k)=xstfsp3(m,l,k)+dsigt_p3(icond,m,icomp)*del(icond)   ! 2014_11_28 . scb
          ENDDO
! 2014_08_06 . SCB         
          DO M=1,NG  
! ADDED END
            
	          xbeta(m,l,k,XDIR)=xsdf(m,l,k)*rhmesh(XDIR,l,k)
            xbeta(m,l,k,YDIR)=xsdf(m,l,k)*rhmesh(YDIR,l,k)
	          xbeta(m,l,k,ZDIR)=xsdf(m,l,k)*rhmesh(ZDIR,l,k)
! 2013_07_15 . scb            
	          xbeta2(m,l,k,XDIR)=xsdf2(m,l,k)*rhmesh(XDIR,l,k)
            xbeta2(m,l,k,YDIR)=xsdf2(m,l,k)*rhmesh(YDIR,l,k)
	          xbeta2(m,l,k,ZDIR)=xsdf2(m,l,k)*rhmesh(ZDIR,l,k)
! added end     
          enddo !m      
! 15_11_16 . scb
      if(rect) then
        iat=iassytyp(la)
        icomp=icompnumz(ka,iat)
        nmx=nmeshx(latoia(la));   nmx1=nmx-1
        nmy=nmeshy(latoja(la));   nmy1=nmy-1
        do ly=1,nmy;  do lx=1,nmx1
          li1=nmx*(ly-1)+lx;     li2=nmx*(ly-1)+lx+1
          l1=latol(la)%fm(li1);  l2=latol(la)%fm(li2)
          do m=1,ng
            xsadf(2,m,l1,k)=1. ;  xsadf(1,m,l2,k)=1.
          enddo
        enddo; enddo
                
        do ly=1,nmy1;  do lx=1,nmx
          li1=nmx*(ly-1)+lx;     li2=nmx*ly+lx
          l1=latol(la)%fm(li1);  l2=latol(la)%fm(li2)
          do m=1,ng
            xsadf(4,m,l1,k)=1. ;  xsadf(3,m,l2,k)=1.
          enddo
        enddo; enddo
      endif      
! add end

1000  continue      
! 2013_09_30 . scb      
      if(flagxesm .and. .not.first) then
        delxea=0.d0
        delsma=0.d0
        delxeat=0.d0
        delsmat=0.d0      
        delxsat=0.d0  
        do m=1,ng
          delxsa=rnxe(l,k)*xsxeaf(m,l,k) + rnsm(l,k)*xssmaf(m,l,k)
          xsaf(m,l,k)=xsaf(m,l,k) + delxsa
          xstf(m,l,k)=xstf(m,l,k) + delxsa 
          !xsaf(m,l,k)=xsaf(m,l,k)*1.e3  ! 2013_10_02 . scb 
          !delxsat=delxsat+delxsa
        enddo      
      endif
! 2012_10_02
!      initxsec=true   ! 2012_09_28 . scb
      if(.not.first)  initxsec=true
! added end   
            
! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
      write (8,'(a,1p,6e16.4)'), 'SUM OF XSEC', sum(xsdf),sum(xsaf),sum(xstf),sum(xsnff),sum(xskpf)
      print '(a,1p,6e16.4)', 'SUM OF XSEC', sum(xsdf),sum(xsaf),sum(xstf),sum(xsnff),sum(xskpf) 
! added end

2000 continue
! Initialize the dxsec variables - 2016_10_17 . pkb
      If(N2AXS) then
        dsigd_tr1 = dsigd_tr10
        dsigt_a1 = dsigt_a10
        dsignf1 = dsignf10
        dsigf1 = dsigf10
        dsigkf1 = dsigkf10
        dsigsm1 = dsigsm10
        if(ifsp3 .and. usesigt_p3)  dsigt_p31 = dsigt_p310        

!        dsigd_tr(:,:,:) = dsigd_tr0(:,:,:)      
!        dsigt_a(:,:,:) = dsigt_a0(:,:,:)
!        dsignf(:,:,:) = dsignf0(:,:,:)
!        dsigf(:,:,:) = dsigf0(:,:,:)
!        dsigkf(:,:,:) = dsigkf0(:,:,:)
!        dsigsm(:,:,:,:) = dsigsm0(:,:,:,:)
!        if(ifsp3 .and. usesigt_p3)  dsigt_p3(:,:,:) = dsigt_p30(:,:,:)
      Endif
      first=FALSE      
      return
END SUBROUTINE