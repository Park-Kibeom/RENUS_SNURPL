!    subroutine updxsec(iffdbk)
!#define XE_SM_DBG  
    subroutine updxsec(iffdbk,iftran)   ! 2012_09_28 . scb
! calculated nodal cross sections at a given t/h condition
      use MASTERXSL  ! 2014_05_22 . pkb
      use param
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
      
      call timeron()     
      
      DEL=0.D0
      
! 2014_08_08 . scb      
      !if(first) then
      nrdirdf=4
      if(hex)  nrdirdf=6      
      !endif
! added end      
      
      if(flagmasxsl) then
        usesigtr=.true.  ! 2014_05_22 . pkb
        usesiga=.true.   ! 2014_05_22 . pkb
      endif      
        
! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
      if(iffdbk) then
!        if(ifmslb) then         ! MSLB  
!          call xsecfb11       ! MSLB        
!          goto 1000           ! MSLB 
        if(iflfr) then
          call updxsecf(iffdbk)
          goto 1000
        !elseif(.not.rect) then      ! 2014_05_22 . scb  commented
        !  call updxsec_h(iffdbk)		
        !  !return
        !  goto 2000
        endif                   ! MSLB  
      else
!        if(ifmslb) then
!          return              ! MSLB
        if(iflfr) then
          call initxsecf	
          goto 1000
!        elseif(.not.rect) then
!! 2012_11_06 . scb             ! 2014_05_22 . scb  commented
!!          call initxsecf	
!          call updxsec_h(iffdbk)	
!! added end          
!          !return
!          goto 2000
        endif
      endif 
! added end     
	   
      if(ixsecver.eq.2) then
!#ifdef MOXBENCH
        call updmoxbchxs(iffdbk,calbchxs2nd)   ! 2015_07_31 . scb recovered MOXBENCH
!#else
        !call terminate('ERROR : Compile RENUS again with a preprocessor "MOXBENCH"')
!#endif                  
        goto 1000
      endif
      
! 2012_09_28 . scb      
      flagdecusp=.false.
      dcspOK=.false.
! added end
      
      sumdkp=0
! update xsec

      If(N2AXS) then
        dsigd_tr10 = dsigd_tr1
        dsigt_a10 = dsigt_a1
        dsignf10 = dsignf1
        dsigf10 = dsigf1
        dsigkf10 = dsigkf1
        dsigsm10 = dsigsm1
        if(ifsp3 .and. usesigt_p3)  dsigt_p310 = dsigt_p31
 
!        dsigd_tr0(:,:,:) = dsigd_tr(:,:,:)      
!        dsigt_a0(:,:,:) = dsigt_a(:,:,:)
!        dsignf0(:,:,:) = dsignf(:,:,:)
!        dsigf0(:,:,:) = dsigf(:,:,:)
!        dsigkf0(:,:,:) = dsigkf(:,:,:)
!        dsigsm0(:,:,:,:) = dsigsm(:,:,:,:)
!        if(ifsp3 .and. usesigt_p3)  dsigt_p30(:,:,:) = dsigt_p3(:,:,:)
      Endif

      do k=1,nz
        ka=ktoka(k)
        if(first) then
          do la=1,nxya
            iat=iassytyp(la)
            icomp=icompnumz(ka,iat)
! 2014_05_22 . scb            
            !nfmsqrt=sqrt(latol(la)%nfm+0.0)
            !nfmm=latol(la)%nfm
            if(rect) then
!              nfmsqrt=sqrt(latol(la)%nfm+0.0)   ! 15_11_16 . scb
              nfmm=latol(la)%nfm
            else
              nfmm=1
              nfmsqrt=1
            endif
! added end            
              
            do li=1,nfmm
! 2014_05_22 . scb              
              !l=latol(la)%fm(li)     
              l=la
              if(rect) l=latol(la)%fm(li)     
! added end              
              xsmax(l,k)=sigsmax(icomp)   ! 2013_07_19 . scb              
              do m=1,ng         
                xschif(m,l,k)=sigchi(m,icomp)
                xsadf(:,m,l,k)=sigadf(:,m,icomp)
!  1      2
!   ------
!  |      |
!  |      |   : cdf index
!   ------
!  3      4         
                
! 2014_05_22 . scb
                !xscdf(:,m,l,k)=sigmdf(:,m,icomp)     
                !  
                !if(li.eq.1) xscdf(1,m,l,k)=sigcdf(1,m,icomp)
                !if(li.eq.nfmsqrt) xscdf(2,m,l,k)=sigcdf(1,m,icomp)
                !if(li.eq.(latol(la)%nfm-nfmsqrt+1)) xscdf(3,m,l,k)=sigcdf(1,m,icomp)
                !if(li.eq.latol(la)%nfm) xscdf(4,m,l,k)=sigcdf(1,m,icomp)
                if(rect) then
                  xscdf(:,m,l,k)=sigmdf(:,m,icomp)     

! 15_11_16 . scb                  
                  !if(li.eq.1) xscdf(1,m,l,k)=sigcdf(1,m,icomp)
                  !if(li.eq.nfmsqrt) xscdf(2,m,l,k)=sigcdf(1,m,icomp)
                  !if(li.eq.(latol(la)%nfm-nfmsqrt+1)) xscdf(3,m,l,k)=sigcdf(1,m,icomp)
                  !if(li.eq.latol(la)%nfm) xscdf(4,m,l,k)=sigcdf(1,m,icomp)
! delete end
                else
                  xscdf(:,m,l,k)=sigcdf(:,m,icomp)     
                endif                
! added end
                xssfs(m,l,k)=sigsms(m,icomp)
                xssfe(m,l,k)=sigsme(m,icomp) 
!                xscdf(:,m,l,k) = xsadf(:,m,l,k)

              enddo ! m
! 2013_09_27 . scb
              if(flagxesm) then  ! 2014_07_29 . scb
                gami(l,k)=gammafp(1,icomp)
                gamxe(l,k)=gammafp(2,icomp)
                gampm(l,k)=gammafp(3,icomp)
              endif              
! added end              
            enddo ! li
          enddo ! la
        endif ! first


        do l=1,nxy
          la=ltola(l)
          iat=iassytyp(la)
          if(rect) then  ! 2014_05_22 . scb
            ia=latoia(la)
            ja=latoja(la)
          endif          
          icomp=icompnumz(ka,iat)
          irodtyp1=abs(irodtyp(la)) ! STUCK ROD

          do m=1,ng         
            if(usesigtr) then
              !xstrfdcsp(m,l,k)=sigtr(m,icomp)  ! 2012_09_28 . scb
              xstrf(m,l,k)=sigtr(m,icomp)   ! 2013_07_17 . scb
!              xstrf(1,l,k)=0.225454700181442
!              xstrf(2,l,k)=0.732923196422403
              xsdf(m,l,k)=1.d0/(3.d0*xstrf(m,l,k))
              xsdf2(m,l,k)=9.d0/7.d0*xsdf(m,l,k)
            else
              xsdf(m,l,k)=sigd(m,icomp)
              xstrf(m,l,k)=1.d0/(3.d0*xsdf(m,l,k))              
              xsdf2(m,l,k)=9.d0/7.d0*sigd(m,icomp)   ! 2013_07_15 . scb
            endif
            
            if(usesiga) then
              xsaf(m,l,k)=siga(m,icomp)
!              xsaf(1,l,k)=8.032042332815962D-003
!              xsaf(2,l,k)=5.149747351622991D-002
            else
              xstf(m,l,k)=sigt(m,icomp)
            endif
            
            xsnff(m,l,k)=signf(m,icomp)
!            xsnff(1,l,k)=4.610836515403681D-003
!            xsnff(2,l,k)=7.586319486475608D-002
            xsff(m,l,k)=sigf(m,icomp)    ! 2014_07_31 . scb
!            xsff(1,l,k)=1.777886346588922D-003
!            xsff(2,l,k)=3.113358792618405D-002
            xsnuf(m,l,k)=xsnff(m,l,k)/(xsff(m,l,k)+1.e-20)   ! 2014_10_24 . scb            
            
            if(ifsp3 .and. usesigt_p3)  xstfsp3(m,l,k)=sigt_p3(m,icomp)   ! 2014_11_28 . scb
            !xsff(m,l,k)=xsnff(m,l,k)/2.3  ! 2013_10_02 . scb
            xskpf(m,l,k)=sigkf(m,icomp)
!            xskpf(1,l,k)=5.886254306165759D-014
!            xskpf(2,l,k)=1.009195402626632D-012
            xssf(sigsms(m,icomp):sigsme(m,icomp),m,l,k) = sigsm(sigsms(m,icomp):sigsme(m,icomp),m,icomp)
!            xssf(1,2,l,k)=1.846712791704523D-002
!            xssf(2,1,l,k)=1.100339853740579D-003
                
! 2014_08_08 . scb            
! 2014_05_22 . pkb            
            if(flagadfmas) then
              do irdir=1,nrdirdf
                xsadf(irdir,m,l,k) = originadf(icomp,m)
                xscdf(irdir,m,l,k) = originadf(icomp,m)
              enddo
            endif                
! added end        
! added end
          enddo
          !xssf(1,2,l,k)=macorigin(5,5,2)   ! 2014_08_05 . scb for dbg 

! control rod            
          if(irodtyp1 .ne. 0) then
            if(rodfrac(k,irodtyp1) .ne. 0) then
!!!!! scv added   
              flagdecusp=decusp.and.rodfrac(k,irodtyp1).lt.1.and.iffuelc(icomp)           
              flagdecusp=flagdecusp .and. rect   ! 2014_05_22 . scb
              if(flagdecusp .and. initxsec) then
                do m=1,ng
                  if(usesigtr) then
                    !xsurod(m)%xstr=xstrfdcsp(m,l,k)
                    xsurod(m)%xstr=xstrf(m,l,k)
                  else
                    xsurod(m)%xstr=xsdf(m,l,k)
                  endif
                  xsdel(m)%xstr=dsigd_tr(drod,m,icomp)
                  
                  if(usesiga) then
                    xsurod(m)%xsa=xsaf(m,l,k)
                  else
                    xsurod(m)%xsa=xstf(m,l,k)
                  endif
                  xsdel(m)%xsa=dsigt_a(drod,m,icomp)
                  
                  xsurod(m)%xsd =0.
                  xsurod(m)%xsf =0.
                  xsurod(m)%xst =0.
                  xsurod(m)%xsnf=xsnff(m,l,k)
                  xsurod(m)%xsf=xsff(m,l,k)   ! 2014_07_31 . scb
                  xsurod(m)%xskp=xskpf(m,l,k)
                  
                  xsdel(m)%xsd =0.
                  xsdel(m)%xsf =0.
                  xsdel(m)%xst =0.
                  xsdel(m)%xsnf=dsignf(drod,m,icomp)
                  xsdel(m)%xsf=dsigf(drod,m,icomp)   ! 2014_07_31 . scb
                  xsdel(m)%xskp=dsigkf(drod,m,icomp)
                  
                  xssurod=0.
                  xssdel =0.
                  do ms=1,ng
                    xssurod(ms,m)=xssf(ms,m,l,k)
                    xssdel(ms,m) =dsigsm(ms,drod,m,icomp)
                  enddo
                enddo
                
                call decusp1n(iftran,l,k,rodfrac(k,irodtyp1),usesigtr,usesiga, &
                             xsurod,xssurod,xsdel,xssdel,xsgen,xssgen,dcspok)
     
              endif
         
              if(dcspok) then
                dcspok=.false.                
                flagdecusp=.false.
                do m=1,ng
                  if(usesigtr) then
                    !xstrfdcsp(m,l,k)=xsgen(m)%xstr
                    !xstrf(m,l,k)=xstrfdcsp(m,l,k)   ! 2013_07_17 . scb
                    !xsdf(m,l,k)=1/(3*xstrfdcsp(m,l,k))
                    xstrf(m,l,k)=xsgen(m)%xstr
                    xsdf(m,l,k)=1.d0/(3.d0*xstrf(m,l,k))
                  else
                    xsdf(m,l,k)=xsgen(m)%xstr
                    xstrf(m,l,k)=1.d0/(3.d0*xsdf(m,l,k))        
                  endif
                  xsdf2(m,l,k)=9.d0/7.d0*xsdf(m,l,k)   ! 2013_07_15 . scb
                  
                  if(usesiga) then
                    xsaf(m,l,k)=xsgen(m)%xsa
                  else
                    xstf(m,l,k)=xsgen(m)%xsa
                  endif
                  
                  xsnff(m,l,k)=xsgen(m)%xsnf
                  xsff(m,l,k)=xsgen(m)%xsf    ! 2014_07_31 . scb
                  xsnuf(m,l,k)=xsnff(m,l,k)/(xsff(m,l,k)+1.e-20)   ! 2014_10_24 . scb                  
                  !xsff(m,l,k)=xsnff(m,l,k)/2.3  ! 2013_10_02 . scb                  
                  xskpf(m,l,k)=xsgen(m)%xskp
                  do ms=1,ng
                    xssf(ms,m,l,k)=xssgen(ms,m)
                  enddo
                enddo
              else
!!!!! added end 
                do m=1,ng
                  if(usesigtr) then
                    !xstrfdcsp(m,l,k)=xstrf(m,l,k)+rodfrac(k,irodtyp1)*dsigd_tr(drod,m,icomp)
                    !xstrf(m,l,k)=xstrfdcsp(m,l,k)   ! 2013_07_17 . scb
                    !xsdf(m,l,k)=1.d0/(3.d0*xstrfdcsp(m,l,k))        
                    xstrf(m,l,k)=xstrf(m,l,k)+rodfrac(k,irodtyp1)*dsigd_tr(drod,m,icomp)
                    xsdf(m,l,k)=1.d0/(3.d0*xstrf(m,l,k))         
                  else
                    xsdf(m,l,k)=xsdf(m,l,k)+rodfrac(k,irodtyp1)*dsigd_tr(drod,m,icomp)
                    xstrf(m,l,k)=1.d0/(3.d0*xsdf(m,l,k))        
                  endif
                  xsdf2(m,l,k)=9.d0/7.d0*xsdf(m,l,k)   ! 2013_07_15 . scb
                
                  if(usesiga) then
                    xsaf(m,l,k)=xsaf(m,l,k)+rodfrac(k,irodtyp1)*dsigt_a(drod,m,icomp)
                  else
                    xstf(m,l,k)=xstf(m,l,k)+rodfrac(k,irodtyp1)*dsigt_a(drod,m,icomp)
                  endif
                  
                  xsnff(m,l,k)=xsnff(m,l,k)+rodfrac(k,irodtyp1)*dsignf(drod,m,icomp)
                  xsff(m,l,k)=xsff(m,l,k)+rodfrac(k,irodtyp1)*dsigf(drod,m,icomp)   ! 2014_07_31 . SCB
                  xsnuf(m,l,k)=xsnff(m,l,k)/(xsff(m,l,k)+1.e-20)   ! 2014_10_24 . scb                  
                  !xsff(m,l,k)=xsnff(m,l,k)/2.3  ! 2013_10_02 . scb                  
                  xskpf(m,l,k)=xskpf(m,l,k)+rodfrac(k,irodtyp1)*dsigkf(drod,m,icomp)
                  do ms=sigsms(m,icomp),sigsme(m,icomp)
                    xssf(ms,m,l,k)=xssf(ms,m,l,k)+dsigsm(ms,drod,m,icomp)*rodfrac(k,irodtyp1)
                  enddo
                  
                  if(ifsp3 .and. usesigt_p3)  xstfsp3(m,l,k)=xstfsp3(m,l,k)+rodfrac(k,irodtyp1)*dsigt_p3(drod,m,icomp)   ! 2014_11_28 . scb
                enddo
              endif ! not decusp      
            endif ! rodfrac(k,irodtyp1) .ne. 0
          endif ! irodtyp1 .ne. 0      
          
! 2014_05_22 . pkb          
          if(flagmasxsl) then
          !if((ixsecver.eq.3 .or. ixsecver.eq.4) .and. iffuelc(icomp)) then   ! 2014_06_16 . scb
! 2014_07_23 . scb            
            !BASECOND(DPPM) = RPPM(ICOMP)
            !BASECOND(DTF)  = RTF(ICOMP)
            !BASECOND(DTM)  = RTM(ICOMP)
            !BASECOND(DDM)  = RDM(ICOMP)
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

          delppm=ppm-basecond(DPPM)
          ice=DPPM
          if(iffdbk) then
            kth=ktokth(k)
            lchan=ltochan(l)
! 2014_04_11 . scb              
            !if(iffuela(iat)) then
            !  deltm=tcool(kth,lchan)-basecond(DTM)
            !  deldm=dcool(kth,lchan)-basecond(DDM)
            !  deltf=tdopl(kth,lchan)-basecond(DTF)
            !else
            !  deltm=tin-basecond(DTM)
            !  deldm=din-basecond(DDM)
            !  deltf=tdopin-basecond(DTF)
            !endif
            !deltm=tcool(kth,lchan)-basecond(DTM)
            !deldm=dcool(kth,lchan)-basecond(DDM)
#ifndef DLL_TASS            
! 2014_08_01 . scb for dbg              
!              tcool(kth,lchan)=275.0
!              tdopl(kth,lchan)=tdopin
!              din=772.08
!              dcool(kth,lchan)=772.08 
!! 3D original              
!            tcool(kth,lchan)=275.01
!            tdopl(kth,lchan)=sqrt(275.15+CKELVIN)
!            din=772.071
!            dcool(kth,lchan)=772.071     
!              
!! 0D RTYA1
!            !tcool(kth,lchan)=275.00
!            !tdopl(kth,lchan)=sqrt(275.31+CKELVIN)
!            !din=772.078
!            !dcool(kth,lchan)=772.078
!! CONTROL ROD TEST
!            tcool(kth,lchan)=275.0
!            tdopl(kth,lchan)=sqrt(275.00+CKELVIN)
!            din=773.132
!            dcool(kth,lchan)=773.132          
!! 3D original              
!            tcool(kth,lchan)=294.26
!            tdopl(kth,lchan)=sqrt(360.95+CKELVIN)
!            din=735.103
!            dcool(kth,lchan)=735.103
!! added end              
            !tcool(kth,lchan)=275.00
            !tdopl(kth,lchan)=sqrt(275.00+CKELVIN)
            !din=772.073
            !dcool(kth,lchan)=772.073
! compare with TASS/ARTOS            
            !tcool(kth,lchan)=567.402-CKELVIN
            !tdopl(kth,lchan)=sqrt(634.154)
            !din=735.147
            !dcool(kth,lchan)=735.147
#endif            

              
            deltm=tcool(kth,lchan)-basecond(DTM)   ! 2014_08_01 . scb removed version dependency
            deldm=dcool(kth,lchan)-basecond(DDM)   ! 2014_05_22 . pkb
            if(flagmasxsl)  deldm=deldm*1.e-3   ! 2014_08_01 . scb
            deltf=tdopl(kth,lchan)-basecond(DTF)               
!#define TMOD_DBG
#ifdef TMOD_DBG
  test =.true.
  deltm =  0._8  ! 8.7  17.64 25.65 33.53
  deldm = 0._8
!  NDcorrect=.true.   ! only using deldm test
  deltf = 0._8  !3.133292 !5.903772 !9.265623 !13.90190


#endif
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
! test ! if there is no problem , below lines should be removed..              
!#ifndef DLL
!              ierr=0
!              if(.not.iffuela(iat)) then
!                if(lchan.ne.nchan+1) ierr=1
!                if(tin.ne.tcool(kth,lchan)) ierr=1
!                if(din.ne.dcool(kth,lchan)) ierr=1
!                if(tdopin.ne.tdopl(kth,lchan)) ierr=1
!              endif
!              if(ierr.eq.1) stop 'there is error in inlet T/H condition'
!#endif              
! added end      
            ice=DTF
          endif
            
          do m=1,ng
! density of boron, temperature of moderator
! density of moderator, temperature of fuel
            
            do icond=DPPM,ice ! NUMOFDXS-1
! 2014_05_22 . scb
              IF(N2AXS) call UpdDen(iffdbk,m,l,k,lchan,kth,icomp,icond,del)  ! 2016_11_23 . pkb                     
              if(flagadfmas) then
                do irdir=1,nrdirdf                  
                  xsadf(irdir,m,l,k) = xsadf(irdir,m,l,k) + DADFACTOR(icond,M,ICOMP)*del(icond)  ! 2014_05_22 . pkb
                  xscdf(irdir,m,l,k) = xscdf(irdir,m,l,k) + DADFACTOR(icond,M,ICOMP)*del(icond)  ! 2014_05_22 . pkb
                enddo
              endif              
! added end              
              if(usesigtr) then
                !xstrfdcsp(m,l,k)=xstrf(m,l,k)+dsigd_tr(icond,m,icomp)*del(icond)
                !xstrf(m,l,k)=xstrfdcsp(m,l,k)   ! 2013_07_17 . scb
                !xsdf(m,l,k)=1.d0/(3.d0*xstrfdcsp(m,l,k))
                xstrf(m,l,k)=xstrf(m,l,k)+dsigd_tr(icond,m,icomp)*del(icond)
                xsdf(m,l,k)=1.d0/(3.d0*xstrf(m,l,k))
              else
                xsdf(m,l,k)=xsdf(m,l,k)+dsigd_tr(icond,m,icomp)*del(icond)
                xstrf(m,l,k)=1.d0/(3.d0*xsdf(m,l,k))        
              endif
              xsdf2(m,l,k)=9.d0/7.d0*xsdf(m,l,k)   ! 2013_07_15 . scb

              if(usesiga) then
                xsaf(m,l,k)=xsaf(m,l,k)+dsigt_a(icond,m,icomp)*del(icond)
              else
                xstf(m,l,k)=xstf(m,l,k)+dsigt_a(icond,m,icomp)*del(icond)
              endif

              xsnff(m,l,k)=xsnff(m,l,k)+dsignf(icond,m,icomp)*del(icond)
              xsff(m,l,k)=xsff(m,l,k)+dsigf(icond,m,icomp)*del(icond)   ! 2014_07_31 . SCB
              xsnuf(m,l,k)=xsnff(m,l,k)/(xsff(m,l,k)+1.e-20)   ! 2014_10_24 . scb              
              !xsff(m,l,k)=xsnff(m,l,k)/2.3  ! 2013_10_02 . scb              
              xskpf(m,l,k)=xskpf(m,l,k)+dsigkf(icond,m,icomp)*del(icond)
              sumdkp=sumdkp+del(icond)
              do ms=sigsms(m,icomp),sigsme(m,icomp)
               xssf(ms,m,l,k)=xssf(ms,m,l,k)+dsigsm(ms,icond,m,icomp)*del(icond)
              enddo
              
              if(ifsp3 .and. usesigt_p3)  xstfsp3(m,l,k)=xstfsp3(m,l,k)+dsigt_p3(icond,m,icomp)*del(icond)   ! 2014_11_28 . scb
              
            enddo ! icond
            
            IF(flagmasxsl .AND. IFFDBK) THEN
                IF(IFMASFUEL(ICOMP)) THEN
!                    CALL DUPDXS1(DEL,DCOOL(KTH,LCHAN),BASECOND(DDM),ICOMP,M,L,K)
                ELSEIF(IFAX(ICOMP)) THEN
                    CALL DUPDXS_AXIAL(DCOOL(KTH,LCHAN),BASECOND(DDM),ICOMP,M,L,K)
                ENDIF
            END IF
          ENDDO
                                  
! 2014_08_22 . scb modified
! 2014_08_21 . PKB
          IF(FLAGCCR) THEN
            L2=IAASS1DINV(L)
            if(ng.ne.2)  stop 'Corner control rod modeling is valid only for 2G'
            
            DO IC=1,6
              FLAGUPDCRXS = FLAGUPDCRXS .OR. (ASSEMCONTROL(L2,IC).NE.0 .AND. RODFRAC(K,ASSEMCONTROL(L2,IC)).GT.0.)
            ENDDO
            
            IF(.NOT.FLAGUPDCRXS) GOTO 822
             
            CALL CRDXS(IFFDBK,DEL,L,K,L2)
               
            DO M=1,NG              
              xstrf(m,l,k)=xstrf(m,l,k)+ORIGINCR(M,L,K,TRAN)&
                                        +PPMCR(M,L,K,TRAN)*DEL(DPPM)&
                                        +DMDCR(M,L,K,TRAN)*DEL(DDM)
              !xstrfdcsp(m,l,k)=xstrf(m,l,k)
              xsdf(m,l,k)=1.d0/(3.d0*xstrf(m,l,k))
              xsdf2(m,l,k)=9.d0/7.d0*xsdf(m,l,k)
              xsaf(m,l,k)=xsaf(m,l,k)+ORIGINCR(M,L,K,ABSO)&
                                      +PPMCR(M,L,K,ABSO)*DEL(DPPM)&
                                      +DMDCR(M,L,K,ABSO)*DEL(DDM)
              xsnff(m,l,k)=xsnff(m,l,k)+ORIGINCR(M,L,K,NUFS)&
                                        +PPMCR(M,L,K,NUFS)*DEL(DPPM)&
                                        +DMDCR(M,L,K,NUFS)*DEL(DDM)
              xsff(m,l,k)=xsff(m,l,k)+ORIGINCR(M,L,K,FISS)&
                                      +PPMCR(M,L,K,FISS)*DEL(DPPM)&
                                      +DMDCR(M,L,K,FISS)*DEL(DDM)
              xsnuf(m,l,k)=xsnff(m,l,k)/(xsff(m,l,k)+1.e-20)   ! 2014_10_24 . scb              
              xskpf(m,l,k)=xskpf(m,l,k)+ORIGINCR(M,L,K,KAPA)&
                                        +PPMCR(M,L,K,KAPA)*DEL(DPPM)&
                                        +DMDCR(M,L,K,KAPA)*DEL(DDM)
            ENDDO
            
            xssf(1,2,l,k)=xssf(1,2,l,k)+ORIGINCR(2,L,K,SCA)    &
                                       +PPMCR(2,L,K,SCA)*DEL(DPPM)&
                                       +DMDCR(2,L,K,SCA)*DEL(DDM)
            
822         CONTINUE            
          END IF          
! ADDED END           
! added end
          
! 2014_08_05 . PKB
            !IF(IXSECVER.EQ.3 .OR. IXSECVER.EQ.4) THEN
            !    IF(IFFDBK .AND. IFMASFUEL(ICOMP)) THEN
            !        CALL DUPDXS1(DEL,DCOOL(KTH,LCHAN),BASECOND(DDM),ICOMP,M,L,K)
            !    ELSEIF(IFFDBK .AND. (IFMASFUEL(ICOMP).EQ. .FALSE.)) THEN
            !        !PRINT*,'TEST'
            !    ENDIF
            !END IF
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
          
        enddo !l
      enddo !k
      
! 15_11_16 . scb
      if(rect) then
        do k=1,nz
          do la=1,nxya
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
          enddo
        enddo
      endif      
! add end
      
      do k=1,nz
        do l=1,nxy
          do m=1,ng
            summ=0.d0
            do md=1,ng
              if(m .ge. xssfs(md,l,k) .and. m.le.xssfe(md,l,k)) then
                summ=summ+xssf(m,md,l,k)
              endif
            enddo
            if(usesiga) then
              xstf(m,l,k)=xsaf(m,l,k)+summ
              !write(723,*) m,l,k,xsaf(m,l,k)  ! 2014_07_23 . scb 
            else
              xsaf(m,l,k)=xstf(m,l,k)-summ
            endif
          enddo
        enddo
      enddo
      
1000  continue
           
! 2013_10_18 . scb for dbg
#ifdef XE_SM_DBG
      if(first) then
        open(unit=1018,file='XEDBG_1D',status='unknown')
        write(1018,'(a5,4(a15))') 'TRAN','DEL_XEA','DEL_SMA','DEL_XSA', 'FNORM'
      endif
#endif      
      
! 2013_09_30 . scb      
      if(flagxesm .and. .not.first) then
        delxea=0.d0
        delsma=0.d0
        delxeat=0.d0
        delsmat=0.d0      
        delxsat=0.d0  
        do k=1,nz
          do l=1,nxy
            do m=1,ng
#ifdef XE_SM_DBG
              delxea=rnxe(l,k)*xsxeaf(m,l,k)
              delsma=rnsm(l,k)*xssmaf(m,l,k)
              delxsa=delxea+delsma
              
              delxeat=delxeat+delxea
              delsmat=delsmat+delsma
              
              delxsat=delxsat+delxsa
              xsaf(m,l,k)=xsaf(m,l,k) + delxsa
              xstf(m,l,k)=xstf(m,l,k) + delxsa 
#else              
              delxsa=rnxe(l,k)*xsxeaf(m,l,k) + rnsm(l,k)*xssmaf(m,l,k)
              xsaf(m,l,k)=xsaf(m,l,k) + delxsa
              xstf(m,l,k)=xstf(m,l,k) + delxsa 
#endif      
              !xsaf(m,l,k)=xsaf(m,l,k)*1.e3  ! 2013_10_02 . scb 
              !delxsat=delxsat+delxsa
            enddo
          enddo
        enddo      
      endif
      
! 2014_10_07 . scb, 2014_11_28 . scb added usesigt_p3
      if(ifsp3) then
        !if(.not. usesigt_p3) then    ! 2015_01_09 . scb commented
        do k=1,nz
          do l=1,nxy
            do m=1,ng
              !xstfsp3(m,l,k)=(xstrf(m,l,k)-xsaf(m,l,k))/xsrsp3(m,1) + xsaf(m,l,k)  
              if(.not. usesigt_p3 .or. xstfsp3(m,l,k) .lt. 1.e-15) xstfsp3(m,l,k)=(xstrf(m,l,k)-xsaf(m,l,k))/xsrsp3(m,1) + xsaf(m,l,k)  ! 2015_01_09 . scb
            enddo
          enddo
        enddo         
        !endif

        do k=1,nz
          do l=1,nxy
            do m=1,ng
              xsdf2(m,l,k)=3./(7.*xstfsp3(m,l,k))
            enddo
          enddo
        enddo    
      endif
! added end      
      
        
! 2013_10_18 . scb for dbg      
#ifdef XE_SM_DBG
      if(.not.first) then
        if(iftran) then
          write(1018,'(a5,4es15.5)') 'T',delxeat,delsmat,delxsat,fnorm
        else
          write(1018,'(a5,4es15.5)') 'F',delxeat,delsmat,delxsat,fnorm
        endif
      endif
#endif      
! added end      

! 2012_10_02
!      initxsec=true   ! 2012_09_28 . scb
      if(.not.first)  initxsec=true
! added end   
            
! added in ARTOS ver. 0.2 . 2012_07_11 by SCB.
      write (8,'(a,1p,6e16.4)'), 'SUM OF XSEC', sum(xsdf),sum(xsaf),sum(xstf),sum(xsnff),sum(xskpf)
      print '(a,1p,6e16.4)', 'SUM OF XSEC', sum(xsdf),sum(xsaf),sum(xstf),sum(xsnff),sum(xskpf) 
! added end

!! 2014_09_05 . scb for dbg
!      xsnffsum=0.d0
!      do k=1,nz
!        do l=1,nxy
!          do m=1,ng
!            xsnffsum=xsnffsum+xsnff(m,l,k)
!            write(905,*) m, l, k, xsnff(m,l,k)
!          enddo
!        enddo
!      enddo
!      write(905,*) xsnffsum
!! added end      

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
! added end                
      first=FALSE
! added end      
     
      call timeroff(txsec)
      
      return
    end subroutine