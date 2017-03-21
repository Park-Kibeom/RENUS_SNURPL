      subroutine xsafp
!
! update nodal microscopic absorption cross sections for Xenon and Samarium
!
      use timer
      use MASTERXSL    ! 2014_08_12 . scb
      use param
      
      include 'global.h'
      include 'geom.h'
      include 'thfdbk.inc'
      include 'thgeom.inc'
      include 'times.h'
      include 'xesm.h'
      include 'xsec.h'
      include 'srchppm.h'   ! 2014_05_26 . scb
      include 'thcntl.inc'  ! 2014_05_26 . scb
      include 'thop.inc'    ! 2014_07_22 . scb
      
      data eps/0.001/
      dimension del(4)
      equivalence (del(1),delppm)
      equivalence (del(2),deltm)
      equivalence (del(3),deldm)
      equivalence (del(4),deltf)
      
      logical,save :: first=.true.
      
      INTEGER :: ID=1     ! 2016.12.22. JJH  

      call timeron()     
      
! 2014_08_12 . scb      
      if(first) then
        first=.false.
        if(flagmasxsl) then
          ice=DPPM
          if(fdbk) ice=DTF
          
          do k=1,nz
            ka=ktoka(k)
            do l=1,nxy
              la=ltola(l)
              iat=iassytyp(la)
              icomp=icompnumz(ka,iat)
              
              do m=1,ng                                               ! 2016.12.22. JJH  
                sigxea(m,icomp)=ORIGINXS(ICOMP,ABSO,IXE,m)*1.e-24
                sigsma(m,icomp)=ORIGINXS(ICOMP,ABSO,ISM,m)*1.e-24
                
                dsigxea(DPPM,m,icomp)=PPMDXS(ICOMP,ABSO,IXE,m,ID)*1.e-24
                dsigsma(DPPM,m,icomp)=PPMDXS(ICOMP,ABSO,ISM,m,ID)*1.e-24
              enddo
              
              if(fdbk) then
                do m=1,ng                
                  dsigxea(DTM,m,icomp)=0.D0
                  dsigsma(DTM,m,icomp)=0.D0
                  
                  dsigxea(DDM,m,icomp)=DMDXS(ICOMP,ABSO,IXE,m,5)*1.e-24
                  dsigsma(DDM,m,icomp)=DMDXS(ICOMP,ABSO,ISM,m,5)*1.e-24
                  
                  dsigxea(DTF,m,icomp)=TFDXS(ICOMP,ABSO,IXE,m,ID)*1.e-24
                  dsigsma(DTF,m,icomp)=TFDXS(ICOMP,ABSO,ISM,m,ID)*1.e-24
                enddo
              endif
            enddo
          enddo
        endif
      endif
! added end

      do k=1,nz
        ka=ktoka(k)
          do l=1,nxy
            la=ltola(l)
            iat=iassytyp(la)
            ia=latoia(la)
            ja=latoja(la)
            icomp=icompnumz(ka,iat)
            icr=abs(irodtyp(la)) ! STUCK ROD
                        
            do m=1,ng
              xsxeaf(m,l,k)= sigxea(m,icomp)
              xssmaf(m,l,k)= sigsma(m,icomp)
            enddo
            
            if(icr.ne.0) then
              do m=1,ng
                xsxeaf(m,l,k)= xsxeaf(m,l,k) + rodfrac(k,icr)*dcontxea(m,icr)
                xssmaf(m,l,k)= xssmaf(m,l,k) + rodfrac(k,icr)*dcontsma(m,icr)
              enddo
            endif            

! 2014_08_12 . scb
            if(flagmasxsl) then
              if(icm2ica(icomp).gt.0) then
                BASECOND(DPPM) = RPPM(icm2ica(ICOMP))
                BASECOND(DTF)  = RTF(icm2ica(ICOMP))
                BASECOND(DTM)  = RTM(icm2ica(ICOMP))
                BASECOND(DDM)  = RDM(icm2ica(ICOMP))
              else
                basecond=0.d0  
              endif   
            END IF
            
! added end            
            ice=DPPM
            delppm=ppm-basecond(DPPM)  ! 2014_05_26 . scb
            
            if(fdbk) then
              kth=ktokth(k)
              lchan=ltochan(l)
              
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
              deltm=tcool(kth,lchan)-basecond(DTM)
              deldm=dcool(kth,lchan)-basecond(DDM)
              if(flagmasxsl)  deldm=deldm*1.e-3   ! 2014_08_12 . scb
              deltf=tdopl(kth,lchan)-basecond(DTF)   
              
! 2014_08_12 . SCB       
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
              !do ifb=1,4
              do ifb=dppm,ice   ! 2014_05_26 . scb
                xsxeaf(m,l,k)=xsxeaf(m,l,k) + dsigxea(ifb,m,icomp)*del(ifb)
                xssmaf(m,l,k)=xssmaf(m,l,k) + dsigsma(ifb,m,icomp)*del(ifb)
              enddo
            enddo
          enddo
      enddo
      !endif
      
      call timeroff(txsec)
      
      return
      end
