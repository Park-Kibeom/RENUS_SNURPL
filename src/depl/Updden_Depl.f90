Subroutine Updden_Depl(iffdbk,m,l,k,lchan,kth,icomp,del)

      use MASTERXSL  ! 2014_05_22 . pkb
      use param
      use timer
      use decusping1n,  only : xsset ! 2012_09_28 . scb
      use xsec,         only : xsrsp3, xstfsp3   ! 2014_10_07 . scb
      use readN2A,       only : assm    ! 2016.10.12   jjh
      use deplmod
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

      Logical :: iffdbk
      Integer :: m,l,k,lchan,kth,icomp
      INTEGER :: ID                          ! 2016.12.22. JJH
      Real :: nH2O,nB10,fact,del(iver)

      fact = 0.998921
      nB10 = ppm*B10_ppm
      nH2O = Atomn(l,k,IH2O)
      If(iffdbk) then
        if (test==.false.) then
          nB10 = dcool(kth,lchan)/basecond(DDM)*nB10
          nH2O = dcool(kth,lchan)/basecond(DDM)*nH2O
        else if (test==.true.) then 
!          nB10 = (basecond(ddm)+del(ddm)*1000)/basecond(DDM)*nB10
!          nH2O = (basecond(ddm)+del(ddm)*1000)/basecond(DDM)*nH2O
!          if (NDcorrect==.true.) then   ! correction of the difference due to ND calculation way between nTRACER and ARTOS
!            nB10 = nB10 * fact
!            nH2O = nH2O * fact
!          end if
        end if
      Endif

      If(iffdbk .eq. .false.) then
        xstrf(m,l,k)=xstrf(m,l,k)+ NodeMicxs(l,k,TRAN,IH2O,m)*nH2O &
                                 + NodeMicxs(l,k,TRAN,IB10,m)*nB10
        xsaf(m,l,k) =xsaf(m,l,k) + NodeMicxs(l,k,ABSO,IH2O,m)*nH2O &
                                 + NodeMicxs(l,k,ABSO,IB10,m)*nB10
        xsnff(m,l,k)=xsnff(m,l,k)+ NodeMicxs(l,k,NUFS,IH2O,m)*nH2O &
                                 + NodeMicxs(l,k,NUFS,IB10,m)*nB10
        xsff(m,l,k) =xsff(m,l,k) + NodeMicxs(l,k,FISS,IH2O,m)*nH2O &
                                 + NodeMicxs(l,k,FISS,IH2O,m)*nB10
        xskpf(m,l,k)=xskpf(m,l,k)+ NodeMicxs(l,k,KAPA,IH2O,m)*nH2O &
                                 + NodeMicxs(l,k,KAPA,IB10,m)*nB10
        do ms=1,ng
          if(ms.eq.m) then
          else
            xssf(ms,m,l,k) = xssf(ms,m,l,k)+NodeMicxs(l,k,SCA,IH2O,m)*nH2O &
                                           +NodeMicxs(l,k,SCA,IB10,m)*nB10
          endif
        enddo

        Do id=1,NBORON(icomp)
          Nodedppm_tr(m,l,k,id)  = Nodedppm_tr(m,l,k,id) +& 
                                   NodeDppm(l,k,TRAN,IH2O,m,ID)*nH2O + NodeDppm(l,k,TRAN,IB10,m,ID)*nB10
          Nodedppm_ab(m,l,k,id)  = Nodedppm_ab(m,l,k,id) +&
                                   NodeDppm(l,k,ABSO,IH2O,m,ID)*nH2O + NodeDppm(l,k,ABSO,IB10,m,ID)*nB10
          Nodedppm_nf(m,l,k,id)  = Nodedppm_nf(m,l,k,id) +&                                       
                                   NodeDppm(l,k,NUFS,IH2O,m,ID)*nH2O + NodeDppm(l,k,NUFS,IB10,m,ID)*nB10
          Nodedppm_kf(m,l,k,id)  = Nodedppm_kf(m,l,k,id) +&                                          
                                   NodeDppm(l,k,KAPA,IH2O,m,ID)*nH2O + NodeDppm(l,k,KAPA,IB10,m,ID)*nB10
          Nodedppm_fi(m,l,k,id)  = Nodedppm_fi(m,l,k,id) +&                                          
                                   NodeDppm(l,k,FISS,IH2O,m,ID)*nH2O + NodeDppm(l,k,FISS,IB10,m,ID)*nB10
          DO MS=1,NG
            IF(MS.EQ.m) THEN
            ELSE
              Nodedppm_sc(m,l,k,id) = Nodedppm_sc(m,l,k,id) +&
                                      NodeDppm(l,k,SCA,IH2O,m,ID)*nH2O + NodeDppm(l,k,SCA,IB10,m,ID)*nB10
            END IF
          END DO
        End do

        call Decide_delta(dppm,boronstep(icomp,:),nboron(icomp),del(dppm), Nodedppm_tr(m,l,k,:), Nodedppm_tr0(m,l,k))
        call Decide_delta(dppm,boronstep(icomp,:),nboron(icomp),del(dppm), Nodedppm_ab(m,l,k,:), Nodedppm_ab0(m,l,k))
        call Decide_delta(dppm,boronstep(icomp,:),nboron(icomp),del(dppm), Nodedppm_nf(m,l,k,:), Nodedppm_nf0(m,l,k))
        call Decide_delta(dppm,boronstep(icomp,:),nboron(icomp),del(dppm), Nodedppm_kf(m,l,k,:), Nodedppm_kf0(m,l,k))
        call Decide_delta(dppm,boronstep(icomp,:),nboron(icomp),del(dppm), Nodedppm_fi(m,l,k,:), Nodedppm_fi0(m,l,k))
        call Decide_delta(dppm,boronstep(icomp,:),nboron(icomp),del(dppm), NppmdADF(l,k,m,:), NppmdADF0(l,k,m))
        Do ms=1,ng
          If(ms .eq. m) then
          Else        
            call Decide_delta(dppm,boronstep(icomp,:),nboron(icomp),del(dppm), Nodedppm_sc(m,l,k,:), Nodedppm_sc0(m,l,k))
          Endif
        Enddo

      Elseif(iffdbk .eq. .true.) then
        xstrf(m,l,k)=xstrf(m,l,k)+ NodeMicxs(l,k,TRAN,IH2O,m)*nH2O &
                                 + NodeMicxs(l,k,TRAN,IB10,m)*nB10
        xsaf(m,l,k) =xsaf(m,l,k) + NodeMicxs(l,k,ABSO,IH2O,m)*nH2O &
                                 + NodeMicxs(l,k,ABSO,IB10,m)*nB10
        xsnff(m,l,k)=xsnff(m,l,k)+ NodeMicxs(l,k,NUFS,IH2O,m)*nH2O &
                                 + NodeMicxs(l,k,NUFS,IB10,m)*nB10
        xsff(m,l,k) =xsff(m,l,k) + NodeMicxs(l,k,FISS,IH2O,m)*nH2O &
                                 + NodeMicxs(l,k,FISS,IH2O,m)*nB10
        xskpf(m,l,k)=xskpf(m,l,k)+ NodeMicxs(l,k,KAPA,IH2O,m)*nH2O &
                                 + NodeMicxs(l,k,KAPA,IB10,m)*nB10
        do ms=1,ng
          if(ms.eq.m) then
          else
            xssf(ms,m,l,k) = xssf(ms,m,l,k)+NodeMicxs(l,k,SCA,IH2O,m)*nH2O &
                                           +NodeMicxs(l,k,SCA,IB10,m)*nB10
          endif
        enddo

        Do id=1,NBORON(icomp)
          Nodedppm_tr(m,l,k,id)  = Nodedppm_tr(m,l,k,id) +& 
                                   NodeDppm(l,k,TRAN,IH2O,m,ID)*nH2O + NodeDppm(l,k,TRAN,IB10,m,ID)*nB10
          Nodedppm_ab(m,l,k,id)  = Nodedppm_ab(m,l,k,id) +&
                                   NodeDppm(l,k,ABSO,IH2O,m,ID)*nH2O + NodeDppm(l,k,ABSO,IB10,m,ID)*nB10
          Nodedppm_nf(m,l,k,id)  = Nodedppm_nf(m,l,k,id) +&                                       
                                   NodeDppm(l,k,NUFS,IH2O,m,ID)*nH2O + NodeDppm(l,k,NUFS,IB10,m,ID)*nB10
          Nodedppm_kf(m,l,k,id)  = Nodedppm_kf(m,l,k,id) +&                                          
                                   NodeDppm(l,k,KAPA,IH2O,m,ID)*nH2O + NodeDppm(l,k,KAPA,IB10,m,ID)*nB10
          Nodedppm_fi(m,l,k,id)  = Nodedppm_fi(m,l,k,id) +&                                          
                                   NodeDppm(l,k,FISS,IH2O,m,ID)*nH2O + NodeDppm(l,k,FISS,IB10,m,ID)*nB10
          DO MS=1,NG
            IF(MS.EQ.m) THEN
            ELSE
              Nodedppm_sc(m,l,k,id) = Nodedppm_sc(m,l,k,id) +&
                                      NodeDppm(l,k,SCA,IH2O,m,ID)*nH2O + NodeDppm(l,k,SCA,IB10,m,ID)*nB10
            END IF
          END DO
        End do

        call Decide_delta(dppm,boronstep(icomp,:),nboron(icomp),del(dppm), Nodedppm_tr(m,l,k,:), Nodedppm_tr0(m,l,k))
        call Decide_delta(dppm,boronstep(icomp,:),nboron(icomp),del(dppm), Nodedppm_ab(m,l,k,:), Nodedppm_ab0(m,l,k))
        call Decide_delta(dppm,boronstep(icomp,:),nboron(icomp),del(dppm), Nodedppm_nf(m,l,k,:), Nodedppm_nf0(m,l,k))
        call Decide_delta(dppm,boronstep(icomp,:),nboron(icomp),del(dppm), Nodedppm_kf(m,l,k,:), Nodedppm_kf0(m,l,k))
        call Decide_delta(dppm,boronstep(icomp,:),nboron(icomp),del(dppm), Nodedppm_fi(m,l,k,:), Nodedppm_fi0(m,l,k))
        call Decide_delta(dppm,boronstep(icomp,:),nboron(icomp),del(dppm), NppmdADF(l,k,m,:), NppmdADF0(l,k,m))
        Do ms=1,ng
          If(ms .eq. m) then
          Else        
            call Decide_delta(dppm,boronstep(icomp,:),nboron(icomp),del(dppm), Nodedppm_sc(m,l,k,:), Nodedppm_sc0(m,l,k))
          Endif
        Enddo

        Do id=1,NTFUEL(icomp)
          Nodedtf_tr(m,l,k,id)  = Nodedtf_tr(m,l,k,id) +& 
                                  NodeDtf(l,k,TRAN,IH2O,m,ID)*nH2O + NodeDtf(l,k,TRAN,IB10,m,ID)*nB10
          Nodedtf_ab(m,l,k,id)  = Nodedtf_ab(m,l,k,id) +&
                                  NodeDtf(l,k,ABSO,IH2O,m,ID)*nH2O + NodeDtf(l,k,ABSO,IB10,m,ID)*nB10
          Nodedtf_nf(m,l,k,id)  = Nodedtf_nf(m,l,k,id) +&                                       
                                  NodeDtf(l,k,NUFS,IH2O,m,ID)*nH2O + NodeDtf(l,k,NUFS,IB10,m,ID)*nB10
          Nodedtf_kf(m,l,k,id)  = Nodedtf_kf(m,l,k,id) +&                                          
                                  NodeDtf(l,k,KAPA,IH2O,m,ID)*nH2O + NodeDtf(l,k,KAPA,IB10,m,ID)*nB10
          Nodedtf_fi(m,l,k,id)  = Nodedtf_fi(m,l,k,id) +&                                          
                                  NodeDtf(l,k,FISS,IH2O,m,ID)*nH2O + NodeDtf(l,k,FISS,IB10,m,ID)*nB10
          DO MS=1,NG
            IF(MS.EQ.m) THEN
            ELSE
              Nodedtf_sc(m,l,k,id) = Nodedtf_sc(m,l,k,id) +&
                                     NodeDtf(l,k,SCA,IH2O,m,ID)*nH2O + NodeDtf(l,k,SCA,IB10,m,ID)*nB10
            END IF
          END DO
        End do
        
        call Decide_delta(dtf,tfuelstep(icomp,:),ntfuel(icomp),del(dtf), Nodedtf_tr(m,l,k,:), Nodedtf_tr0(m,l,k))
        call Decide_delta(dtf,tfuelstep(icomp,:),ntfuel(icomp),del(dtf), Nodedtf_ab(m,l,k,:), Nodedtf_ab0(m,l,k))
        call Decide_delta(dtf,tfuelstep(icomp,:),ntfuel(icomp),del(dtf), Nodedtf_nf(m,l,k,:), Nodedtf_nf0(m,l,k))
        call Decide_delta(dtf,tfuelstep(icomp,:),ntfuel(icomp),del(dtf), Nodedtf_kf(m,l,k,:), Nodedtf_kf0(m,l,k))
        call Decide_delta(dtf,tfuelstep(icomp,:),ntfuel(icomp),del(dtf), Nodedtf_fi(m,l,k,:), Nodedtf_fi0(m,l,k))
        call Decide_delta(dtf,tfuelstep(icomp,:),ntfuel(icomp),del(dtf), NtfdADF(l,k,m,:), NtfdADF0(l,k,m))
        Do ms=1,ng
          If(ms .eq. m) then
          Else        
            call Decide_delta(dtf,tfuelstep(icomp,:),ntfuel(icomp),del(dtf), Nodedtf_sc(m,l,k,:), Nodedtf_sc0(m,l,k))
          Endif
        Enddo
        
        Do id=1,NDMOD(icomp)
          Nodeddm_tr(m,l,k,id)  = Nodeddm_tr(m,l,k,id) +& 
                                  NodeDdm(l,k,TRAN,IH2O,m,ID)*nH2O + NodeDdm(l,k,TRAN,IB10,m,ID)*nB10
          Nodeddm_ab(m,l,k,id)  = Nodeddm_ab(m,l,k,id) +&
                                  NodeDdm(l,k,ABSO,IH2O,m,ID)*nH2O + NodeDdm(l,k,ABSO,IB10,m,ID)*nB10
          Nodeddm_nf(m,l,k,id)  = Nodeddm_nf(m,l,k,id) +&                                       
                                  NodeDdm(l,k,NUFS,IH2O,m,ID)*nH2O + NodeDdm(l,k,NUFS,IB10,m,ID)*nB10
          Nodeddm_kf(m,l,k,id)  = Nodeddm_kf(m,l,k,id) +&                                          
                                  NodeDdm(l,k,KAPA,IH2O,m,ID)*nH2O + NodeDdm(l,k,KAPA,IB10,m,ID)*nB10
          Nodeddm_fi(m,l,k,id)  = Nodeddm_fi(m,l,k,id) +&                                          
                                  NodeDdm(l,k,FISS,IH2O,m,ID)*nH2O + NodeDdm(l,k,FISS,IB10,m,ID)*nB10
          DO MS=1,NG
            IF(MS.EQ.m) THEN
            ELSE
              Nodeddm_sc(m,l,k,id) = Nodeddm_sc(m,l,k,id) +&
                                     NodeDdm(l,k,SCA,IH2O,m,ID)*nH2O + NodeDdm(l,k,SCA,IB10,m,ID)*nB10
            END IF
          END DO
        End do
        
        call Decide_delta(ddm,modstep(icomp,:),ndmod(icomp),del(ddm), Nodeddm_tr(m,l,k,:), Nodeddm_tr0(m,l,k))
        call Decide_delta(ddm,modstep(icomp,:),ndmod(icomp),del(ddm), Nodeddm_ab(m,l,k,:), Nodeddm_ab0(m,l,k))
        call Decide_delta(ddm,modstep(icomp,:),ndmod(icomp),del(ddm), Nodeddm_nf(m,l,k,:), Nodeddm_nf0(m,l,k))
        call Decide_delta(ddm,modstep(icomp,:),ndmod(icomp),del(ddm), Nodeddm_kf(m,l,k,:), Nodeddm_kf0(m,l,k))
        call Decide_delta(ddm,modstep(icomp,:),ndmod(icomp),del(ddm), Nodeddm_fi(m,l,k,:), Nodeddm_fi0(m,l,k))
        call Decide_delta(ddm,modstep(icomp,:),ndmod(icomp),del(ddm), NdmdADF(l,k,m,:), NdmdADF0(l,k,m))
        Do ms=1,ng
          If(ms .eq. m) then
          Else        
            call Decide_delta(ddm,modstep(icomp,:),ndmod(icomp),del(ddm), Nodeddm_sc(m,l,k,:), Nodeddm_sc0(m,l,k))
          Endif
        Enddo

        Do id=1,NTMOD(icomp)
          Nodedtm_tr(m,l,k,id)  = Nodedtm_tr(m,l,k,id) +& 
                                  NodeDtm(l,k,TRAN,IH2O,m,ID)*nH2O + NodeDtm(l,k,TRAN,IB10,m,ID)*nB10
          Nodedtm_ab(m,l,k,id)  = Nodedtm_ab(m,l,k,id) +&
                                  NodeDtm(l,k,ABSO,IH2O,m,ID)*nH2O + NodeDtm(l,k,ABSO,IB10,m,ID)*nB10
          Nodedtm_nf(m,l,k,id)  = Nodedtm_nf(m,l,k,id) +&                                       
                                  NodeDtm(l,k,NUFS,IH2O,m,ID)*nH2O + NodeDtm(l,k,NUFS,IB10,m,ID)*nB10
          Nodedtm_kf(m,l,k,id)  = Nodedtm_kf(m,l,k,id) +&                                          
                                  NodeDtm(l,k,KAPA,IH2O,m,ID)*nH2O + NodeDtm(l,k,KAPA,IB10,m,ID)*nB10
          Nodedtm_fi(m,l,k,id)  = Nodedtm_fi(m,l,k,id) +&                                          
                                  NodeDtm(l,k,FISS,IH2O,m,ID)*nH2O + NodeDtm(l,k,FISS,IB10,m,ID)*nB10
          DO MS=1,NG
            IF(MS.EQ.m) THEN
            ELSE
              Nodedtm_sc(m,l,k,id) = Nodedtm_sc(m,l,k,id) +&
                                     NodeDtm(l,k,SCA,IH2O,m,ID)*nH2O + NodeDtm(l,k,SCA,IB10,m,ID)*nB10
            END IF
          END DO
        End do
        
        call Decide_delta(dtm,tmodstep(icomp,:),ntmod(icomp),del(dtm), Nodedtm_tr(m,l,k,:), Nodedtm_tr0(m,l,k))
        call Decide_delta(dtm,tmodstep(icomp,:),ntmod(icomp),del(dtm), Nodedtm_ab(m,l,k,:), Nodedtm_ab0(m,l,k))
        call Decide_delta(dtm,tmodstep(icomp,:),ntmod(icomp),del(dtm), Nodedtm_nf(m,l,k,:), Nodedtm_nf0(m,l,k))
        call Decide_delta(dtm,tmodstep(icomp,:),ntmod(icomp),del(dtm), Nodedtm_kf(m,l,k,:), Nodedtm_kf0(m,l,k))
        call Decide_delta(dtm,tmodstep(icomp,:),ntmod(icomp),del(dtm), Nodedtm_fi(m,l,k,:), Nodedtm_fi0(m,l,k))
        call Decide_delta(dtm,tmodstep(icomp,:),ntmod(icomp),del(dtm), NtmdADF(l,k,m,:), NtmdADF0(l,k,m))
        Do ms=1,ng
          If(ms .eq. m) then
          Else        
            call Decide_delta(dtm,tmodstep(icomp,:),ntmod(icomp),del(dtm), Nodedtm_sc(m,l,k,:), Nodedtm_sc0(m,l,k))
          Endif
        Enddo

      End if

End Subroutine