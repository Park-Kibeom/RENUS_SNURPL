
  
SUBROUTINE N2A_SCANXSL(indev, icomp)              ! 2016.9.9. jjh

USE MASTERXSL
!USE TEMP_BURNUP
USE PARAM
use readN2A

INCLUDE 'GLOBAL.H'
INCLUDE 'BENCH.H'
INCLUDE 'XSEC.H'
INCLUDE 'SRCHPPM.H'
INCLUDE 'CARDS.H'

CHARACTER NAME*20
INTEGER :: CNUM, NAXIAL,NRADIAL
integer :: indev, icomp
integer :: io=11, i, g, a=1, m
logical :: first_alloc1=.true., first_alloc2=.true.


CNUM    = ncomp
NAXIAL = 2 
NRADIAL = 3

if (first_alloc1==.true.) then
  CALL N2A_ALLOC1(CNUM,NAXIAL,NRADIAL)
  first_alloc1=.false.
end if
REWIND(INDEV)
MASCOM = CNUM

call toupper(assmlist(icomp))
NADF(icomp)=nburn(icomp)*2
RPPM(icomp)=basecond(dppm)
RTM(icomp)=basecond(dtm)
RDM(icomp)=basecond(ddm)
RTF(icomp)=basecond(dtf)
RPRESS(icomp)=pressure
cname(icomp)=assmlist(icomp)
FLAGMASXSL=.TRUE.
!ixsecver=4            ! 2016.9.23. jjh

if (first_alloc2==.true.) then
  CALL N2A_ALLOC2(CNUM)
  CALL NUTYP               ! declare nuclide type and amass(i)
  first_alloc2=.false.
end if

!boronstep(icomp,1:nboron(icomp)) = boronstep1(1:nboron(icomp))
!tfuelstep(icomp,1:ntfuel(icomp)) = tfuelstep1(1:ntfuel(icomp))
!tmodstep(icomp,1:ntmod(icomp)) = tmodstep1(1:ntmod(icomp))
!modstep(icomp,1:ndmod(icomp)) = dmodstep1(1:ndmod(icomp))





REWIND(INDEV)
do g=1, ng
  do i=1, nburn(icomp)
    baseadf(icomp, i, g) = assm(icomp)%adf(i,g)
  end do
  do i=1, nder(icomp)
    do m=1, nboron(icomp)
      ppmadf(icomp, i, g,m) = assm(icomp)%dADF(i,g,m, brn)
    end do
  end do
  do i=1, nder(icomp)
    do m=1, ntfuel(icomp)
      tfadf(icomp, i, g,m) = assm(icomp)%dADF(i,g,m, tf)
    end do
  end do
  do i=1, nder(icomp)
    do m=1, ntmod(icomp)
      tmadf(icomp, i, g,m) = assm(icomp)%dADF(i,g,m, tm)
    end do
  end do
  do i=1, nder(icomp)
    do m=1, ntmod(icomp)
      dmadf(icomp, i, g,m) = assm(icomp)%dADF(i,g, m, dm)
    end do
  end do
  
!  do i=1, NUCNUM
!    KAPPA_T(icomp,i,g)=basexs(icomp,kap,1,i,g)
!    KAPPA_T(icomp,i,g)=assm(icomp)%nuclei(i)%xsec(kap)%basexs(1,g)
!  end do 
end do

do i=1, NUCNUM
  nucnd_t(icomp, i) = assm(icomp)%nuclei(i)%nd(a)    ! a=1
end do

do m=1, nburn(icomp)
  do i=1, NUCNUM
    nucnd_depl(icomp, i,m) = assm(icomp)%nuclei(i)%nd(m)    ! a=1
   ! write(*,*) 'a', nucnd_depl(icomp, i,m)
  end do
  !write(*,*) '----'
end do


!do g=1, ng  
!  do k=1, ntype2
!    do i=1, NUCNUM
!!      do j=1, nburn(icomp)
!!         basexs(icomp,k,j,i,g) = assm(icomp)%nuclei(i)%xsec(k)%basexs(j,g)
!!      end do
!      do j=1, nder(icomp)
!         do m=1, nboron(icomp)
!!           ppmxs(icomp,k,j,i,g,m)  = assm(icomp)%nuclei(i)%xsec(k)%dXS(j, g, m, brn)
!         end do
!         do m=1, ntfuel(icomp)
!!           tfxs(icomp,k,j,i,g,m)   = assm(icomp)%nuclei(i)%xsec(k)%dXS(j,g, m, tf)
!         end do
!         do m=1, ntmod(icomp)
!!          tmxs(icomp,k,j,i,g,m)   = assm(icomp)%nuclei(i)%xsec(k)%dXS(j, g, m, tm)
!         end do
!         do m=1, ntmod(icomp)
!!          dmxs(icomp,k,j,i,g,m)   = assm(icomp)%nuclei(i)%xsec(k)%dXS(j, g, m, dm)
!         end do
!      end do
!    end do
!  end do    
!end do

END SUBROUTINE  
  
  
  
!===============================================================================  
  
SUBROUTINE N2A_ALLOC1(CNUM,NAXIAL,NRADIAL)      ! 2016.9.9. jjh

USE MASTERXSL
USE TEMP_BURNUP
USE PARAM
use readN2A

INCLUDE 'GLOBAL.H'
INCLUDE 'BENCH.H'
INCLUDE 'XSEC.H'
INCLUDE 'SRCHPPM.H'
INCLUDE 'CARDS.H'

INTEGER,INTENT(IN)::CNUM,NAXIAL,NRADIAL

ALLOCATE(CNAME(CNUM))
ALLOCATE(NADF(CNUM))

ALLOCATE(RPPM(NCOMP))
ALLOCATE(RTF(NCOMP))
ALLOCATE(RTM(NCOMP))
ALLOCATE(RDM(NCOMP))
! added end
ALLOCATE(RPRESS(CNUM))
ALLOCATE(RFND(NAXIAL,NRC))
ALLOCATE(RFPPM(NRADIAL))
ALLOCATE(RFTEMP(NRADIAL))
ALLOCATE(RFMOD(NRADIAL))   
ALLOCATE(RFPRES(NRADIAL))
ALLOCATE(AXIAL(NAXIAL,NRC,NTYPE2,NG))
ALLOCATE(RADIAL(NRADIAL,NTYPE2,NG))
ALLOCATE(DTMRADIAL(NRADIAL,NTYPE2,NG))
ALLOCATE(MACAXIAL(NAXIAL,NRC,NTYPE2,NG))

cname=' '
nadf=0
RFND(:,:)         = 0.
AXIAL(:,:,:,:)    = 0.
MACAXIAL(:,:,:,:) = 0.
RADIAL(:,:,:)     = 0.
DTMRADIAL(:,:,:)  = 0.


! 2014_07_31 . scb
rppm=0.d0
rtf=0.d0
rtm=0.d0
rdm=0.d0
rpress=0.d0
! added end

! 2014_08_07 . SCB
ALLOCATE(IXSDERV(NRADIAL))
IXSDERV=0
! ADDED END

END SUBROUTINE


!==============================================================================


SUBROUTINE N2A_ALLOC2(CNUM)          ! 2016.9.9. jjh

USE MASTERXSL
USE TEMP_BURNUP
USE PARAM
use readN2A

INCLUDE 'GLOBAL.H'
INCLUDE 'BENCH.H'
INCLUDE 'XSEC.H'
INCLUDE 'SRCHPPM.H'
INCLUDE 'CARDS.H'

INTEGER,INTENT(IN)::CNUM

ALLOCATE(SBFLAG(CNUM))
ALLOCATE(SBURNSTEP(CNUM,MAXBURN,NUCNUM))
ALLOCATE(SDERSTEP(CNUM,MAXDER,NUCNUM))
ALLOCATE(NUCND_T(CNUM,NUCNUM))
ALLOCATE(NUCND_depl(CNUM,NUCNUM, maxburn))
ALLOCATE(NUCND_depl2(CNUM,NUCNUM, maxburn))

! 2014_07_28 . pkb
ALLOCATE(BASECRDXS(CNUM,NTYPE2,MAXBURN,NG))
ALLOCATE(PPMCRDXS(CNUM,NTYPE2,MAXDER,NG))
ALLOCATE(DMCRDXS(CNUM,NTYPE2,MAXDER,NG))
ALLOCATE(PPMDXS_H(CNUM,NTYPE2,NG,3))
ALLOCATE(DMDXS_H(CNUM,NTYPE2,NG,3,MAXDMOD))

ALLOCATE(CRDKAPPA(CNUM,NG,3))
ALLOCATE(BASEADF_CR(CNUM,MAXBURN,NG,3))
ALLOCATE(PPMADF_CR(CNUM,MAXDER,NG,3))
ALLOCATE(TMADF_CR(CNUM,MAXDER,NG,3))
ALLOCATE(DMADF_CR(CNUM,MAXDER,NG,3,MAXDMOD))
ALLOCATE(BASECRDXS_H(CNUM,NTYPE1,MAXDER,NG,3))
ALLOCATE(PPMCRDXS_H(CNUM,NTYPE1,MAXDER,NG,3))
ALLOCATE(DMCRDXS_H(CNUM,NTYPE1,MAXDER,NG,3,MAXDMOD))
! added end

ALLOCATE(NU2ID(NUCNUM))
ALLOCATE(AMASS(NUCNUM))
ALLOCATE(N2N_T(MAXBURN))
ALLOCATE(SBNUC(CNUM,NUCNUM))

ALLOCATE(ORIGINXS(NCOMP,NTYPE2,NUCNUM,NG))
ALLOCATE(PPMDXS(NCOMP,NTYPE2,NUCNUM,NG,MAXBORON))
ALLOCATE(TFDXS(NCOMP,NTYPE2,NUCNUM,NG,MAXTFUEL))
ALLOCATE(TMDXS(NCOMP,NTYPE2,NUCNUM,NG,MAXTMOD))
ALLOCATE(DMDXS(NCOMP,NTYPE2,NUCNUM,NG,MAXDMOD))

ALLOCATE(MACORIGIN(NCOMP,NTYPE2,NG))
ALLOCATE(MACDPPM(NCOMP,NTYPE2,NG,MAXBORON))
ALLOCATE(MACDTF(NCOMP,NTYPE2,NG,MAXTFUEL))
ALLOCATE(MACDTM(NCOMP,NTYPE2,NG,MAXTMOD))
ALLOCATE(MACDDM(NCOMP,NTYPE2,NG,MAXDMOD))

ALLOCATE(TEMPD(CNUM,NTYPE1,MAXBURN,NUCNUM,NG))

ALLOCATE(BASEADF(CNUM,MAXBURN,NG))
ALLOCATE(PPMADF(CNUM,MAXDER,NG,maxboron))
ALLOCATE(TFADF(CNUM,MAXDER,NG,maxtfuel))
ALLOCATE(TMADF(CNUM,MAXDER,NG,maxtmod))
ALLOCATE(DMADF(CNUM,MAXDER,NG,MAXDMOD))

ALLOCATE(ORIGINADF(NCOMP,NG))
ALLOCATE(PPMDADF(NCOMP,NG,maxboron))
ALLOCATE(TFDADF(NCOMP,NG,maxtfuel))
ALLOCATE(TMDADF(NCOMP,NG,maxtmod))
ALLOCATE(DMDADF(NCOMP,NG,MAXDMOD))

ALLOCATE(dsigt_a1(NUMOFDXS,max(maxboron,maxtfuel,maxtmod,maxdmod),NG,NCOMP))
ALLOCATE(dsigd_tr1(NUMOFDXS,max(maxboron,maxtfuel,maxtmod,maxdmod),NG,NCOMP))
ALLOCATE(dsignf1(NUMOFDXS,max(maxboron,maxtfuel,maxtmod,maxdmod),NG,NCOMP))
ALLOCATE(dsigkf1(NUMOFDXS,max(maxboron,maxtfuel,maxtmod,maxdmod),NG,NCOMP)) 
ALLOCATE(dsigs1(NUMOFDXS,max(maxboron,maxtfuel,maxtmod,maxdmod),NG,NCOMP)) 
ALLOCATE(dsigchi1(NUMOFDXS,max(maxboron,maxtfuel,maxtmod,maxdmod),NG,NCOMP))
ALLOCATE(dsigf1(NUMOFDXS,max(maxboron,maxtfuel,maxtmod,maxdmod),NG,NCOMP))
ALLOCATE(dsigt_p31(NUMOFDXS,max(maxboron,maxtfuel,maxtmod,maxdmod),NG,NCOMP))
ALLOCATE(dsigsm1(NG,NUMOFDXS,max(maxboron,maxtfuel,maxtmod,maxdmod),NG,NCOMP))

ALLOCATE(dsigt_a10(NUMOFDXS,max(maxboron,maxtfuel,maxtmod,maxdmod),NG,NCOMP))
ALLOCATE(dsigd_tr10(NUMOFDXS,max(maxboron,maxtfuel,maxtmod,maxdmod),NG,NCOMP))
ALLOCATE(dsignf10(NUMOFDXS,max(maxboron,maxtfuel,maxtmod,maxdmod),NG,NCOMP))
ALLOCATE(dsigkf10(NUMOFDXS,max(maxboron,maxtfuel,maxtmod,maxdmod),NG,NCOMP)) 
ALLOCATE(dsigs10(NUMOFDXS,max(maxboron,maxtfuel,maxtmod,maxdmod),NG,NCOMP)) 
ALLOCATE(dsigchi10(NUMOFDXS,max(maxboron,maxtfuel,maxtmod,maxdmod),NG,NCOMP))
ALLOCATE(dsigf10(NUMOFDXS,max(maxboron,maxtfuel,maxtmod,maxdmod),NG,NCOMP))
ALLOCATE(dsigt_p310(NUMOFDXS,max(maxboron,maxtfuel,maxtmod,maxdmod),NG,NCOMP))
ALLOCATE(dsigsm10(NG,NUMOFDXS,max(maxboron,maxtfuel,maxtmod,maxdmod),NG,NCOMP))

ALLOCATE(ADFACTOR(NG,NCOMP))
ALLOCATE(ADFACTOR0(NG,NCOMP))
ALLOCATE(DADFACTOR(NUMOFDXS,NG,NCOMP))
ALLOCATE(DADFACTOR0(NUMOFDXS,max(maxboron,maxtfuel,maxtmod,maxdmod),NG,NCOMP))

ALLOCATE(OCRDXS(NCOMP,NTYPE2,NG,3)) ! 2014_07_28 . PKB
ALLOCATE(ORIGINADF_CR(NCOMP,NG,3))
ALLOCATE(PPMDADF_CR(NCOMP,NG,3))
ALLOCATE(TMDADF_CR(NCOMP,NG,3))
ALLOCATE(DMDADF_CR(NCOMP,NG,MAXDMOD,3))

! 2014_08_05 . PKB
ALLOCATE(IFMASFUEL(NCOMP))  
ALLOCATE(IFAX(NCOMP))       
! ADDED END

SBURNSTEP(:,:,:)  = 0.d0
SDERSTEP(:,:,:)   = 0.d0
NUCND_T(:,:)      = 0.d0

N2N_T(:)          = 0.d0
SBNUC(:,:)        = 0.d0
ORIGINXS(:,:,:,:) = 0.d0
PPMDXS(:,:,:,:,:)   = 0.d0
TFDXS(:,:,:,:,:)    = 0.d0
TMDXS(:,:,:,:,:)    = 0.d0
DMDXS(:,:,:,:,:)  = 0.d0
MACORIGIN(:,:,:)  = 0.d0
MACDPPM(:,:,:,:)    = 0.d0
MACDTF(:,:,:,:)     = 0.d0
MACDTM(:,:,:,:)     = 0.d0
MACDDM(:,:,:,:)   = 0.d0
BASEADF(:,:,:)    = 1.d0
PPMADF(:,:,:,:)     = 0.d0
TFADF(:,:,:,:)      = 0.d0
TMADF(:,:,:,:)      = 0.d0
DMADF(:,:,:,:)    = 0.d0
ORIGINADf(:,:)    = 1.d0
PPMDADF(:,:,:)      = 0.d0
TFDADF(:,:,:)       = 0.d0
TMDADF(:,:,:)       = 0.d0
DMDADF(:,:,:)     = 0.d0
ADFACTOR(:,:)     = 1.d0
ADFACTOR0(:,:)    = 1.d0
DADFACTOR(:,:,:)  = 0.d0
DADFACTOR0(:,:,:,:)  = 0.d0


dsigt_a1 = 0.d0
dsigd_tr1 = 0.d0
dsignf1 = 0.d0
dsigkf1 = 0.d0 
dsigs1 = 0.d0 
dsigchi1 = 0.d0
dsigf1 = 0.d0
dsigt_p31 = 0.d0
dsigsm1 = 0.d0

dsigt_a10 = 0.d0
dsigd_tr10 = 0.d0
dsignf10 = 0.d0
dsigkf10 = 0.d0 
dsigs10 = 0.d0 
dsigchi10 = 0.d0
dsigf10 = 0.d0
dsigt_p310 = 0.d0
dsigsm10 = 0.d0



! 2014_07_28 . pkb   - CONTORL ROD CROSS-SECTION
BASECRDXS(:,:,:,:)     = 0.d0
PPMCRDXS(:,:,:,:)      = 0.d0
DMCRDXS(:,:,:,:)       = 0.d0
BASECRDXS_H(:,:,:,:,:) = 0.d0
PPMCRDXS_H(:,:,:,:,:)  = 0.d0
DMCRDXS_H(:,:,:,:,:,:) = 0.d0
OCRDXS(:,:,:,:)        = 0.d0
CRDKAPPA(:,:,:)        = 0.d0
BASEADF_CR(:,:,:,:)    = 0.d0
PPMADF_CR(:,:,:,:)     = 0.d0
TMADF_CR(:,:,:,:)     = 0.d0
DMADF_CR(:,:,:,:,:)    = 0.d0
PPMDXS_H(:,:,:,:)      = 0.d0
DMDXS_H(:,:,:,:,:)     = 0.d0
ORIGINADF_CR(:,:,:)    = 0.d0
PPMDADF_CR(:,:,:)      = 0.d0
TMDADF_CR(:,:,:)      = 0.d0
DMDADF_CR(:,:,:,:)     = 0.d0
! added end

! 2014_07_23 . scb
allocate(icm2ica(ncomp))
icm2ica=0
! added end

! 2014_08_04 . scb
allocate(ica2icm(cnum+5))   ! 2014_08_06 . scb
ica2icm=0
! added end

END SUBROUTINE
