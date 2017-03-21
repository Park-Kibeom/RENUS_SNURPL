SUBROUTINE REFCALMAC(ICOMP,INDEX)

USE MASTERXSL
!USE TEMP_BURNUP
USE PARAM

INCLUDE 'GLOBAL.H'
INCLUDE 'BENCH.H'
INCLUDE 'XSEC.H'
INCLUDE 'SRCHPPM.H'
INCLUDE 'CARDS.H'

CHARACTER,INTENT(IN)::INDEX*3
INTEGER,INTENT(IN)::ICOMP
INTEGER::NUM
REAL::FACT,FACT2

FACT  = 3.2992605E-07
FACT2 = PPM * FACT

IF(INDEX.EQ.'TOP') THEN
    NUM=2
    ! 2014_08_06 . scb
    if(naxial.eq.1) stop 'The number of axial reflector types should be 2 when "TOP" reflector is used.'
ELSEIF(INDEX.EQ.'BOT') THEN
    NUM=1
END IF

! 2014_08_06 . scb
RDM(icm2ica(ICOMP)) = RFND(NUM,RH2O) * AH2O / AVOG *1.E3
! added end

RFND(NUM,RB10) = RFND(NUM,RH2O)*FACT2
DO IG=1,NG
    MACAXIAL(NUM,RB10,ABSO,IG)  = AXIAL(NUM,RB10,ABSO,IG)*RFND(NUM,RB10)
    MACAXIAL(NUM,RH2O,ABSO,IG)  = AXIAL(NUM,RH2O,ABSO,IG)*RFND(NUM,RH2O)
    MACAXIAL(NUM,RSTRM,ABSO,IG) = AXIAL(NUM,RSTRM,ABSO,IG)
    MACAXIAL(NUM,RB10,TRAN,IG)  = AXIAL(NUM,RB10,TRAN,IG)*RFND(NUM,RB10)
    MACAXIAL(NUM,RH2O,TRAN,IG)  = AXIAL(NUM,RH2O,TRAN,IG)*RFND(NUM,RH2O)
    MACAXIAL(NUM,RSTRM,TRAN,IG) = AXIAL(NUM,RSTRM,TRAN,IG)
    MACAXIAL(NUM,RB10,SCA,IG)   = AXIAL(NUM,RB10,SCA,IG)*RFND(NUM,RB10)
    MACAXIAL(NUM,RH2O,SCA,IG)   = AXIAL(NUM,RH2O,SCA,IG)*RFND(NUM,RH2O)
    MACAXIAL(NUM,RSTRM,SCA,IG)  = AXIAL(NUM,RSTRM,SCA,IG)
END DO

END SUBROUTINE