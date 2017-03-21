SUBROUTINE READREFL_N2A_xsc(INDEV, icomp,rf)

USE MASTERXSL
!USE TEMP_BURNUP
USE PARAM

INCLUDE 'GLOBAL.H'
INCLUDE 'BENCH.H'
INCLUDE 'XSEC.H'
INCLUDE 'SRCHPPM.H'
INCLUDE 'CARDS.H'
INCLUDE 'geom.H'

INTEGER,INTENT(IN)::INDEV, icomp
CHARACTER, intent(in)::RF*6
INTEGER::NBLTAB(5), INDEX=0, i=0, g, bu
real(8) :: sum1

  if(rect==.true.) then
    nrdirdf=4
  else
    nrdirdf=6
  endif
    

   IF(RF.EQ.'CORNER') THEN
    INDEX=1
   ELSEIF(RF.EQ.'EDGE  ') THEN
    INDEX=2
   ELSEIF(RF.EQ.'NOOK  ') THEN
    INDEX=3
   ENDIF     
    
    NBLTAB(2)=1
    SELECT CASE(rf)
        CASE('BOTTOM')
            NBLTAB(1)=1
!            CALL READAXIAL(INDEV,RBOTTOM)
        CASE('CORNER')
            NBLTAB(1)=2
            DO IG=1, NG
              RADIAL(INDEX,TRAN,IG)=SIGTR(IG,ICOMP)
              RADIAL(INDEX,SCA,IG)=SIGS(IG,ICOMP)
              RADIAL(INDEX,ABSO,IG)=SIGA(IG,ICOMP)              
            END DO
!            CALL READRADIAL(INDEV,RCORNER,NBLTAB(2))  ! 2014_08_07 . SCB
        CASE('EDGE  ')
            NBLTAB(1)=3            
            DO IG=1, NG
              RADIAL(INDEX,TRAN,IG)=SIGTR(IG,ICOMP)
              RADIAL(INDEX,SCA,IG)=SIGS(IG,ICOMP)
              RADIAL(INDEX,ABSO,IG)=SIGA(IG,ICOMP)              
            END DO
            !CALL READRADIAL(INDEV,REDGE)
!            CALL READRADIAL(INDEV,REDGE,NBLTAB(2))  ! 2014_08_07 . SCB
        CASE('TOP   ')
            NBLTAB(1)=4
            CALL READAXIAL(INDEV,RTOP)
        CASE('NOOK  ')
            NBLTAB(1)=5
            DO IG=1, NG
              RADIAL(INDEX,TRAN,IG)=SIGTR(IG,ICOMP)
              RADIAL(INDEX,SCA,IG)=SIGS(IG,ICOMP)
              RADIAL(INDEX,ABSO,IG)=SIGA(IG,ICOMP)              
            END DO            
            !CALL READRADIAL(INDEV,RNOOK)
!            CALL READRADIAL(INDEV,RNOOK,NBLTAB(2))  ! 2014_08_07 . SCB
        CASE DEFAULT
    END SELECT
BACKSPACE(INDEV)

do g=1, ng
  sum1=0.d0; i=0
  do irdir=1,nrdirdf 
    if (sigadf(irdir,g, icomp)/=1) then
      sum1=sum1 + sigadf(irdir,g, icomp)
      i=i+1
    end if
  end do
  if (i==0) then
    baseadf(icomp, :, g)=1
  else 
    baseadf(icomp, :, g)=sum1/real(i)
  end if
end do

DO G=1,NG
    ORIGINADF(ICOMP,G) = BASEADF(IComp,1,G)
END DO

DO L=1,3
    DO G=1,NG
        ORIGINADF_CR(ICOMP,G,L) = BASEADF_CR(IComp,1,G,L)
    END DO
END DO

do l=1,nxy
  xsadf(:,:,l,k)=sigadf(:,:,icomp)
end do

END SUBROUTINE READREFL_N2A_xsc