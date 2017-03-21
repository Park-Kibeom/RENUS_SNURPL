Subroutine ReadND(ixy,iz)
USE MASTERXSL
!USE TEMP_BURNUP
USE PARAM
!USE TIMER
!USE DECUSPING1N,  ONLY : XSSET ! 2012_09_28 . SCB
use readN2A
use deplmod

INCLUDE 'GLOBAL.H'
INCLUDE 'TIMES.H'
INCLUDE 'XSEC.H'
INCLUDE 'GEOM.H'
INCLUDE 'SRCHPPM.H'
INCLUDE 'THFDBK.INC'
!INCLUDE 'THFUEL.INC'
INCLUDE 'THCNTL.INC'
INCLUDE 'THGEOM.INC'
INCLUDE 'THOP.INC'
! ADDED IN ARTOS VER. 0.2 . 2012_07_11 BY SCB.
INCLUDE 'FILES.H'
INCLUDE 'FFDM.H'

INTEGER,INTENT(IN)::ixy,iz
REAL :: KLO, KHI, a

!-------------------------------------
!        EXTRAPOLATION START
!-------------------------------------

icomp = icompnumz(iz,ixy)
ic = icompnumz(iz,ixy)

IF (BUP.GT.BURNSTEP(IC,NBURN(IC))) THEN
    KLO = NBURN(IC)-2
    AFRAC1 = 0.
    AFRAC3 = (BUP-BURNSTEP(IC,NBURN(IC)-1))/(BURNSTEP(IC,NBURN(IC))-BURNSTEP(IC,NBURN(IC)-1))
    AFRAC2 = 1.-AFRAC3
    RETURN
ELSE IF (BUP.LT.BURNSTEP(IC,1)) THEN
    KLO = 1.
    AFRAC2 = (BUP-BURNSTEP(IC,1))/(BURNSTEP(IC,2)-BURNSTEP(IC,1))
    AFRAC1 = 1.-AFRAC2
    AFRAC3 = 0.
    RETURN
END IF

!-------------------------------------
!        INTERPOLATION START
!-------------------------------------

IF(BUP .LT. BURNSTEP(IC,2)) THEN
    KLO = 1
    KHI = 2
    GOTO 280
ELSE IF(BUP .GT. BURNSTEP(IC,NBURN(IC)-2)) THEN
    KLO = NBURN(IC)-2
    KHI = NBURN(IC)-1
    GOTO 280
END IF

KLO = 2
KHI = NBURN(IC)-2

1000 IF(KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(BURNSTEP(IC,K).GT.BUP) THEN
           KHI=K
        ELSE
           KLO=K
        ENDIF
        GOTO 1000       
     ENDIF 

280 CONTINUE
H=BURNSTEP(IC,KHI)-BURNSTEP(IC,KLO)
IF(H.EQ.0) THEN
   WRITE(*,*) 'SOMETHING IS WRONG'
ENDIF

!if (nburn(ic)>=3) then    ! add 2017.2.1 jjh    

AB1 = BUP
A1  = BURNSTEP(IC,KLO)
A2  = BURNSTEP(IC,KHI)
A3  = BURNSTEP(IC,KHI+1)
 
A12 = A1 - A2
A21 = -A12
A23 = A2 - A3
A32 = -A23
A31 = A3 - A1
A13 = -A31
XA1 = AB1-A1
XA2 = AB1-A2
XA3 = AB1-A3
AFRAC1 = XA2*XA3/(A12*A13)
AFRAC2 = XA1*XA3/(A21*A23)
AFRAC3 = XA1*XA2/(A31*A32)

DO in=1,NUCNUM
    atomn(ixy,iz,in)= afrac1*assm(icomp)%nuclei(in)%ND(klo)&
                    + afrac2*assm(icomp)%nuclei(in)%ND(klo+1)&
                    + afrac3*assm(icomp)%nuclei(in)%ND(klo+2)
END DO

End Subroutine