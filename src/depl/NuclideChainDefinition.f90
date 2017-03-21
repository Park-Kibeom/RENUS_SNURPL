SUBROUTINE NuclideChainDefinition

USE DEPLMOD
USE MASTERXSL

!***********************************************************************
!*                                                                     *
!*    ROUTINE NAME      : NuclideChainDefinition                                       *
!*    ROUTINE TYPE      : SUBROUTINE                                   *
!*                                                                     *
!***********************************************************************
!*                                                                     *
!*    DESCRIPTION       : NUCLIDE CHAIN DESCRIPTION                    *
!*                        IN THIS ROUTINE COMMON BLOCK HCHAIN, FCHAIN  *
!*                        DCYCON AND FYIELD ARE DEFINED                *
!*                                                                     *
!***********************************************************************
!
!     NCHNH     : NUMBER OF HEAVY NUCLIDE CHAINS
!     NCHNF     : NUMBER OF FISSION PRODUCT CHAINS
!     NTNH      : TOTAL NUCLIDES IN HEAVY CHAIN
!     NTNF      : TOTAL NUCLIDES IN FISSION PRODUCT CHAIN
!     NDCY      : NUMBER OF NUCLIDES DECAY CONSIDERED
!     NFCNT     : NUMBER OF NUCLIDES GENERATE FISSION PRODUCT
!     NFP       : NUMBER OF FISSION PRODUCTS
! 
INTEGER NUCSRCH
DIMENSION ICHNH(NCHNH),ICHNF(NCHNF)
DIMENSION LCHAINH(NTNH),LPTYPE(NTNH),LCHAINF(NTNF),LDPLCT(NTNH)
DIMENSION IDCY(NDCY),DECAY(NDCY)
DIMENSION IFH(NFCNT),IFF(NFP),FYLD(NFP,NFCNT)
CHARACTER*4 LCHAINH,LCHAINF,IDCY,IFH,IFF
!
!***********************************************************************
!                                                                      *
!     HEAVY NUCLIDE CHAIN DESCRIPTION                                  *
!                                                                      *
!***********************************************************************
!
!       1 : U235  - U236  - NP237 - PU238 - PU239 - PU240 - PU241 - PU242 - AM243                     
!       2 : U238  - NP239 - PU239 - PU240 - PU241 - PU242 - AM243
!       3 : U238  - NP237 - PU238 - PU239 - PU240 - PU241 - PU242 - AM243
!
      DATA ICHNH   /     9,     7,     8/
      DATA LCHAINH /'U235','U236','NP37','PU48','PU49','PU40','PU41','PU42','AM43'&
                   ,'U238','NP39','PU49','PU40','PU41','PU42','AM43' &
                   ,'U238','NP37','PU48','PU49','PU40','PU41','PU42','AM43'/

!     LPTYPE  1 : PRODUCED BY CAPTURE OF I-1 NUCLIDE
!             2 : PRODUCED BY BETA DECAY OF I-1 NUCLIDE
!             3 : PRODUCED BY (N,2N) OF I-1 NUCLIDE

      DATA LPTYPE  /     1,     1,     1,     1,     1,     1,     1,     1,     1&
                   ,     1,     1,     2,     1,     1,     1,     1 &
                   ,     1,     3,     1,     1,     1,     1,     1,     1/

!     LDPLCT 0 : DUPLICATED CHAIN, INITIAL DENSITY SET TO 0.0
!            1 : USE INITIAL NUMBER DENSITY
!            2 : NEGLECT DEPLETED DENSITY

      DATA LDPLCT  /     1,     1,     1,     1,     0,     0,     0,     0,     0&
                   ,     1,     1,     1,     1,     1,     1,     1&
                   ,     2,     0,     0,     0,     0,     0,     0,     0/
!
!***********************************************************************
!                                                                      *
!     FISSION PRODUCT CHAIN DESCRIPTION                                *
!                                                                      *
!***********************************************************************
!
!       1 : PM149 - SM149
!       2 : I135  - XE135

      DATA ICHNF   /     2,     2/
      DATA LCHAINF /'PM49','SM49'&
                   ,'I135','XE35'/
!
!***********************************************************************
!                                                                      *
!     DECAY CONSTANT ASIGN (/SEC)                                      *
!                                                                      *
!***********************************************************************
!
      DATA IDCY /'NP39'    ,'PU41'    ,'PM49'    ,'I135'    ,'XE35'/

      DATA DECAY/3.441E-06 ,1.536E-09 ,3.626E-06 ,2.924E-05 ,2.100E-05/
!
!***********************************************************************
!                                                                      *
!     FISSION YIELD MATRIX                                             *
!                                                                      *
!***********************************************************************

      DATA IFH /'U235','U236','NP37','U238','PU49','PU40','PU41','PU42'/
      DATA IFF /'PM49','SM49','I135','XE35','LFP'/

      DATA FYLD/0.01067, 0.0, 0.06298, 0.00242, 1.0,&
                0.01067, 0.0, 0.06298, 0.00242, 1.0,&
                0.01067, 0.0, 0.06298, 0.00242, 0.0,&
                0.01608, 0.0, 0.06827, 0.00028, 1.0,&
                0.01239, 0.0, 0.06447, 0.01152, 1.0,&
                0.01239, 0.0, 0.06447, 0.01152, 1.0,&
                0.01524, 0.0, 0.07068, 0.00231, 1.0,&
                0.01524, 0.0, 0.07068, 0.00231, 1.0/

!***********************************************************************
!     HEAVY NUCLIDE CHAIN DEFINE

NHCHN = NCHNH
DO 100 I=1,NHCHN
    J1=0
    NIHCHN(I) = ICHNH(I)
    IF(I.LE.1) GOTO 60
    DO K=1,I-1
        J1=J1+ICHNH(K)
    ENDDO
60  J1=J1+1
    J2=J1+ICHNH(I)-1
    DO 100 J=J1,J2
        K=J-J1+1
        IHCHN(I,K)=NUCSRCH(LCHAINH(J),MNUCL)
        IPTYP(I,K)=LPTYPE(J)
        IDPCT(I,K)=LDPLCT(J)
100 CONTINUE
   
!     FISSION PRODUCT CHAIN DEFINE

NFCHN = NCHNF
DO 200 I=1,NFCHN
    J1=0
    NIFCHN(I) = ICHNF(I)
    IF(I.LE.1) GOTO 160
    DO K=1,I-1
        J1=J1+ICHNF(K)
    ENDDO
160 J1=J1+1
    J2=J1+ICHNF(I)-1
    DO 200 J=J1,J2
        IFCHN(I,J-J1+1)=NUCSRCH(LCHAINF(J),MNUCL)
200 CONTINUE

!     DECAY CONSTANT ASSIGN

DO I=1,NUCNUM
    DCNST(I)=0.0
ENDDO
DO I=1,NDCY
    K = NUCSRCH(IDCY(I),MNUCL)
    DCNST(K)  =  DECAY(I)
ENDDO

!      FISSION YIELD MATRIX ASSIGN

NHFCNT=NFCNT
DO I=1,NUCNUM
    DO J=1,NFCNT
        FYFRC(I,J)=0.0
    ENDDO
ENDDO

DO J=1,NFCNT
    IHFCNT(J) = NUCSRCH(IFH(J),MNUCL)
ENDDO

!  --- FOR YIELD FRACTION MODIFICATION, IF IXSMOD = 1
ixsmod = 0
iyld = 0
DO I=1,NFP
    IF (IXSMOD .GE. 1 .AND. IYLD .GE. 1) THEN
        IF     (I .EQ. 1) THEN
           FACTMOD = FYPM 
        ELSEIF (I .EQ. 3) THEN
           FACTMOD = FYIOD 
        ELSEIF (I .EQ. 4) THEN
           FACTMOD = FYXEN
        ELSE
           FACTMOD = 1.0
        ENDIF
    ELSE
       FACTMOD = 1.0
    ENDIF

    K = NUCSRCH(IFF(I),MNUCL)
    
    DO J=1,NFCNT
        FYFRC(K,J) = FYLD(I,J) * FACTMOD
    ENDDO            
ENDDO

DO I= 1,NFP
    DO J= 1,NFCNT
       FYLDP(I,J)=FYLD(I,J)
    ENDDO
ENDDO

END