SUBROUTINE READCOMP(INDEV,CNUM)

USE MASTERXSL
!USE TEMP_BURNUP
USE PARAM

INCLUDE 'GLOBAL.H'
INCLUDE 'BENCH.H'
INCLUDE 'XSEC.H'
INCLUDE 'SRCHPPM.H'
INCLUDE 'CARDS.H'

LOGICAL::FIRST=.TRUE.   ! 2014_07_28 . pkb
INTEGER,INTENT(IN)::INDEV,CNUM
INTEGER::LIN1,LIN2,LIN3
INTEGER::IOPT(4),NBLTAB(5)
REAL::ND

FIRST = .TRUE.   ! 2014_08_20 . pkb
SBFLAG(CNUM) = .FALSE.

! 2014_10_02 . PKB ADDED
NXY2 = NXY
NZ2 = NZ
!ADDEND


LIN1=NBURN(CNUM)/6
IF(MOD(NBURN(CNUM),6).GT.0) THEN
    LIN1=LIN1+1
END IF
LIN2=NDER(CNUM)/6
IF(MOD(NDER(CNUM),6).GT.0) THEN
    LIN2=LIN2+1
END IF
LIN3=NDER(CNUM)/6
IF(MOD(NDER(CNUM),6).GT.0) THEN
    LIN3=LIN3+1
END IF

! PREVIOUSLY READ LINE 7 WHICH CONTAIN REFERENCE VALUES AND THE NUMBER OF VALUES
CALL SKIP(INDEV,7)

! READ THE BURNUP STEP OF EACH COMPOSITION
CALL READLIN(NBURN(CNUM),LIN1,BURNSTEP(CNUM,:),INDEV)

! READ THE BRANCH BURNUP STEP OF EACH COMPOSITION
CALL READLIN(NDER(CNUM),LIN2,DERSTEP(CNUM,:),INDEV)

! READ THE MODERATOR STEP OF EACH COMPOSITION
CALL READLIN(NDMOD(CNUM),LIN3,MODSTEP(CNUM,:),INDEV)

BACKSPACE(INDEV)
BACKSPACE(INDEV)
IF(NADF(CNUM).GT.0) THEN
    FLAGADFMAS = .TRUE.   ! 2014_07_28 . pkb
    if(iver.eq.4) then   ! 2014_08_01 . scb
      DO IA=1,NADF(CNUM)
          DO IG=1,NG
              CALL READLIN(NBURN(CNUM),LIN1,BASEADF(CNUM,:,IG),INDEV)
          END DO
          DO IM=1,NBORON(CNUM)
            DO IG=1,NG
              CALL READLIN(NDER(CNUM),LIN2,PPMADF(CNUM,:,IG,IM),INDEV)
            END DO
          END DO
          DO IM=1,NTFUEL(CNUM)
            DO IG=1,NG
              CALL READLIN(NDER(CNUM),LIN2,TFADF(CNUM,:,IG,IM),INDEV)
            END DO
          END DO
          DO IM=1,NDMOD(CNUM)
              DO IG=1,NG
                  CALL READLIN(NDER(CNUM),LIN2,DMADF(CNUM,:,IG,IM),INDEV)
              END DO
          END DO
      END DO
    else
      stop 'iver should be 4 for ADF treatment !'
    endif
    
END IF

! READ NUCLIDE WISE CROSS SECTION
BACKSPACE(INDEV)
DO WHILE(.TRUE.)
    READ(INDEV,'(A512)',END=4000) ONELINE
    IF(ONELINE(1:4).EQ.'COMP') GOTO 4000
    IF(PROBE.EQ.DOT .OR. PROBE.EQ.SLASH) EXIT
    IF(PROBE.EQ.'-' .OR. PROBE.EQ.BANG .OR. ONELINE.EQ.BLANK .OR. IFNUMERIC(ONELINE)) CYCLE

! READ A NUMBER DENSITY AND XS OPTIONS WHICH ARE IOPT AND NBLTAB
    READ(ONELINE(6:),*) ND,(IOPT(I),I=1,4),(NBLTAB(J),J=1,5)

! FIND A SPECIFIED NUCLIDE ID NUMBER
    DO IN=1,NUCNUM
        IF(ONELINE(1:4).EQ.NU2ID(IN)) EXIT
    END DO
    NUCND_T(CNUM,IN)=ND

! READ A KAPPA VALUE USED TO CALCULATE A POWER
    !IF(IN.LE.MAXHEAVY .OR. IN.GE.30) READ(INDEV,*) (KAPPA_T(CNUM,IN,I),I=1,NG) 
    IF(NBLTAB(1) .NE. 0) READ(INDEV,*) (KAPPA_T(CNUM,IN,I),I=1,NG)  ! 2014_07_28 . pkb
    IF(IOPT(2).NE.0) THEN
        CALL READLIN(NBURN(CNUM),LIN1,N2N_T(:),INDEV)
    END IF
    SELECT CASE(IOPT(1))

! READ A NORMAL TYPE OF CROSS SECTION
        CASE(0)
            DO IT=1,NTYPE1
                IF(IT.NE.NTYPE1) THEN
                    NUM=NG
                ELSE 
                    NUM=NG-1
                ENDIF
! DO NOT READ ANY XS
                IF(NBLTAB(IT).EQ.0) THEN
                    CONTINUE

! READ THE XS SETS CONTAINING THE BASE CONDITION AND DERIVATIVES TERM
                ELSEIF(NBLTAB(IT).EQ.15) THEN
                    DO IG=1,NUM
                        CALL READLIN(NBURN(CNUM),LIN1,BASEXS(CNUM,IT,:,IN,IG),INDEV)
                    END DO
                    DO IM=1,NBORON(CNUM)
                      DO IG=1,NUM
                        CALL READLIN(NDER(CNUM),LIN2,PPMXS(CNUM,IT,:,IN,IG,IM),INDEV)
                      END DO
                    END DO
                    DO IM=1,NTFUEL(CNUM)
                      DO IG=1,NUM
                        CALL READLIN(NDER(CNUM),LIN2,TFXS(CNUM,IT,:,IN,IG,IM),INDEV)
                      END DO
                    END DO
                    DO IM=1,NDMOD(CNUM)
                        DO IG=1,NUM
                            CALL READLIN(NDER(CNUM),LIN2,DMXS(CNUM,IT,:,IN,IG,IM),INDEV)
                        END DO
                    END DO
! 2014_07_28 . pkb                    
                    IF(IVER.EQ.3) THEN
                      stop  ' stop in line 128 in readcomp subroutine ' 
                        DO IG=1,NUM
                            !CALL READLIN(NDER(CNUM),LIN2,TFXS(CNUM,IT,:,IN,IG),INDEV)
                            !CALL READLIN(NDER(CNUM),LIN2,TMXS(CNUM,IT,:,IN,IG),INDEV)   ! 2014_08_01 . scb
                        END DO
                    END IF
! added end                    

! READ THE XS SETS CONTAINING THE BASE CONDITION ONLY
                ELSEIF(NBLTAB(IT).EQ.16) THEN
                    DO IG=1,NUM
                        CALL READLIN(NBURN(CNUM),LIN1,BASEXS(CNUM,IT,:,IN,IG),INDEV)
                    END DO                   
                END IF
            END DO

! 2014_07_28 . pkb
! CASE 1 IS SEPARATED BY MASTER XSD VERSION
! IT SHOULD BE SUMMED UP WITH PROPER CONDITION
            
! SPECIAL BURNUP STEP IS CONTAINED
        CASE(1)          
! 2014_07_28 . pkb          
          IF(IVER .EQ. 3) THEN
! CONTROL ROD CROSS-SECTION IS READ FROM HERE
            DO IT=1,NTYPE1
                IF(IT.NE.NTYPE1) THEN
                    NUM=NG
                ELSE 
                    NUM=NG-1
                ENDIF
! DO NOT READ ANY XS
                IF(NBLTAB(IT).EQ.0) THEN
                    CONTINUE
                ELSEIF(NBLTAB(IT).EQ.15) THEN
                    DO IG=1,NUM
                        CALL READLIN(NDER(CNUM),LIN2,BASECRDXS(CNUM,IT,:,IG),INDEV)
                    END DO
                    DO IG=1,NUM
                        CALL READLIN(NDER(CNUM),LIN2,PPMCRDXS(CNUM,IT,:,IG),INDEV)
                    END DO
                    DO IG=1,NUM
                        CALL READLIN(NDER(CNUM),LIN2,DMCRDXS(CNUM,IT,:,IG),INDEV)
                    END DO
                ELSEIF(NBLTAB(IT).EQ.16) THEN
                    DO IG=1,NUM
                        CALL READLIN(NDER(CNUM),LIN2,BASECRDXS(CNUM,IT,:,IG),INDEV)
                    END DO
                ENDIF
            END DO
          ELSEIF(IVER .EQ. 4) THEN
! added end
              
            SBFLAG(CNUM) = .TRUE.
            SBNUC(CNUM,IN)=IN
! READ THE BURNUP STEP OF EACH COMPOSITION
            CALL READLIN(NBURN(CNUM),LIN1,SBURNSTEP(CNUM,:,IN),INDEV)

! READ THE BRANCH BURNUP STEP OF EACH COMPOSITION
            CALL READLIN(NDER(CNUM),LIN2,SDERSTEP(CNUM,:,IN),INDEV)

            DO IT=1,NTYPE1
                IF(IT.NE.NTYPE1) THEN
                    NUM=NG
                ELSE 
                    NUM=NG-1
                ENDIF
! DO NOT READ ANY XS
                IF(NBLTAB(IT).EQ.0) THEN
                    CONTINUE

! READ THE XS SETS CONTAINING THE BASE CONDITION AND DERIVATIVES TERM
                ELSEIF(NBLTAB(IT).EQ.15) THEN
                    DO IG=1,NUM
                        CALL READLIN(NBURN(CNUM),LIN1,BASEXS(CNUM,IT,:,IN,IG),INDEV)
                    END DO
                    DO IM=1,NBORON(CNUM)
                    DO IG=1,NUM
                        CALL READLIN(NDER(CNUM),LIN2,PPMXS(CNUM,IT,:,IN,IG,IM),INDEV)
                    END DO
                    END DO
                    DO IM=1,NTFUEL(CNUM)
                    DO IG=1,NUM
                        CALL READLIN(NDER(CNUM),LIN2,TFXS(CNUM,IT,:,IN,IG,IM),INDEV)
                    END DO
                    END DO
                    DO IM=1,NDMOD(CNUM)
                        DO IG=1,NUM
                            CALL READLIN(NDER(CNUM),LIN2,DMXS(CNUM,IT,:,IN,IG,IM),INDEV)
                        END DO
                    END DO

! READ THE XS SETS CONTAINING THE BASE CONDITION ONLY
                ELSEIF(NBLTAB(IT).EQ.16) THEN
                    DO IG=1,NUM
                        CALL READLIN(NBURN(CNUM),LIN1,BASEXS(CNUM,IT,:,IN,IG),INDEV)
                    END DO                    
                END IF
            END DO
          END IF ! 2014_07_28 . PKB ADDED IF STATEMENT

! CONTROL ROD XS IS CONTAINED
        CASE(3)
          if(iver.eq.3) stop 'In version 3, iopt(1) should not be 3 ! '   ! 2014_08_01 . scb
! 2014_07_28 . PKB          
            DO L=1,3
! 2014_08_20 . pkb              
                IF(FIRST) THEN
                    BACKSPACE(INDEV)
                    FIRST = .FALSE.
                END IF
                READ(INDEV,*) (CRDKAPPA(CNUM,IG,L), IG=1,NG)
                IF(FLAGADFMAS) THEN
                    DO IG=1,NG
                        CALL READLIN(NDER(CNUM),LIN2,BASEADF_CR(CNUM,:,IG,L),INDEV)
                    END DO
                    DO IG=1,NG
                        CALL READLIN(NDER(CNUM),LIN2,PPMADF_CR(CNUM,:,IG,L),INDEV)
                    END DO
                    DO IM=1,NDMOD(CNUM)
                        DO IG=1,NG
                            CALL READLIN(NDER(CNUM),LIN2,DMADF_CR(CNUM,:,IG,IM,L),INDEV)
                        END DO
                    END DO
                END IF       
! added end                
                DO IT=1,NTYPE1
                    IF(IT.NE.NTYPE1) THEN
                        NUM=NG
                    ELSE 
                        NUM=NG-1
                    ENDIF

! 2014_08_20 . pkb                    
                    !IF(FLAGADFMAS .AND. FIRST) THEN
                    !    CALL SKIP(INDEV,48)
                    !    FIRST=.FALSE.
                    !END IF
! commented end                    
! DO NOT READ ANY XS
                    IF(NBLTAB(IT).EQ.0) THEN
                        CONTINUE
                    ELSEIF(NBLTAB(IT).EQ.15) THEN
                        DO IG=1,NUM
                            CALL READLIN(NDER(CNUM),LIN2,BASECRDXS_H(CNUM,IT,:,IG,L),INDEV)
                        END DO
                        DO IG=1,NUM
                            CALL READLIN(NDER(CNUM),LIN2,PPMCRDXS_H(CNUM,IT,:,IG,L),INDEV)
                        END DO
                        DO IM=1,NDMOD(CNUM)
                            DO IG=1,NUM
                               CALL READLIN(NDER(CNUM),LIN2,DMCRDXS_H(CNUM,IT,:,IG,L,IM),INDEV)
                            END DO
                        END DO
                    ELSEIF(NBLTAB(IT).EQ.16) THEN
                        DO IG=1,NUM
                            CALL READLIN(NDER(CNUM),LIN2,BASECRDXS_H(CNUM,IT,:,IG,L),INDEV)
                        END DO
                    ENDIF
                END DO
            END DO
!          
            !GOTO 400
! ADDED END            
        CASE DEFAULT
            GOTO 400
    END SELECT
400 CONTINUE
END DO

4000 CONTINUE

END SUBROUTINE