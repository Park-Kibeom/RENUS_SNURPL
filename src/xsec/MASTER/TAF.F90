! 2014_08_07 . SCB
  REAL FUNCTION TAF( I1, P1, P2  )
!                            12.11.1985/KOE
!     FILE: TAF.FOR
!  ***********************************************************
!
!  ***************************
!     F U N C T I O N  TAF
!  ***************************
!
!
!
!     0. STAND
!
!        AUSGABESTAND SIEHE DATA - ANWEISUNG
!
!
!
!     1. FUNKTION
!
!        W A S S E R - D A M P F - T A F E L
!
!
!        TAF ERRECHNET DIE        TAF CALCULATES THE PROPERTIES
!        ZUSTANDSVARIABLEN        OF LIGHT WATER  ACCORDING TO IFC 1967
!        FUER 'LEICHTES' WASSER   FORMULATION OF PROPERTIES FOR
!        NACH DER IFC-FORMU-      WATER AND STEAM.
!        LIERUNG VON 1967.
!
!
!        DER ERSTE UBERGABE-      THE FIRST FORMAL PARAMETER I1
!        PARAMETER I1 BESTIMMT    CONTROLS THE PROPERTY AND THE MODE
!        DIE ZU BERECHNENDE       OF CALCULATION. IT MUST BE DETERMINED
!        ZUSTANDSVARIABLE UND     FOR THE ACTUAL CALL ACCORDING TO
!        DEN RECHEN-MODUS NACH    THE FOLLOWING RULES:
!        DEN FOLGENDEN REGELN:
!
!
!        I1= IA+IB                CONTROLPARAMETER FOR PROPERTY
!                                 AND MODE OF CALCULATION;
!                                 I1 POSITIVE FOR LIQUID PROPERTIES,
!                                 I1 NEGATIVE FOR VAPOR PROPERTIES.
!        I1 POSITIV FUER
!        DEN ZUSTAND 'WASSER'.
!
!        I1 NEGATIV FUER DEN
!        DAMPF-ZUSTAND.
!
!
!        IA IST,WIE FOLGT         IA IS DETERMINED
!        AUSZUWAEHLEN:            AS FOLLOWS:
!
!        IA= 1        SPEC. VOLUME       (M3/KG)      V
!        IA= 2        SPEC. ENTHALPY     (KJ/KG)      H
!        IA= 3        SPEC. ENTROPY      (KJ/KG/K)    S
!        IA= 4        THERM. EXPANSION   (1/K)       BETA
!        IA= 5        SPEC. HEAT         (KJ/KG/K)    CP
!        IA= 6        THERM. CONDUCTIVITY(KW/M/K)    LAMBDA
!        IA= 7        VISKOSITY          (KG/M/S)     ETA
!        IA= 8        SURFACE TENSION    (N/M)       SIGMA
!
!
!        IB BESTIMMT DEN          IB DETERMINES THE MODE OF CALCULATION
!        RECHEN-MODUS UND DIE     DEPENDING WHAT PARAMETERS ARE INPUT.
!        EINGANGS-VARIABLEN NACH  IN THE TABLE THE VARIABLES ARE
!        DEN FOLGENDEN FEST-      USED AS FOLLOWS:
!        LEGUNGEN:
!
!        V    ZUSTANDSVARIABLE    PROPERTY
!        P    DRUCK               PRESSURE
!        T    TEMPERATUR          TEMPERATURE
!
!        DER DRUCK IST ALS        THE PRESSURE MUST BE CHOSEN AS
!        ABSOLUT-DRUCK IN         ABSOLUT PRESSURE.
!        N/M2 (PASCAL) ANZUGEBEN.
!        BEISPIEL:                EXAMPLE:
!        10**5 PASCAL  =  1 BAR
!
!        DER MAXIMALE DRUCK       THE MAXIMUM PRESSURE IS 1000 BAR
!        BETRAEGT 1000 BAR BZW.   OR 10**8 PASCAL.
!        10**8 PASCAL.
!
!        DIE TEMPERATUR IST IN    THE TEMPERATURE MUST BE CHOSEN IN
!        GRAD CELSIUS ANZUGEBEN.  CELSIUS DEGREES.
!        BEISPIEL:                EXAMPLE:
!        0.0  0C  =  273.15 K
!
!        DIE MAXIMALE TEMPERATUR  THE MAXIMUM TEMPERATURE IS 800 0C.
!        BETRAEGT 800 0C.
!
!        IB   AUSGABE(OUTPUT) EINGABE(INPUT)        BEISPIEL(EXAMPLE)
!        100  PROPERTY        P,T                   V= TAF(I1,P,T)
!        200  PRESSURE        PROPERTY,T            P= TAF(I1,V,T)
!        300  TEMPERATURE     P,PROPERTY            T= TAF(I1,P,V)
!
!
!        FUER ALLE ZUSTANDS-      IB MUST BE CHOSEN  IB = 100
!        VARIABLEN MIT IA > 3     FOR ALL PROPERTIES WITH  IA > 3.
!        IST AUS KONVERGENZ-
!        GRUENDEN NUR IB = 100
!        ZULAESSIG.
!
!        BEACHTEN SIE ANORDNUNG   PLEASE RECOGNIZE ORDER AND UNITS
!        UND EINHEITEN DER        OF THE ACTUAL PARAMETERS.
!        UEBERGABE-PARAMETER.
!
!        AUS KONVERGENZ-GRUENDEN  IB = 200  CANNOT BE CHOSEN FOR
!        IST  IB = 200  FUER      LIQUID PROPERTIES BECAUSE OF
!        'WASSER'-ZUSTAENDE       CONVERGENCE REASONS.
!        NICHT ZULAESSIG.
!
!        IN DEN FAELLEN           FOR ALL CASES  IB > 100  THE OUTPUT
!        IB > 100  WERDEN DIE     IS DETERMINED ITERATIVLY.
!        RECHEN-WERTE ITERATIV    THEREFORE, MORE COMPUTER TIME
!        BESTIMMT. DESHALB WIRD   IS NEEDED IN SUCH CASES.
!        MEHR RECHENZEIT
!        VERBRAUCHT.
!
!        BEI EINER                IF THE LIQUID TEMPERATURE EXCEEDS
!        WASSER-TEMPERATUR VON    T = 350 0C , TAF IS USED OUTSIDE
!        T = 350 0C  WIRD DIE     THE RANGE OF VALIDITY.
!        GUELTIGKEITS-GRENZE
!        FUER DIE TAFEL UEBER-
!        SCHRITTEN.
!
!        BEI EINER                IF THE LIQUID TEMPERATURE EXCEEDS
!        WASSER-TEMPERATUR VON    APPROXIMATELY  T = 359 0C ,
!        ETWA  T = 359 0C  WIRD   TAF TERMINATES WITH A
!        DIE RECHNUNG MIT DER     NEGATIVE SQRT.
!        FEHLER-MELDUNG
!        'NEGATIVES WURZEL-
!        ARGUMENT' ABGEBROCHEN.
!
!        GESAETTIGTE ZUSTAENDE    SATURATED PROPERTIES ARE CALCULATED
!        MUESSEN MIT DEM SAETTI-  USING SATURATION PRESSURE AND
!        GUNGS-DRUCK ODER DER     SATURATION TEMPERATURE AS ACTUAL
!        SAETTIGUNGS-TEMPERATUR   PARAMETERS.
!        ALS UEBERGABE-PARAMETER  PLEASE USE THE SUBROUTINES TSA
!        GERECHNET WERDEN.AUS     OR PSA LIKE THEY ARE USED IN TAF
!        STETIGKEITS-GRUENDEN     BECAUSE OF AVOIDING DISCONTINUITIES.
!        SOLLTEN, WIE IM UNTER-
!        PROGRAMM TAF, DIE
!        UNTERPROGRAMME  TSA
!        BZW. PSA  BENUTZT
!        WERDEN.
!
!        DER AUFRUF VON TAF       DO NOT CALL TAF WITH CONDITIONS
!        UEBER DIE SAETTIGUNGS-   BEYOND SATURATION.
!        BEDINGUNGEN HINAUS,      TWO PHASE MIXTURE CONDITIONS
!        ALSO : FORMEL FUER DIE   MUST BE IDENTIFIED AND HANDLED
!               WASSER-PHASE AN-  IN THE CALLING CODE APPLYING
!               GEWANDT IN DER    SATURATED PROPERTIES.
!               DAMPF-PHASE UND
!               UMGEKEHRT,
!        IST NICHT ZULAESSIG.
!        ZWEI-PHASEN-GEMISCH-
!        ZUSTAENDE SIND NUR MIT
!        AUF DEN SAETTIGUNGS-
!        BEREICH ZUGESCHNITTENEN
!        RECHEN-PROGRAMMEN SINN-
!        VOLL BERECHENBAR.
!
!
!
!     2. VERSORGUNG
!
!        PARAMETER
!            I1
!            P1
!            P2
!
!
!
!     3. ERGEBNIS
!
!        FUNKTIONS-WERT TAF
!
!
!
!     4. ANWENDER - MELDUNGEN UND - FEHLERMELDUNGEN
!
!        2140  KEINE RECHENFUNKTION ZUR ANGEGEBENEN STEUERZAHL VORHANDEN
!        2141  PHYSIKALISCHE EINGANGSGROESSE NEGATIV
!        2142  ITERATION DIVERGENT
!
!
!
!     5. VERWENDETE UP'S :
!
!        REAL FUNCTION PSA                BERECHNUNG DER
!                                         SAETTIGUNGS-TEMPERATUR
!        REAL FUNCTION TSA                BERECHNUNG DES
!                                         SAETTIGUNGS-DRUCKS
!
!
!
!     6. SPRACHNORM/-UMFANG, BENUTZTE SPRACHERWEITERUNGEN
!        FORTRAN 77
!
!
!
!  ************************************************************
!
!
!  ************************************************************
!
!                      D A T E N K A T A L O G
!
!  ************************************************************
!
!     NAME    TYP      DIM      BEDEUTUNG
!
!
!     AHD     R*8               KRIT. VOL. * KRIT. DRUCK (DAMPFPHASE)
!     AHW     R*8               KRIT. VOL. * KRIT. DRUCK (WASSERPHASE)
!     ARG     R*8               HILFSVARIABLE FUER DIE WAERMELEITFAEIG.
!     ASD     R*8               KRIT. VOL. * KRIT. DRUCK /
!                               KRIT. TEMP.(DAMPFPHASE)
!     ASW     R*8               KRIT. VOL. * KRIT. DRUCK /
!                               KRIT. TEMP.(WASSERPHASE)
!     B       R*8    3          HILFSFELD FUER ITERATION
!     BEZ     C*20              DYN. ZEICHENSTRING FUER FEHLERMELDUNG
!     BEZBL   C*20              BLANK-VORBESETZUNG FUER BEZ
!     BIJ     R*8    6,5        KOEFFIZIENTEN FUER DYN. VISKOSITAET
!     C       R*8    3          HILFSFELD FUER ITERATION
!     CPK     R*8               KRITISCHER DRUCK
!     CT      R*8               KRITISCHE TEMPERATUR
!     CVD     R*8               KRITISCHES VOLUMEN (DAMPFPHASE)
!     CVW     R*8               KRITISCHES VOLUMEN (WASSERPHASE)
!     DD      R*8               SCHRITTWEITE FUER ITERATION
!     DELLAM  R*8               DELTA-LAMBDA-TERM FUER WAERMELEITFAEIGK.
!     DELTST  R*8               DELTA-T-STERN-TERM FUER WAERMELEITF.
!     DLAMB0  R*8               DELTA-LAMBDA0-TERM FUER WAERMELEITF.
!     DLAMBQ  R*8               DELTA-LAMBDAQUER-TERM FUER WAERMELEITF.
!     ETA0    R*8               ETA0-TERM FUER DYN. VISKOSITAET
!     FENR    C*4               FEHLERNUMMER
!     FEORT   C*4               FEHLERORT
!     I       I*4               LAUFVARIABLE
!     I1      I*4               STEUERPARAMETER
!     IA      I*4               STEUERZAHL FUER DIE DIV. ZUSTAENDE
!     IB      I*4               STEUERZAHL FUER AUSGABEWERT DER TAF
!     J       I*4               LAUFVARIABLE
!     JJ      I*4               STEUERZAHL FUER AUSGABEWERT DER TAF
!     K       I*4               LAUFVARIABLE
!     KK      I*4               STEUERZAHL FUER WEITERE ITERATION(C)
!     LL      I*4               STEUERZAHL FUER WEITERE ITERATION(C)
!     NITERA  I*4               ANZAHL DER ITERATIONS-SCHRITTE
!     NN      I*4               STEUERZAHL FUER WEITERE ITERATION
!                               (STEUERUNG DES SCHRITT-ABLAUFS)
!     P       R*8            E  REDUZIERTER DRUCK
!     P1      R*4               PARAMETER
!     P2      R*4               PARAMETER
!     POH     R*8               HILFSVARIABLE F. SPEZ. ENTROPIE(DAMPF)
!     POLH    R*8               HILFSVARIABLE F. SPEZ. ENTROPIE(WASSER)
!                                                SPEZ. ENTHALPIE(DAMPF)
!                                                SPEZ. ENTHALPIE(WASSER)
!     QSCHLA  R*8               Q-SCHLANGE-TERM F. WAERMELEITFAEHIKEIT
!     R       R*8               ALLGEMEINE RECHENVARIABLE (IM ALLGEMEI-
!                               NEN IDENTISCH MIT TAF)
!     RHOM1   R*8               (V-STERN/V - 1)-TERM F. DYN. VISKOSITAET
!     RSCHLA  R*8               R-TERM FUER WAERMELEITFAEHIGKEIT
!     SA      R*8               HILFSVARIABLE FUER WAERMELEITFAEHIGKEIT
!     SB      R*8               HILFSVARIABLE FUER WAERMELEITFAEHIGKEIT
!     T       R*8            E  REDUZIERTE TEMPERATUR
!     TAF     R*4               FUNKTIONS-WERT
!     TAFALT  R*8               ALT-WERT VON TAF ZUR BERECHNUNG VON
!                               SPEZ. WAERME  UND
!                               THERMISCHER EXPANSION
!     TAFEL   R*8               FUNKTIONS-WERT ALS  R*8-KONSTANTE
!     TC      R*8               P**5  FUER DAMPF
!     TETA    R*8               REDUZIERTE TEMPERATUR F. WAERMELEITF.
!                               UND DYN. VISKOSITAET
!     TETAM1  R*8               (TETA - 1)-TERM FUER DYN. VISKOSITAET
!     TFUN    R*8               REDUZIERTE TEMPERATUR FUER
!                               OBERFLAECHENSPANNUNG
!     UA      R*8               HILFSVARIABLE
!     UB      R*8               HILFSVARIABLE (DAMPF)
!     UBL     R*8               BETA-L-TERM FUER DAMPF-PHASE
!     UC      R*8               HILFSVARIABLE
!     UD      R*8               HILFSVARIABLE (DAMPF)
!     UE      R*8               HILFSVARIABLE (DAMPF)
!     UF      R*8               HILFSVARIABLE (DAMPF)
!     UG      R*8               HILFSVARIABLE (DAMPF)
!     UH      R*8               HILFSVARIABLE (DAMPF)
!     UI      R*8               HILFSVARIABLE (DAMPF)
!     UP      C*8               UNTERPROGRAMMNAME (UEBERGABEPARAMETER)
!     UX      R*8               X-EXPONENTIAL-TERM FUER DAMPFPHASE
!     V       R*8               HILFSVARIABLE FUER PARAMETERUEBERGABE
!                               UND ITERATION(ALTWERT)
!     VA      R*8               HILFSVARIABLE
!     VB      R*8               HILFSVARIABLE
!     VC      R*8               HILFSVARIABLE
!     VD      R*8               HILFSVARIABLE
!     VE      R*8               HILFSVARIABLE
!     VERSIO  C*4               AUSGABESTAND
!     VERSNR  C*4               AUSGABESTAND (UEBERGABEPARAMETER)
!     VF      R*8               HILFSVARIABLE (WASSER)
!     VH      R*8               HILFSVARIABLE (WASSER)
!     W       L*4               LOGISCHE VAR. .TRUE. : WASSER-PHASE
!                                             .FALSE.: DAMPF-PHASE
!     WA      R*8               HILFSVARIABLE
!     WB      R*8               HILFSVARIABLE
!     WC      R*8               HILFSVARIABLE
!     WD      R*8               HILFSVARIABLE
!     WE      R*8               HILFSVARIABLE (DAMPF)
!     WF      R*8               HILFSVARIABLE (DAMPF)
!     WG      R*8               HILFSVARIABLE
!     WH      R*8               HILFSVARIABLE
!     WY      R*8               Y-TERM FUER WASSER-PHASE
!     WZ      R*8               Z-TERM FUER WASSER-PHASE
!     X       R*8    2          FELD FUER ITERATION:
!                               X(1) = T (REDUZIERTE TEMPERATUR)
!                               X(1) = P (REDUZIERTER DRUCK)
!     X1      R*8               HILFSVARIABLE FUER ITERATION
!     X3      R*8               HILFSVARIABLE FUER ITERATION
!     XA      R*8               HILFSVARIABLE (DAMPF)
!     XB      R*8               HILFSVARIABLE (DAMPF)
!     XC      R*8               HILFSVARIABLE (DAMPF)
!     ZVZ     R*8               VORZEICHEN VON I1 (PARAMETER)
!                               ZVZ = +1  : WASSER-PHASE
!                               ZVZ = -1  : DAMPF-PHASE
!     ZW1     R*8               HILFSVARIABLE F. WAERMELEITFAEHIGKEIT,
!                               SPEZ. WAERME UND THERMISCHE EXPANSION
!     ZW2     R*8               HILFSVARIABLE F. DYN. VISKOSITAET
!     ZW3     R*8               HILFSVARIABLE F. DYN. VISKOSITAET
!
!
!
!  ************************************************************
!
!
!     IMPLIZITE TYPANWEISUNGEN
!  ************************
!
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
!
!
!     EXPLIZITE TYPANWEISUNGEN
!  ************************
!
      REAL            P1
      REAL            P2
      REAL            FLOAT
      REAL            PSA
      REAL            SNGL
      REAL            TSA
!
      LOGICAL W
!
      CHARACTER * 20 BEZ
      CHARACTER * 20 BEZBL
      CHARACTER * 4  FENR
      CHARACTER * 4  FEORT
      CHARACTER * 8  UP
      CHARACTER * 4  VERSIO
      CHARACTER * 4  VERSNR
!
!
!
!
!     FELDANWEISUNGEN
!     ***************
!
      DIMENSION B(3)
      DIMENSION BIJ(6,5)
      DIMENSION C(3)
      DIMENSION X(2)
!
!
!     SPEICHERBEREICHSANWEISUNGEN FUER LOKALE COMMON - BEREICHE
!  *********************************************************
!
      COMMON  / IFEHL / VERSNR,FENR,UP,BEZ
!
!
!     INCLUDE - ANWEISUNG FUER GLOBALE COMMON - BEREICHE
!     **************************************************
!
!     KEINE
!
!
!     AEQUIVALENZANWEISUNGEN
!     **********************
!
      EQUIVALENCE (T,X(1))
      EQUIVALENCE (P,X(2))
!
!     ANFANGSWERTANWEISUNGEN
!     **********************
!
!
!     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!                            AUSGABESTAND
!CC   DATA VERSIO         /'A001'/
!CC   OPTION (PRODUCT_ID = 'A001')
!
!     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
      DATA BEZBL  / '                    ' /
!
      DATA        CVW,        AHW,             ASW    &
          /0.00317D0,  70.1204D0,  0.1083275143D0/
!
      DATA        CVD,        AHD,             ASD    &
          /0.00317D0,  70.1204D0,  0.1083275143D0/
!
      DATA      CPK,       CT    &
          /221.2D5,  647.3D0/
!
      DATA ((BIJ( I , J ), I = 1, 6 ), J = 1, 5 )    &
          / 0.501938D0 ,  0.162888D0 , -0.130356D0 ,    &
            0.907919D0 , -0.551119D0 ,  0.146543D0 ,    &
            0.235622D0 ,  0.789393D0 ,  0.673665D0 ,    &
            1.207552D0 ,  0.0670665D0, -0.084337D0 ,    &
           -0.274637D0 , -0.743539D0 , -0.959456D0 ,    &
           -0.687343D0 , -0.497089D0 ,  0.195286D0 ,    &
            0.145831D0 ,  0.263129D0 ,  0.347247D0 ,    &
            0.213486D0 ,  0.100754D0 , -0.032932D0 ,    &
           -0.0270448D0, -0.0253093D0, -0.0267758D0,    &
           -0.0822904D0,  0.0602253D0, -0.0202595D0/
!
!
!
!     FORMATANWEISUNGEN
!     *****************
!
!CC 1 FORMAT (I4)
!
!
!
!  ***************************************************************
!
!
!
!     AUSFUEHRBARE ANWEISUNGEN
!     ************************
!
!
!
!
!                            FEHLER-VORBESETZUNG
      BEZ    = BEZBL
      FEORT  = '   1'
      FENR   = '    '
!
!                            ABFRAGE, OB P1 ODER P2 FEHLERHAFT
!                            FEHLERORT -------------------------->>>   1
      IF(P1.LE.0.0 .OR. P2.LE.0.0) GOTO 90100
!
!                            ABFRAGE, OB STEUERZAHL FEHLERHAFT
      IF(I1.GE.-303 .AND. I1.LE.-301   .OR.    &
        I1.GE.-203 .AND. I1.LE.-201   .OR.    &
        I1.GE.-107 .AND. I1.LE.-101   .OR.    &
        I1.GE. 101 .AND. I1.LE. 108   .OR.    &
        I1.GE. 301 .AND. I1.LE. 303       ) GOTO 5
      FEORT  = '   2'
!                            FEHLERORT -------------------------->>>   2
      GOTO 90000
!
    5 CONTINUE
!
!
!
!
!
!
!                            INVESTIGATION OF I1 AND
!                            -----------------------------
!                            INITIAL VALUES FOR ITERATIONS
!                            -----------------------------
!
!                            VORBESETZUNG STEUERZAHL JJ=2 BZW.IB=200
      JJ     = 2
!                            VORBESETZUNG W FUER FLUESSIGE PHASE
      W      = .TRUE.
!                            VORBESETZUNG NN=1 (KEINE ITERATION)
      NN     = 1
!                            VORBESETZUNG FUER BETA,CP
      TAFALT = 0.D0
!                            VORBESETZUNG W FUER DAMPF-PHASE
      IF (I1.LT.0) W = .FALSE.
!
      I      = IABS(I1)
      P      = DBLE(P1)/CPK
      T      = (DBLE(P2)+273.15D0)/CT
      ZVZ    = DBLE(FLOAT(ISIGN(1,I1)))
      I      = I - 100
!
!                            VERZWEIGUNG FUER IB=100,200,300
      IF (I.LT.100) GO TO 30
      IF (I.GT.200) GO TO 10
!
!
      I      = I - 100
      IA     = I
      IB     = 200
      FEORT  = '   3'
!
!
!                            ABFRAGEN,OB P1,P2 GUELTIG
!                            FEHLERORT -------------------------->>>   3
      IF(P2.GT. 800.0) GOTO 90200
      IF(IA .EQ. 1  .AND.  P1.GT.   1.0    .OR.    &
        IA .EQ. 2  .AND.  P1.GT.6000.0    .OR.    &
        IA .EQ. 3  .AND.  P1.GT.  10.0        ) GOTO 90200
      P = 1.D0
      V = DBLE(P2)
!                            TEMPERATUR KLEINER T-KRIT
      IF (T.LE.1.D0) GOTO 7
    6 V = DBLE(P1)
!                            VORBESETZUNG SCHRITTWEITE
      DD = 0.6D0
      GO TO 20
!
    7 STOP 'FUNCTION TSA DOES NOT EXIST'
!                            FEHLERABFRAGE FUER TSA
      IF(FENR .NE. '    ') GOTO 90500
      GOTO 6
!
!
   10 CONTINUE
      I      = I - 200
      IA     = I
      IB     = 300
      FEORT  = '   4'
!
!                            ABFRAGEN,OB P1,P2 GUELTIG
!                            FEHLERORT -------------------------->>>   4
      IF(P1.GT.1.0E8 ) GOTO 90200
      IF(IA .EQ. 1  .AND.  P2 .GT.   1.0   .OR.    &
        IA .EQ. 2  .AND.  P2 .GT.6000.0   .OR.    &
        IA .EQ. 3  .AND.  P2 .GT.  10.0       ) GOTO 90200
!
      T      = 1.D0
      V      = DBLE(P1)
!                            DRUCK KLEINER P-KRIT
      IF (P.LE.1.D0) GOTO 16
   15 V      = DBLE(P2)
!                            SONDERBEHANDLUNG FLUESSIGE PHASE
!                            TEMPERATUR IN DER NAEHE VON  T-KRIT
      IF (T.GE.0.9822D0 .AND. W) T = .0876D0*P + 0.907D0
!                            VORBESETZUNG SCHRITTWEITE
      DD     = 1.5D0
!                            STEUERZAHL JJ=1 BZW. IB=300
      JJ     = 1
      GOTO 20
   16 STOP 'FUNCTION PSA DOES NOT EXIST'
!                            FEHLERABFRAGE FUER PSA
      IF(FENR .NE. '    ') GOTO 90500
      GOTO 15
!
   20 CONTINUE
!                            ITERATION EINLEITEN
!                            SCHRITT-ZAEHLER VORBESETZEN
      NITERA = 1
!                            STEUERZAHL FUER ITERATIONS-START VORBESETZ.
      NN     = 2
!                            SCHRITTWEITE FUER FLUESSIGE ZUSTAENDE
!                            VORBESETZEN
      IF (W) DD = 1.D0/DD
      GOTO 35
!
!
   30 CONTINUE
      IA     = I
      IB     = 100
      FEORT  = '   5'
!
!                            ABFRAGE,OB P1,P2 GUELTIG
!                            FEHLERORT -------------------------->>>   5
      IF(P1.GT.1.0E8   .OR.  P2.GT. 800.0) GOTO 90200
   35 CONTINUE
!                            SONDERFAELLE BETA,CP
!                            BERECHNUNG VON ZWEI FUNKTIONSWERTEN
!                            MIT EINEM DELTA-T VON 5.0D-4
      IF (IA.EQ.4 .OR. IA.EQ.5) T = T - 5.0D-4
!
!
!
   40 CONTINUE
!                            VERZWEIGUNG ZU DEN DAMPF-ZUSTAENDEN
      IF (.NOT.W) GO TO 180
!                            SONDERFALL SIGMA
      IF (IA.EQ.8) GO TO 70
!
!
!
!
!
!
!                            COMMON EQUATIONS FOR ALL LIQUID PROPERTIES
!                            ------------------------------------------
!
!
      WA = T*T
      WC = WA**4
      WD = WC*WC
      WY = (-5.362162162D-4/WC-0.8438375405D0)*WA + 1.D0
      WG = (0.65371543D0-T)**9*242.1647003D0
      WH = WD*WA*T + 115.D-8
      VA = WC*WA*T + 15108.D-9
      VB = WA + 0.14188D0
      VC = 1269716088.D-19/WH/WH
      VD = 6.047626338D-14*P**3/(WD*WA*WA*T)
      VE = WD*T*12.93441934D0
      VF = 7.002753165D0 + P
      VH = (1.105710498D-9*P+2.17402035D-8)*P + 2.074838328D-7
      UA = 1.308119072D-5*P*P
      WB = 1.72D0*WY*WY - T*1.4684556978D-1 + P*.0995171774D0
!                            KONVERGENZ-ABFRAGE: WURZEL DES  Z-TERMS
!                                                KLEINER/GLEICH O.
      IF (WB.GE.0.0D0) GO TO 50
      UC = WB
      WB = 0.D0
      IF (NITERA.LT.3 .AND. NN.GT.1) GO TO 50
!
!
!                            ITERATION DIVERGENT
!                            FEHLERORT -------------------------->>>   6
      FEORT  = '   6'
      GOTO 90200
!
!
!
   50 CONTINUE
      WZ = DSQRT(WB) + WY
!                            KONVERGENZ-ABFRAGE UEBER  Z-TERM
      IF (WZ.GE.1.D-8) GO TO 60
      UC = WZ
      WZ = 1.D-8
      IF (NITERA.LT.3 .AND. NN.GT.1) GO TO 60
!
!
!                            ITERATION DIVERGENT
!                            FEHLERORT -------------------------->>>   7
      FEORT  = '   7'
      GOTO 90200
!
!
!
   60 CONTINUE
      UC = 7.982692717D0/WZ**0.294117647D0
!
!
!
   70 CONTINUE
!                            VERZWEIGUNG FUER FLUESSIGE ZUSTAENDE
      GO TO (90, 100, 170, 90, 100, 90, 90, 80), IA
!
!
!
!
!                            PROPERTIES OF LIQUID LIGHT WATER
!                            --------------------------------
!
!
!
   80 CONTINUE
!
!                            SURFACE TENSION IN N/M2
!                            -----------------------
!
!
      TFUN  = (374.136D0-DBLE(P2))/647.286D0
      TAFEL = 235.8D-3*(TFUN**1.256D0)*(1.D0-0.625D0*TFUN)
!
!
!
      GOTO 99999
!
!
!
   90 CONTINUE
!
!                            SPEC. VOLUME OF LIQUID IN  M3/KG
!                            --------------------------------
!
!
      R = (((P*7.2467345D-9+1.383225552D-7-VH)/    &
       VA+(0.204D0-T)*UA)*3.D0+UC*    &
       0.0497585887D0+(0.02284279054D0*T+    &
       (3.D0/VF**4-2.995284926D-4)*VB*    &
       VE+VD*4.D0+1.52241179D-3)*T-0.02616571843D0+    &
       WG*(0.65371543D0-T)+VC*WH)*CVW
!
!
!                            SONDERFALL BETA
      IF (IA.EQ.4) GO TO 110
!                            SONDERFALL LAMBDA
      IF (IA.EQ.6) GO TO 130
!                            SONDERFALL ETA
      IF (IA.EQ.7) GO TO 140
!                            ITERATIONSABFRAGE
      GO TO 220
!
!
!
  100 CONTINUE
!
!                            SPEC. ENTHALPY OF LIQUID IN  KJ/KG
!                            ----------------------------------
!
!
      R = (1.340540541D-3/WC-0.70319795D0)*WA
      POLH = ((((((((16138.168904D0*T-99269.72482D0)*T+270670.12452D0)*    &
       T-429542.08335D0)*T+437564.7096D0)*T-    &
       297071.4308D0)*T+134665.55478D0)*    &
       T-39412.86787D0)*T+6824.687741D0)*T - 542.2063673D0
      R = (((WZ*0.586206897D0-WY*1.4166666667D0+R)*WZ-1.728D0*WY*    &
       R+7.342278489D-2*T)*UC+POLH-((WC*WA*T*12.D0+1.5108D-5)/VA/VA*    &
       VH-VD*T*21.D0-(WD*WA*T*20.D0+1.15D-6)*VC-(T*9.D0+0.65371543D0)*    &
       WG+2.284279054D-2*WA+2.616571843D-2-UA*0.204D0)*P+(1.D0/VF**3+P*    &
       2.995284926D-4)*(VB*17.D0+WA+WA)*VE*T)*AHW
!
!
!                            SONDERFALL CP
      IF (IA.EQ.5) GOTO 110
!                            ITERATIONSABFRAGE
      GOTO 220
!
!
!
  110 CONTINUE
!
!                            SPEC. HEAT CP OF LIQUID AND VAPOR
!                            ---------------------------------
!                            IN KJ/KG/K
!                            ----------
!
!                            THERMAL EXPANSION IN  1/K
!                            -------------------------
!
!
      IF (TAFALT.NE.0.D0) GO TO 120
      TAFALT = R
      T      = T + 1.D-3
      GO TO 40
  120 CONTINUE
      ZW1    = (R+TAFALT)*0.5D0
!                            AUSGABE CP
      TAFEL  = (R-TAFALT)*1.D3/CT
!                            AUSGABE BETA
      IF (IA.EQ.4) TAFEL = TAFEL/ZW1
!
!
!
      GOTO 99999
!
!
!
  130 CONTINUE
!
!                            THERMAL CONDUCTIVITY IN KW/M/K
!                            ------------------------------
!
!                            T.MINAMIYAMA U. J.YATA:"A FORMULATION
!                            OF THERMAL CONDUCTIVITY OF WATER"
!                            SUBMITTED TO SPECIAL COMITTEE OF I.A.P.S.
!                            KYOTO, SEPTEMBER 1976.
!                            LAST REVISION: JUNE 1977.
!                            RANGE OF VALIDITY:
!                            TEMPERATURE : 0- 800 0C     (0-1500 0C )
!                            PRESSURE    : 0-1000 BAR    (0-3000 BAR)
!                            VALID FOR LIQUID AND STEAM
!
      TETA   = T
      R      = 1.D0/(R*317.7D0)
      SB     = DSQRT(TETA)
      DLAMB0 = (((-4.22464D-3*TETA+1.56146D-2)*TETA+2.99621D-2)*    &
                 TETA+0.0102811D0)*SB
      SA     = R + 2.39219D0
      DLAMBQ = DEXP(-0.171587D0*SA*SA)*1.06D0 + 0.400302D0*R - 0.39707D0
      DELTST = DABS(TETA-1.D0) + 0.00308976D0
      ZW1    = DELTST**0.6D0
      QSCHLA = 0.0822994D0/ZW1 + 2.D0
      RSCHLA = QSCHLA*TETA + 1.D0
      ZW1    = 10.0932D0/ZW1
      IF (TETA.GE.1.D0) ZW1 = 1.D0/DELTST
      ARG = -4.11717D0 *SB*TETA - 6.17937D0/(R*R*R*R*R)
      IF (ARG.LE.-6.D2) ARG = -6.D2
      SA     = R**1.8D0
      DELLAM = (0.0701309D0*TETA**(-10)+0.0118520D0)*SA*    &
               DEXP((1.D0-SA*R)*0.642857D0) +    &
              0.00169937D0*ZW1*R**QSCHLA*DEXP(QSCHLA/RSCHLA*    &
               (1.D0-R**RSCHLA)) - 1.02D0*DEXP(ARG)
      TAFEL = (DLAMB0+DLAMBQ+DELLAM)*1.0D-3
!
!
!
      GOTO 99999
!
!
!
  140 CONTINUE
!
!                            DYNAMIC VISCOSITY IN KG/M/S
!                            ---------------------------
!
!                            ACCORDING TO A  A. ALEXANDROV RELEASE
!                            ON THE DYNAMIC VISCOSITY OF WATER SUBSTANCE
!                            - AUG.1ST,1975,APPENDIX C
!                            RANGE OF VALIDITY:
!                            TEMPERATURE :
!                            0- 800 0C     (0- 560 0C )  (0-  100 0C )
!                            PRESSURE    :
!                            0-1000 BAR    (0-3500 BAR)  (0-10000 BAR)
!                            VALID FOR LIQUID AND VAPOR
!
      TETA   = 647.27D0/(DBLE(P2)+273.15D0)
      TETAM1 = TETA - 1.D0
      RHOM1  = 1.D0/(R*0.317763D3) - 1.D0
      ETA0   = (((-3.6744D-3*TETA+1.05287D-2)*TETA+1.77624D-2)*TETA+    &
                  1.81583D-2)*DSQRT(TETA)
      ZW2    = 0.D0
      ZW3    = 1.D0
      DO 160 I=1,6
        SA   = 0.D0
        DO 150 J=1,5
          K    = 6 - J
          SA   = SA*RHOM1 + BIJ(I,K)
  150   CONTINUE
        ZW2 = ZW3*SA + ZW2
        ZW3 = ZW3*TETAM1
  160 CONTINUE
      TAFEL = DEXP(ZW2/(R*0.317763D3))/ETA0*1.D-6
!
!
!
      GOTO 99999
!
!
!
  170 CONTINUE
!
!                            SPEC. ENTROPY OF LIQUID IN  KJ/KG/K
!                            -----------------------------------
!
!
      POLH = (((((((18155.440017D0*T-113451.11408D0)*T+315781.81194D0)*    &
       T-515450.50002D0)*T+546955.887D0)*T-396095.24112D0)*T+    &
       201998.33217D0)*    &
       T-78825.73574D0)*T + 20966.66205D0 + DLOG(T)*6824.687741D0
      R = (((WZ*0.5787037037D0-WY)*(3.2172972972D-3/WC-    &
       1.687675081D0)*0.72D0*    &
       T+0.07342278489D0)*UC+POLH+(1.D0/VF**3+P*2.995284926D-4)*    &
       (VB*9.D0+WA)*    &
       2.D0*VE+(VD*20.D0+(WD*VC*19.D0-WC*VH*11.D0/VA/VA)*WA+WG*    &
       10.D0-4.568558108D-2*T-1.52241179D-3+UA)*P)*ASW
!
!
!                            ITERATIONSABFRAGE
      GO TO 220
!
!
!
!
!
!
  180 CONTINUE
!
!                            COMMON EQUATIONS FOR ALL VAPOR  PROPERTIES
!                            ------------------------------------------
!
!
      UX = DEXP((1.-T)*0.763333333333D0)
      UC = UX**3
      UE = UC**3*UX
      UF = UE*UC*UX
      UG = UF*UC*UX
      UH = UF*UE
      WA = (UE*6.670375918D-2+1.388983801D0)*UC
      VA = WA*3.D0 + 6.670375918D-1*UE*UC
      WB = UG*1.6780208656D-1 + (UX*5.229341786D-2-6.746878906D-2)*UX
      VB = UG*2.8526354715D0 + UX*5.229341786D-2*UX + WB
      WC = (0.2364644526D0*UE+UG)*1.3562756712D0
      VC = (UG*1.085020537D0+WC)*10.D0
      WD = (-UH*UX*6.753673384D0-UF)*0.35390143216D0
      UI = UH*UC
      VD = (-UH*UX*1.877962965D0+WD)*14.D0
      WE = UH*1.037510561D0+UF*UG*2.9790258045D0-UI*UX*2.5796516865D0
      VE = (UF*UG*0.9930086015D0-UI*UX*0.42994194775D0+WE)*24.D0
      WF = (UX*0.4762441084D0-0.39468696528D0)*UX*UE
      WG = (UH*2.8987292793D0-UG)*0.29047190005D0
      WH = (UH*11.474849789D0+UF)*3.4261311894D-3
      TC = P**5
      XC = UF*0.4006073948D0*TC
      XB = UG*UX*8.636081627D-2*TC
      XA = (-0.8532322921D0*UI+0.3460208861D0)*UI
      UD = P/(XC+P)
      XC = (-XC*14.D0*UD/P+11.D0)*WF + UE*UX*UX*0.4762441084D0
      UC = 1.D0/(XB+1.D0)
      XB = (-XB*19.D0*UC+18.D0)*WG + UH*5.051996409D0
      UB = 1.D0/(TC*XA*P+1.D0)
      XA = ((46.07454378D0*UI-9.342563925D0)*UI*UB*TC*P+14.D0)*WH +    &
       0.39314340756D0*UH
      UI = (((((UX*523.5718623D0-2693.088365D0)*UX+5745.984054D0)*    &
       UX-6508.211677D0)*UX+4126.607219D0)*UX-1388.522425D0)*UX +    &
       193.6587558D0
      UH = (((((3141.4311738D0*UX-13465.441825D0)*UX+22983.936216D0)*    &
       UX-19524.635031D0)*UX+8253.214438D0)*UX-1388.522425D0)*UX
      UG = T*0.7633333333D0
      UBL= (19.31380707D0*T-34.17061978D0)*T + 15.74373327D0
      UE = (386.2761414D0*T-341.7061978D0)/UBL
      UA = (P/UBL)**10
!
!
!                            VERZWEIGUNG FUER DAMPF-ZUSTANDSGROESSEN
      GO TO (190, 200, 210, 190, 200, 190, 190), IA
!
!
!
!
!                            PROPERTIES OF LIGHT WATER VAPOR
!                            -------------------------------
!
!
!
  190 CONTINUE
!
!                            SPEC. VOLUME OF VAPOR IN  M3/KG
!                            -------------------------------
!
!
      R = (((((-WH*UB*UB*P-WG*UC*UC-WE)*P-WF*UD*UD-WD)*P-WC)*P-WB)*    &
       P-WA+4.260321148D0*T/P+11.D0*UI*UA)*CVD
!
!
!
!                            SONDERFALL LAMBDA
      IF (IA.EQ.6) GO TO 130
!                            SONDERFALL ETA
      IF (IA.EQ.7) GO TO 140
!                            SONDERFALL BETA
      IF (IA.EQ.4) GO TO 110
!                            ITERATIONSABFRAGE
      GO TO 220
!
!
!
  200 CONTINUE
!
!                            SPEC. ENTHALPY OF VAPOR IN  KJ/KG
!                            ---------------------------------
!
!
      POLH = ((((-18.01781976D0*T+91.82563266D0)*T-30.36678102D0)*    &
       T+1180.5465453D0)*T+2002.686163D0)/AHD
      R = (((((((-XA*UG-WH)*UB*.8333333333D0*P-    &
       (XB*UG+WG)*UC-VE*UG-WE)*P*    &
       .8D0-(XC*UG+WF)*UD-VD*UG-WD)*P*.25D0-(VC*UG+WC)*0.333333333D0)*    &
       P-(VB*UG+WB)*.5D0)*P-VA*UG-WA+(T*UE*UI+UI+UG*UH)*UA)*P+POLH)*AHD
!
!
!                            SONDERFALL CP
      IF (IA.EQ.5) GO TO 110
!                            ITERATIONS-ABFRAGE
      GO TO 220
!
!
!
  210 CONTINUE
!
!                            SPEC. ENTROPY OF VAPOR IN  KJ/KG/K
!                            ----------------------------------
!
!
      POH = DLOG(T)*16.83599274D0+54.38923329D0+((-0.34260728232D0*    &
       T+1.9643135091D0)*T-0.8661325668D0)*T
      R = (((((((-XA*UB*P*.8333333333D0-XB*UC-VE)*P*.8D0-XC*UD-VD)*P*    &
       .25D0-VC*.3333333333D0)*P-VB*.5D0)*P-VA+UH*UA)*.7633333333D0+    &
       UE*UI*UA)*P-DLOG(P)*4.260321148D0+POH)*ASD
!
!
!
!
!
!
!                            PROPERTY ITERATION, IF NEEDED
!                            -----------------------------
!
!
  220 CONTINUE
!                            STEUER-VERZWEIGUNG FUER ITERATION
!                            NN=1  : KEINE ITERATION
!                            NN=2  : START ITERATION (SCHRITT 1)
      GO TO (310, 230, 260, 280), NN
!
!
  230 CONTINUE
!                            SCHRITT 2 DER ITERATION
      NN     = 3
      C(3)   = R - V
  240 CONTINUE
      C(1)   = C(3)
  250 CONTINUE
      B(1)   = B(3)
      B(3)   = X(JJ)
      X(JJ)  = DMIN1(X(JJ)*DD,4.52D0)
      IF (JJ.EQ.1) T = (1.D0-DD)*0.42D0+ T
!                            WEITERER ITERATIONS-SCHRITT
      GO TO 290
!
!
  260 CONTINUE
!                            SCHRITT 3 DER ITERATION
      C(3)   = R - V
      IF (C(1)*C(3).GT.0.D0) GO TO 270
      NN     = 4
      DD     = 0.5D0/DD + 0.5D0
      C(2)   = 2.D0*C(1)
!
      GO TO 250
!
!
  270 CONTINUE
      IF ((C(3)-C(1))*C(3).LT.0.D0) GO TO 240
      DD     = 1.D0/DD/DD
!
      GO TO 240
!
!
  280 CONTINUE
      LL     = 1
      KK     = 3
      IF (C(2)*C(3).GT.0.D0) LL = 2
      IF (C(1)*C(2).GT.0.D0) KK = 2
      IF (DABS(C(KK)).GT.DABS(C(LL))) LL = KK
      B(LL)  = X(JJ)
      C(LL)  = R - V
!
!                            ABFRAGE: ITERATIONS-ENDE
      IF (DABS(R/V-1.D0)                .LT.1.D-5    .OR.    &
         NITERA                        .GT.   20    .OR.    &
         DABS((C(2)-C(1))*(C(2)-C(3))) .LT.1.D-50       )GOTO 300
      X3     = (B(2)-B(1))/(C(2)-C(1))
      X1     = (B(3)-B(1))/(((B(3)-B(2))/(C(2)-C(3))+X3)*C(2)+B(1)-B(3))
      X(JJ)  = C(1)*X3*X1 + B(1)
!
  290 CONTINUE
!                            ZAEHLINDEX FUER ITERATION ERHOEHEN,
!                            WEITERER ITERATIONSSCHRITT
      NITERA = NITERA + 1
      GO TO 40
!
!                            ITERATIONSENDE
  300 CONTINUE
!                            FUER JJ=1 BZW. IB=300  TEMPERATUR AUSGEBEN
      R      = T*CT - 273.15D0
!                            FUER JJ=2 BZW. IB=200  DRUCK AUSGEBEN
      IF (JJ.EQ.2) R = P*CPK
  310 CONTINUE
      TAFEL  = R
!
!
!
      GOTO 99999
!
!
!
!
!
!
!                            FEHLERBEHANDLUNG:
!                            -----------------
!
!
!                            KEINE RECHENFUNKTION ZUR ANGEGEBENEN
!                            STEUERZAHL VORHANDEN
90000 CONTINUE
!CC   VERSNR     = VERSIO
      UP         = 'TAF     '
      FENR       = '2140'
      BEZ(10:13) = FEORT
      TAFEL      = 1.0D0
      GOTO 99999
!
!                            PHYSIKALISCHE EINGANGSGROESSE NEGATIV
90100 CONTINUE
!CC   VERSNR     = VERSIO
      UP         = 'TAF     '
      FENR       = '2141'
!CC   WRITE(BEZ(1:4)  ,1) I1
      BEZ(10:13) = FEORT
      TAFEL      = 1.0D0
      GOTO 99999
!
!                            ITERATION DIVERGENT
90200 CONTINUE
!CC   VERSNR     = VERSIO
      UP         = 'TAF     '
      FENR       = '2142'
!CC   WRITE(BEZ(1:4)  ,1) I1
      BEZ(10:13) = FEORT
      TAFEL      = 1.0D0
      GOTO 99999
!
!                            FEHLER IM UNTERPROGRAMM
90500 CONTINUE
      TAFEL    = 1.0D0
      GOTO 99999
!
!
!
!
!
!
!                            ENDE UNTERPROGRAMM
!                            ------------------
!                            RUNDUNG
99999 TAF = SNGL(TAFEL*1.00000005D0)
!
!
!
      RETURN
      END
