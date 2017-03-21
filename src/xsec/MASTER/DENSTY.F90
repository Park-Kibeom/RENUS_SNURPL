! 2014_08_07 . SCB
  REAL FUNCTION DENSTY(TM,P0)
      REAL :: TM, P0, T0
!     ******************************************************************
!     TEMPERATUR IN GRAD CELSIUS, DRUCK IN BAR, DICHTE IN KG/DM**3
!     TEMPERATUR 1 BIS 350 GRAD C, DRUCK 0.2 BIS 300 BAR
!     ******************************************************************
      T0 = TM
      IF(T0.GT.350.0) T0=350.0
!
!     BERECHNUNG DER SAETTIGUNGSTEMPERATUR T0 (GR C) BEIM DRUCK P0 (BAR)
!
      P = P0
      DENSTY=0.001/TAF(101,P*1.E5,T0)
      RETURN
      END
