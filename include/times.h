! file name = times.h - cpu time data
	  real :: tinit,tth,ttotal,txsec,tppr
      common /times/ tinit,tth,ttotal,txsec,tppr

	  real :: tsteady,ttransient
	  common /times2 / tsteady, ttransient ! 2015_06_22 . scb
