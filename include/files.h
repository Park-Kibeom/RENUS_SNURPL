! file name = files.h - input and ouput file names
      character*80 caseid,filename(40),localfn,args
      common /filenames/caseid,filename,localfn,args
      common /iodev/ibase,ixslib,ioutp,irstin,irstout,isum,icrho,ifbrho
      common /iodev/ io1, io2, io3, io4, io5, io6, io7, io8, io9,io10, &
                    io11,io12,io13,io14,io15,io16,io17,io18,io19,io20, &
                    io21,io22,io23,io24,io25,io26,io27,io28,io29,io30

	  logical  :: filexist   ! 2014_10_06 . scb

	  logical :: flagout(8)   ! 2014_10_15 . scb
	  common / ioflag / flagout   ! 2014_10_15 . scb 
! 1 : outl
! 2 : plt
! 3 : peak
! 4 : rod
! 5 : dt
! 6 : flux
! 7 : xsf
! 8 : nu


! filename(1) / io1 : plt           
! filename(2) / io2 : dt
! filename(3) / io3 : rod
! filename(4) / io4 : peak
! io8 : out
! io9 : outl
! filename(9) / io19 : rsi
! filename(10) / io20 : rst
! filename(17) / io17 : FBdata writing
! filename(21) / io21 : READ_DOPL
! filename(22) / io22 : READ_TMDM
! 