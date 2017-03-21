    subroutine default
!
      use param
      use allocs
!
      include 'global.h'
      include 'xsec.h'
      include 'times.h'
      include 'files.h'
      include 'itrcntl.h'
      include 'geom.h'
      include 'nodal.h'
      include 'pinpr.h'
      include 'cntl.h'
      include 'srchppm.h'
      include 'pow.h'
      include 'thgeom.inc'
      include 'thop.inc'
      include 'thcntl.inc'
      include 'thfuel.inc'
      include 'perturb.inc'
      include 'ff.h'
      include 'frzfdbk.inc'
      include 'trancntl.inc'
      include 'thexpan.inc' ! added in ARTOS ver. 0.2 (THERMAL EXPANSION). 2012_07_05 by SCB 
      include 'thlink.inc' ! added in ARTOS ver. 0.2 . 2012_07_25 by SCB 
      include 'xesm.h'    ! 2013_09_26 . scb

! global.h
      data iflfr /FALSE/    ! added in ARTOS ver. 0.2 (LFR). 2012_07_05 by SCB 
	
! xsec.h
      data nprec /6/
      data ifcompname /FALSE/ ! added in ARTOS ver. 0.2 (TRINX). 2012_07_05 by SCB 
      data decusp, initxsec /FALSE,FALSE/   ! 2012_08_23 . scb
      data lfb /FALSE/ !2016_02_11 . BYS
! cntl.h
      data prtscrn,tracinp,transient /TRUE,FALSE,FALSE/
      data qtransient /FALSE/
      data nthread /1/

! itrcntl.h defaults
! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
!      data ifcmfd/TRUE/
      data flagcmfd/TRUE/
      data noutmax,ninmax,nupmax,ninitcmfd,nmaxcmfd2g,nmaxcmfdmg  &
          /1000,   1,    1,       4,        4,      3/
      data epseig   ,epsl2    ,epsin    ,epserf,  &
           epserfmg, epscmfd2g,epscmfdmg    &
          /1.e-6,1.e-5,1.e-3,1.e-1,   &
           1.e-2,1.e-1,5.e-1/
      data eshift0,eshift  &
          /0.1,0.04/
          !/0.0,0.00/
      data innermg/'bicg'/

! times.h
      data tinit,tth,ttotal,txsec,tppr  &
          /  0.0,  0.0,  0.0, 0.0,  0.0 /

! nodal.h
      data nodalonly /FALSE/
      data is1n /TRUE/
      data epsnodal / 0. /
      data underrelx /100./
      data nsrcexp /4/
      
! geom.h defaults
      data isymang,isymloc,symopt,nrodtype /360,4,'REFL',0/
!
! files.h defaults
      data io1 ,io2 ,io3 ,io4 ,io5 ,io6 ,io7 ,io8 ,io9 ,io10 ,  &          
           io11,io12,io13,io14,io15,io16,io17,io18,io19,io20 ,  &  
           io21,io22,io23,io24,io25,io26,io27,io28,io29,io30    &
           / 1, 2, 3, 4, 5, 6, 7, 8, 9,10, &    
            11,12,13,14,15,16,17,18,19,20, & 
            21,22,23,24,25,26,27,28,29,30/
!
! pinpr.h defaults
      data pinpower /FALSE/

! srchppm.h
      data ppm,targetk /0.0,1.0/
      data srchppm / FALSE /  ! 2013_01_27 . scb
      
! pow.h
      data plevel /1./

! trancntl.inc
      data texpand, texpand2, tswitch, tswitch2  &
          /   1   ,    1    ,   BIG  ,  BIG   /     
      data texpand3, texpand4, tswitch3, tswitch4  &       ! 2012_08_22 . scb
          /   1   ,    1    ,   BIG  ,  BIG   /             ! 2012_08_22 . scb
! 2012_08_23 . scb            
      data bdforder, bdforder0, autodt, errtol, safefac, retry, unit_tm, maxdeltm  &      
          /   5    , 5, FALSE,  1.e-4 ,  0.8  , FALSE , 0.0001 ,    1.   /           
      data maxbound(1),maxbound(2),maxbound(3)   &      
          /   1.25,       1.5,        2.0     /         
      data scrmflag,scram,trip,isvel / FALSE, FALSE, FALSE, FALSE /
      data tripfirst / TRUE /
      data tscrmend / 10000. /
      
      data flagouttr / .false. / ! 2015_06_22 . scb
! added end      
! added in ARTOS ver. 0.2 . 2012_07_05 by SCB           
      data decayht /FALSE/ ! DECAY HEAT
      data decalpha /2.35402E-02,1.89077E-02,1.39236E-02,6.90315E-03,3.56888E-03,3.31633E-03/ ! DECAY HEAT
      data deczeta /1.05345E-01,8.37149E-03,5.20337E-04,4.73479E-05,3.28153E-06,1.17537E-11/  ! DECAY HEAT	
      data epsxsec, epstemp /0. , 0./ ! conditional nodal update
      data epsl2tr, epserftr /1.e-4, 1.e-3/
! added end      
	
! thgeom.inc
      data nperfa,npint,ngt,    rs,    rw,   tw,   rgt  &
          /     1,  264, 25,4.1195,4.7585,0.571,6.1295 /
! thop.inc
      data   pfa0,  powfa0,  tin, cfrperfa, hgap   &
          /21.606,17.67516,286.0,82.12102,10000./
! frzfdbk.inc
      data freezetf,freezedm,writefbv,readdopl,readtmdm   &
          / FALSE, FALSE, FALSE, FALSE, FALSE /     
! thfuel.inc
      data ifeffdopl,wfsurf,wfcl   &
!          /   TRUE,   0.7, 0.3  /          !   original value
          /   TRUE,   0.64, 0.36  /         !   NEA w=0.64 
!         /   TRUE,   0.78, 0.22  /         !   Studsvik w=0.78 
! thcntl.inc
      data fdbk,isflatpower,ntth,epstf / FALSE, TRUE, 0, 1.e-3 /      
      data thetaf,thetac / 0.5, 0.5 /  !parameters of theta-method in integrator used in TH transient

! perturb.inc
      data npbank / 0 /   

! ff.h
      data printff /FALSE/  
 
! added in ARTOS ver. 0.2 . 2012_07_05 by SCB     
! thexpan.inc
      data thexpan, ifreact /FALSE, FALSE/ ! THERMAL EXPANSION
! added end      

      data ilink /-1/  ! added in ARTOS ver. 0.2 . 2012_07_25 by SCB  
      
! 2013_09_16 . scb
      data filemap /'artos.map'/
! added end      

! 2013_09_26 . scb
      data flagxesm / .false. /  ! 2014_07_29 . scb
      data ixeopt, ismopt / 0, 0 /
      data rlambi, rlambxe, rlambpm / 0.28750e-4, 0.209167e-4, 0.355568e-5 /   ! 2014_08_12 . scb . PM data was wrong
      data gammafpi, gammafpxe, gammafppm / 0.06386, 0.00228, 0.0113 /
      call dmalloc(sigxea0,ng2)
      call dmalloc(sigsma0,ng2)
      sigxea0(1) = 1.05279e2*1.0e-24
      sigxea0(2) = 1.45710e6*1.0e-24
      sigsma0(1) = 0*9.07729D-23
      sigsma0(2) = 0*5.52864D-20
! added end

! 2014_04_10 . scb for rst
#ifndef DLL_TASS   ! 2014_04_10 . scb
      data rstrt, rsted / .false. , .false. /
#endif      
      data irstfreq / 0 /
! added end
            !print *, 'epseig : ',epseig   ! 2014_07_24 . scb    
            
      data ifcmfd2g / .true. /   ! 2014_12_05 . scb
	
    end subroutine