! file name = cntl.h

! tracinp == option for trac input use (default=FALSE)
! prtscrn - print on screen
      logical tracinp,prtscrn,transient
      logical qtransient                   ! flag for quasi-transient calculations
      integer nthread
      common /cntli/ tracinp,prtscrn,transient,qtransient
      common /cntli/ nthread

! 2013_05_13 . scb for restart
      logical rstrt
	  logical ssinit, ssfirst, ssdata, trrst

	  common /rstrtbl/ rstrt, ssinit, ssfirst, ssdata, trrst 

      integer itimeadv, irstadv, irstfreq                               !rstrt
      common /irsttr/ itimeadv, irstadv, irstfreq                       !rstrt
	  
      real(8) rsstime                                                    !ssinit
      common /stdy/ rsstime                                              !ssinit

      real(8) rsttime                                                   !rstrt
      common /rstcntl/ rsttime                                          !rstrt
      common /irstcntl/ irstbeg                                         !rstrt

      logical rsted
      common /rstedits/ rsted
! added end

