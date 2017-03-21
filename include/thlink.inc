! added in ARTOS ver. 0.2 ( for System-Neutronics Coupling ). 2012_07_03 by SCB
! thlink.h - T/H link control data

	  character*72 :: filemap
	  common /marsmap/ filemap 

	  character*12 :: thcode
      parameter (icouple=1, thcode='MARS')

!      logical ifcobra
      common /thcntlx/ilink,iconv,nssstep,ntrstep,nmars,ifcobra
      common /thgeodim/nchanmars,nlevel,nthn,nbank,nchancob
      common /fthlink/dtsys,qpavg0,qpavg,rqpavg,tfws,tfwc,qtripavg,sstime,trtime
      common /ccobra/fnormhin,fnormgin,fnormpexit !normalization fac
	  common /pow0/ power


