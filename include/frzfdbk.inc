! file name = frzfdbk.inc - freeze xsec feedback
      logical freezetf,freezedm,writefbv,readdopl,readtmdm
      real :: unifdm,uniftf,uniftm
      common /lfreeze/freezetf,freezedm,writefbv,readdopl,readtmdm
      common /ffreeze/frozentf,frozendm
      common /rfreeze/ unifdm,uniftf,uniftm
