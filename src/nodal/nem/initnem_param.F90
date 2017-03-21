! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
      subroutine initnem_param(beta0, beta2, t)
      implicit none
      
      real :: beta0, beta2, t(16)
      real :: rxj, rx0, rxf
      real :: f(14)

      rx0 = 1.d0/(5.d0+128.d0*beta2+32.d0*beta0*(7.d0+160.d0*beta2))
      rxf = 1.d0/(25.d0+384.d0*beta2+96.d0*beta0*(7.d0+96.d0*beta2))
      rxj = rx0*rxf

      f(1) = 7.d0+96.d0*beta2
      f(2) = 7.d0+160.d0*beta2
      f(3) = 5.d0+128.d0*beta2
      f(4) = 25.d0+384.d0*beta2
      f(5) = 7.d0+128.d0*beta2
      f(6) = 1.d0+24.d0*beta0
      f(7) = 1.d0+40.d0*beta0
      f(8) = 5.d0+48.d0*beta2
      f(9) = -5.d0+3072.d0*beta0*beta2
      f(10) = 5.d0+84.d0*beta0
      f(11) = 5.d0+128.d0*beta2
      f(12) = 5.d0+224.d0*beta0
      f(13) = 25.d0+384.d0*beta2
      f(14) = 25.d0+672.d0*beta0

      t(1) = -rxj*(-f(3)*f(4)+3072.d0*f(1)*f(2)*beta0*beta0)
      t(2) = -64.d0*rxj*beta0*(35.d0+256.d0*f(8)*beta2)
      t(3) = -768.d0*rxj*beta0*(5.d0+24.d0*f(5)*beta0+96.d0*beta2)
      t(4) = 192.d0*f(9)*rxj*beta0
      t(5) = -256.d0*rxj*beta2*(5.d0+24.d0*f(5)*beta0+96.d0*beta2)
      t(6) = 64.d0*f(9)*rxj*beta2

      t(7) = rxj*(125.d0+8960.d0*beta0+150528.d0*beta0*beta0-49152.d0*f(6)*f(7)*beta2*beta2)
      t(8) = -256.d0*rxj*(5.d0+64.d0*f(10)*beta0)*beta2
      t(9) = 20.d0*f(3)*rx0*beta0
      t(10) = -640.d0*rx0*beta0*beta2
      t(11) = -1920.d0*rx0*beta0*beta2
      t(12) = 20.d0*f(12)*rx0*beta2

      t(13) = 30.d0*f(4)*rxf*beta0
      t(14) = -2880.d0*rxf*beta0*beta2
      t(15) = -8640.d0*rxf*beta0*beta2
      t(16) = 30.d0*f(14)*rxf*beta2

      return
      end subroutine

