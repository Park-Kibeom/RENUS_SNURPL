! added in ARTOS ver. 0.2 . 2012_07_24 by SCB  
      subroutine initnem_mat(invh, xsr, xst, xsd, xsd2, t, &
                             ca0, ca2, b0, &
                             cs0, cs2, b2, &
                             cf, b1)
      implicit none
! input
      real :: invh, t(16), xsr, xst, xsd, xsd2
! output
      real :: ca0(4), ca2(4), b0(4)
      real :: cs0(4), cs2(4), b2(4)
      real :: cf(4), b1(4)

      real :: hxsd, hxsd2

      hxsd = xsd*invh*invh
      hxsd2 = xsd2*invh*invh

! nodal balance:         ca0 * phi + ca2 * phi_s = b0 * jin + q 
! second moment balance: cs0 * phi + cs2 * phi_s = b2 * jin + q_s 
! first moment balance:  cf * phi_f = b1 * jin + q_f

! nodal balance
      ca0(1) = 2.d0*t(9)*invh +xsr
      ca0(2) = 2.d0*t(11)*invh -2.d0*xsr
      ca0(3) = 2.d0*t(10)*invh -2.d0/3.d0*xsr
      ca0(4) = 2.d0*t(12)*invh +4.d0/3.d0*xsr +5.d0/3.d0*xst

      ca2(1) = -7.d0*t(9)*invh
      ca2(2) = -7.d0*t(11)*invh
      ca2(3) = -7.d0*t(10)*invh
      ca2(4) = -7.d0*t(12)*invh

      b0(1) = -(-1.d0 +t(1) +t(2))*invh
      b0(2) = -(t(3) +t(4))*invh
      b0(3) = -(t(5) +t(6))*invh
      b0(4) = -(-1.d0 +t(7) +t(8))*invh

! first moment balance
      cf(1) = (60.d0 -224.d0/5.d0*t(13) -96.d0/5.d0*t(14))*hxsd +xsr
      cf(2) = (-224.d0/5.d0*t(15) -96.d0/5.d0*t(16))*hxsd -2.d0*xsr
      cf(3) = (-32.d0/5.d0*t(13) -128.d0/5.d0*t(14))*hxsd2 -2.d0/3.d0*xsr
      cf(4) = (60.d0 -32.d0/5.d0*t(15) -128.d0/5.d0*t(16))*hxsd2 +4.d0/3.d0*xsr +5.d0/3.d0*xst
      
      b1(1) = 7.d0 +7.d0*t(1) -7.d0*t(2) +3.d0*t(5) -3.d0*t(6)
      b1(2) = 3.d0 +7.d0*t(3) -7.d0*t(4) +3.d0*t(7) -3.d0*t(8)
      b1(3) = 1.d0 +t(1) -t(2) +4.d0*t(5) -4.d0*t(6)
      b1(4) = 4.d0 +t(3) -t(4) +4.d0*t(7) -4.d0*t(8)
      
      b1(1) = 16.d0/5.d0*hxsd*b1(1)    
      b1(2) = 16.d0/5.d0*hxsd*b1(2)   
      b1(3) = 16.d0/5.d0*hxsd2*b1(3)   
      b1(4) = 16.d0/5.d0*hxsd2*b1(4) 
      
! second moment balance
      cs0(1) = -28.d0*hxsd +1568.d0/25.d0*t(9)*hxsd +672.d0/25.d0*t(10)*hxsd
      cs0(2) = 1568.d0/25.d0*t(11)*hxsd +672.d0/25.d0*t(12)*hxsd          
      cs0(3) = 224.d0/25.d0*t(9)*hxsd2 +896.d0/25.d0*t(10)*hxsd2
      cs0(4) = -28.d0*hxsd2 +224.d0/25.d0*t(11)*hxsd2 +896.d0/25.d0*t(12)*hxsd2

      cs2(1) = 140.d0*hxsd -5488.d0/25.d0*t(9)*hxsd -2352.d0/25.d0*t(10)*hxsd +xsr
      cs2(2) = -5488.d0/25.d0*t(11)*hxsd -2352.d0/25.d0*t(12)*hxsd -2.d0*xsr
      cs2(3) = -784.d0/25.d0*t(9)*hxsd2 -3136.d0/25.d0*t(10)*hxsd2 -2.d0/3.d0*xsr
      cs2(4) = 140.d0*hxsd2 -784.d0/25.d0*t(11)*hxsd2 -3136.d0/25.d0*t(12)*hxsd2 &
               +4.d0/3.d0*xsr +5.d0/3.d0*xst

      b2(1) = 7.d0 +7.d0*t(1) +7.d0*t(2) +3.d0*t(5) +3.d0*t(6)
      b2(2) = 3.d0 +7.d0*t(3) +7.d0*t(4) +3.d0*t(7) +3.d0*t(8)
      b2(3) = 1.d0 +t(1) +t(2) +4.d0*t(5) +4.d0*t(6)
      b2(4) = 4.d0 +t(3) +t(4) +4.d0*t(7) +4.d0*t(8)

      b2(1) = -112.d0/25.d0*hxsd*b2(1)    
      b2(2) = -112.d0/25.d0*hxsd*b2(2)   
      b2(3) = -112.d0/25.d0*hxsd2*b2(3)   
      b2(4) = -112.d0/25.d0*hxsd2*b2(4) 

      return
      end subroutine
