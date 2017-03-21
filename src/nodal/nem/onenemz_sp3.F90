! added in ARTOS ver. 0.2 . 2012_07_24 by SCB        
   subroutine onenemz_sp3(h, xsd, xsd2, xsr, xst, &
                           jinl, jinr, q0, q1, q2, &
                           joutl, joutr, &
                           phibar, phi1st, phi2nd, &
                           phisl, phisr)
	   implicit none
! input      
      real :: xsd, xsd2, xsr, xst, h
      real :: jinl(2), jinr(2)
      real :: q0(2), q1(2), q2(2)
! output
      real :: joutl(2), joutr(2)
      real :: phibar(2), phi1st(2), phi2nd(2)
! local
      real :: invh, beta0, beta2, t(16)
      real :: ca0(4), ca2(4), b0(4)
      real :: cs0(4), cs2(4), b2(4)
      real :: cf(4), b1(4)
      
      real :: rcs0(4), rdet, rcf(4)
      real :: ab(4), abc(4), x(2), y(2), z(2)
      real :: q0t(2), q1t(2), q2t(2)
      real :: tlhs(4), trhs(2), rtlhs(4)

      real :: tjn(4), tjf(4), t0(4), tf(4)

      real :: rca2(4)

      real :: phisl(2), phisr(2), jp(2)


#ifdef DEBUG
      xst = 0.3
      h = 10.
      xsr = 0.5

      q0(1) = 0.5
      q1(1) = 1.5
      q2(1) = 2.5

      q0(2) = -2./3.*q0(1)
      q1(2) = -2./3.*q1(1)
      q2(2) = -2./3.*q2(1)

      jinl(1) = 2.5
      jinl(2) = 1.0

      jinr(1) = 3.1
      jinr(2) = 1.5

      xsd = 1./3./xst
      xsd2 = 3./7./xst
#endif


      invh = 1.d0/h
      beta0 = xsd*invh
      beta2 = xsd2*invh
      call initnem_param(beta0, beta2, t)
      call initnem_mat(invh, xsr, xst, xsd, xsd2, t, &
                        ca0, ca2, b0, cs0, cs2, b2, cf, b1)

! first moment balance equation
      x(1) = -jinl(1) +jinr(1)
      x(2) = -jinl(2) +jinr(2)
      call rhs(b1,x,q1,q1t)
      call rmat(cf,rcf)

      call m2xn1(rcf,q1t,phi1st)

! nodal balance equation
! second moment balance equation
      call rmat(cs0,rcs0) 
      call m2xn2(ca0,rcs0,ab)
      call m2xn2(ab,cs2,abc) 

      x(1) = jinl(1) +jinr(1)
      x(2) = jinl(2) +jinr(2)
      call rhs(b0,x,q0,q0t) ! nodal balance 
      call rhs(b2,x,q2,q2t) ! second moment balance

      call m2xn1(ab,q2t,y)

      tlhs = ca2 -abc
      trhs = q0t -y
      call rmat(tlhs,rtlhs) 
      call m2xn1(rtlhs,trhs,phi2nd)
!

      call m2xn1(cs2,phi2nd,y)
      trhs = q2t -y
      call m2xn1(rcs0,trhs,phibar)

! jout
      tjn(1) = t(1)
      tjn(2) = t(3)
      tjn(3) = t(5)
      tjn(4) = t(7)

      tjf(1) = t(2)
      tjf(2) = t(4)
      tjf(3) = t(6)
      tjf(4) = t(8)

      t0(1) = t(9)
      t0(2) = t(11)
      t0(3) = t(10)
      t0(4) = t(12)

      tf(1) = t(13)
      tf(2) = t(15)
      tf(3) = t(14)
      tf(4) = t(16)

      joutr(1) = tjn(1)*jinr(1) +tjn(2)*jinr(2) &
                  +tjf(1)*jinl(1) +tjf(2)*jinl(2) &
                  +t0(1)*phibar(1) +t0(2)*phibar(2) &
                  +tf(1)*phi1st(1) +tf(2)*phi1st(2) &
                  -3.5d0*(t0(1)*phi2nd(1) +t0(2)*phi2nd(2)) 
      joutr(2) = tjn(3)*jinr(1) +tjn(4)*jinr(2) &
                  +tjf(3)*jinl(1) +tjf(4)*jinl(2) &
                  +t0(3)*phibar(1) +t0(4)*phibar(2) &
                  +tf(3)*phi1st(1) +tf(4)*phi1st(2) &
                  -3.5d0*(t0(3)*phi2nd(1) +t0(4)*phi2nd(2)) 

      joutl(1) = tjn(1)*jinl(1) +tjn(2)*jinl(2) &
                  +tjf(1)*jinr(1) +tjf(2)*jinr(2) &
                  +t0(1)*phibar(1) +t0(2)*phibar(2) &
                  -tf(1)*phi1st(1) -tf(2)*phi1st(2) &
                  -3.5d0*(t0(1)*phi2nd(1) +t0(2)*phi2nd(2)) 
      joutl(2) = tjn(3)*jinl(1) +tjn(4)*jinl(2) &
                  +tjf(3)*jinr(1) +tjf(4)*jinr(2) &
                  +t0(3)*phibar(1) +t0(4)*phibar(2) &
                  -tf(3)*phi1st(1) -tf(4)*phi1st(2) &
                  -3.5d0*(t0(3)*phi2nd(1) +t0(4)*phi2nd(2)) 


      jp = joutl + jinl
      phisl(1) = 8.d0/25.d0*(7.d0*jp(1)+3.d0*jp(2))
      phisl(2) = 8.d0/25.d0*(1.d0*jp(1)+4.d0*jp(2))

      jp = joutr + jinr
      phisr(1) = 8.d0/25.d0*(7.d0*jp(1)+3.d0*jp(2))
      phisr(2) = 8.d0/25.d0*(1.d0*jp(1)+4.d0*jp(2))
      !
      !  pause
      !
      !print *, h
      !print *, xsd
      !print *, xsd2
      !print *, xsr
      !print *, xst
      !print *, jinl
      !print *, jinr
      !print *, q0
      !print *, q1
      !print *, q2
      !print *, joutl
      !print *, joutr
      !print *, phibar
      !print *, phi1st
      !print *, phi2nd
      !print *, phisl
      !print *, phisr
      !
      !  pause
        
      return
      contains
         subroutine rmat(a,ra)
            implicit none
            real :: a(4), ra(4)
            real :: rdet

            rdet=1.d0/(a(1)*a(4)-a(2)*a(3))
            ra(1)=rdet*a(4)
            ra(2)=-rdet*a(2)
            ra(3)=-rdet*a(3)
            ra(4)=rdet*a(1)
            return
         end subroutine

         subroutine m2xn2(a,b,c)   ! a*b=c
            implicit none
            real :: a(4), b(4), c(4)

            c(1)=a(1)*b(1)+a(2)*b(3)
            c(2)=a(1)*b(2)+a(2)*b(4)
            c(3)=a(3)*b(1)+a(4)*b(3)
            c(4)=a(3)*b(2)+a(4)*b(4)
         
            return
         end subroutine

         subroutine m2xn1(a,b,c)   ! a*b=c
            implicit none
            real :: a(4), b(2), c(2)

            c(1)=a(1)*b(1)+a(2)*b(2)
            c(2)=a(3)*b(1)+a(4)*b(2)
        
            return
         end subroutine

         subroutine rhs(b,y,qi,qo)   ! qi + b[2x2]*y[2x1] = qo
            implicit none
            real :: b(4), y(2), qi(2), qo(2)
            qo(1) = qi(1) +b(1)*y(1) +b(2)*y(2)
            qo(2) = qi(2) +b(3)*y(1) +b(4)*y(2)
            return
         end subroutine
                                                                 
   end subroutine


