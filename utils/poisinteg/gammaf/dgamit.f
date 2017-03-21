      double precision function dgamit (a, x)
c july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
c
c evaluate tricomi-s incomplete gamma function defined by
c
c gamit = x**(-a)/gamma(a) * integral t = 0 to x of exp(-t) * t**(a-1.)
c
c and analytic continuation for a .le. 0.0.  gamma(x) is the complete
c gamma function of x.  gamit is evaluated for arbitrary real values of
c a and for non-negative values of x (even though gamit is defined for
c x .lt. 0.0).
c
c      a slight deterioration of 2 or 3 digits accuracy will occur when
c gamit is very large or very small in absolute value, because log-
c arithmic variables are used.  also, if the parameter a is very close
c to a negative integer (but not a negative integer), there is a loss
c of accuracy, which is reported if the result is less than half
c machine precision.
c
c ref. -- w. gautschi, an evaluation procedure for incomplete gamma
c functions, acm trans. math. software.
c
      double precision a, x, aeps, ainta, algap1, alneps, alng, alx,
     1  bot, h, sga, sgngam, sqeps, t, d1mach, dgamr, d9gmit, d9lgit,
     2  dlngam, d9lgic, dint, dexp, dlog, dsqrt
      external d1mach, d9gmit, d9lgic, d9lgit, dexp, dgamr, dint,
     1  dlngam, dlog, dsqrt
c
      data alneps, sqeps, bot / 3*0.d0 /
c
      if (alneps.ne.0.d0) go to 10
      alneps = -dlog (d1mach(3))
      sqeps = dsqrt (d1mach(4))
      bot = dlog (d1mach(1))
c
 10   if (x.lt.0.d0) call seteru (21hdgamit  x is negative, 21, 2, 2)
c
      if (x.ne.0.d0) alx = dlog (x)
      sga = 1.0d0
      if (a.ne.0.d0) sga = dsign (1.0d0, a)
      ainta = dint (a + 0.5d0*sga)
      aeps = a - ainta
c
      if (x.gt.0.d0) go to 20
      dgamit = 0.0d0
      if (ainta.gt.0.d0 .or. aeps.ne.0.d0) dgamit = dgamr(a+1.0d0)
      return
c
 20   if (x.gt.1.d0) go to 30
      if (a.ge.(-0.5d0) .or. aeps.ne.0.d0) call dlgams (a+1.0d0, algap1,
     1  sgngam)
      dgamit = d9gmit (a, x, algap1, sgngam, alx)
      return
c
 30   if (a.lt.x) go to 40
      t = d9lgit (a, x, dlngam(a+1.0d0))
      if (t.lt.bot) call erroff
      dgamit = dexp (t)
      return
c
 40   alng = d9lgic (a, x, alx)
c
c evaluate dgamit in terms of dlog (dgamic (a, x))
c
      h = 1.0d0
      if (aeps.eq.0.d0 .and. ainta.le.0.d0) go to 50
c
      call dlgams (a+1.0d0, algap1, sgngam)
      t = dlog (dabs(a)) + alng - algap1
      if (t.gt.alneps) go to 60
c
      if (t.gt.(-alneps)) h = 1.0d0 - sga * sgngam * dexp(t)
      if (dabs(h).gt.sqeps) go to 50
c
      call erroff
      call seteru (32hdgamit  result lt half precision, 32, 1, 1)
c
 50   t = -a*alx + dlog(dabs(h))
      if (t.lt.bot) call erroff
      dgamit = dsign (dexp(t), h)
      return
c
 60   t = t - a*alx
      if (t.lt.bot) call erroff
      dgamit = -sga * sgngam * dexp(t)
      return
c
      end
