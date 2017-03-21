      double precision function dgamma (x)
c jan 1984 edition.  w. fullerton, c3, los alamos scientific lab.
c jan 1994 wpp@ips.id.ethz.ch, ehg@research.att.com   declare xsml
      double precision x, gamcs(42), dxrel, pi, sinpiy, sq2pil, xmax,
     1  xmin, y, d9lgmc, dcsevl, d1mach, dexp, dint, dlog,
     2  dsin, dsqrt, xsml
      external d1mach, d9lgmc, dcsevl, dexp, dint, dlog, dsin, dsqrt,
     1  initds
c
c series for gam        on the interval  0.          to  1.00000e+00
c                                        with weighted error   5.79e-32
c                                         log weighted error  31.24
c                               significant figures required  30.00
c                                    decimal places required  32.05
c
      data gam cs(  1) / +.8571195590 9893314219 2006239994 2 d-2      /
      data gam cs(  2) / +.4415381324 8410067571 9131577165 2 d-2      /
      data gam cs(  3) / +.5685043681 5993633786 3266458878 9 d-1      /
      data gam cs(  4) / -.4219835396 4185605010 1250018662 4 d-2      /
      data gam cs(  5) / +.1326808181 2124602205 8400679635 2 d-2      /
      data gam cs(  6) / -.1893024529 7988804325 2394702388 6 d-3      /
      data gam cs(  7) / +.3606925327 4412452565 7808221722 5 d-4      /
      data gam cs(  8) / -.6056761904 4608642184 8554829036 5 d-5      /
      data gam cs(  9) / +.1055829546 3022833447 3182350909 3 d-5      /
      data gam cs( 10) / -.1811967365 5423840482 9185589116 6 d-6      /
      data gam cs( 11) / +.3117724964 7153222777 9025459316 9 d-7      /
      data gam cs( 12) / -.5354219639 0196871408 7408102434 7 d-8      /
      data gam cs( 13) / +.9193275519 8595889468 8778682594 0 d-9      /
      data gam cs( 14) / -.1577941280 2883397617 6742327395 3 d-9      /
      data gam cs( 15) / +.2707980622 9349545432 6654043308 9 d-10     /
      data gam cs( 16) / -.4646818653 8257301440 8166105893 3 d-11     /
      data gam cs( 17) / +.7973350192 0074196564 6076717535 9 d-12     /
      data gam cs( 18) / -.1368078209 8309160257 9949917230 9 d-12     /
      data gam cs( 19) / +.2347319486 5638006572 3347177168 8 d-13     /
      data gam cs( 20) / -.4027432614 9490669327 6657053469 9 d-14     /
      data gam cs( 21) / +.6910051747 3721009121 3833697525 7 d-15     /
      data gam cs( 22) / -.1185584500 2219929070 5238712619 2 d-15     /
      data gam cs( 23) / +.2034148542 4963739552 0102605193 2 d-16     /
      data gam cs( 24) / -.3490054341 7174058492 7401294910 8 d-17     /
      data gam cs( 25) / +.5987993856 4853055671 3505106602 6 d-18     /
      data gam cs( 26) / -.1027378057 8722280744 9006977843 1 d-18     /
      data gam cs( 27) / +.1762702816 0605298249 4275966074 8 d-19     /
      data gam cs( 28) / -.3024320653 7353062609 5877211204 2 d-20     /
      data gam cs( 29) / +.5188914660 2183978397 1783355050 6 d-21     /
      data gam cs( 30) / -.8902770842 4565766924 4925160106 6 d-22     /
      data gam cs( 31) / +.1527474068 4933426022 7459689130 6 d-22     /
      data gam cs( 32) / -.2620731256 1873629002 5732833279 9 d-23     /
      data gam cs( 33) / +.4496464047 8305386703 3104657066 6 d-24     /
      data gam cs( 34) / -.7714712731 3368779117 0390152533 3 d-25     /
      data gam cs( 35) / +.1323635453 1260440364 8657271466 6 d-25     /
      data gam cs( 36) / -.2270999412 9429288167 0231381333 3 d-26     /
      data gam cs( 37) / +.3896418998 0039914493 2081663999 9 d-27     /
      data gam cs( 38) / -.6685198115 1259533277 9212799999 9 d-28     /
      data gam cs( 39) / +.1146998663 1400243843 4761386666 6 d-28     /
      data gam cs( 40) / -.1967938586 3451346772 9510399999 9 d-29     /
      data gam cs( 41) / +.3376448816 5853380903 3489066666 6 d-30     /
      data gam cs( 42) / -.5793070335 7821357846 2549333333 3 d-31     /
c
      data pi / 3.1415926535 8979323846 2643383279 50 d0 /
c sq2pil is 0.5*alog(2*pi) = alog(sqrt(2*pi))
      data sq2pil / 0.9189385332 0467274178 0329736405 62 d0 /
      data ngam, xmin, xmax, xsml, dxrel / 0, 4*0.d0 /
c
      if (ngam.ne.0) go to 10
      ngam = initds (gamcs, 42, 0.1*sngl(d1mach(3)) )
c
      call d9gaml (xmin, xmax)
      xsml = dexp (dmax1 (dlog(d1mach(1)), -dlog(d1mach(2)))+0.01d0)
      dxrel = dsqrt (d1mach(4))
c
 10   y = dabs(x)
      if (y.gt.10.d0) go to 50
c
c compute gamma(x) for -xbnd .le. x .le. xbnd.  reduce interval and find
c gamma(1+y) for 0.0 .le. y .lt. 1.0 first of all.
c
      n = x
      if (x.lt.0.d0) n = n - 1
      y = x - dble(float(n))
      n = n - 1
      dgamma = 0.9375d0 + dcsevl (2.d0*y-1.d0, gamcs, ngam)
      if (n.eq.0) return
c
      if (n.gt.0) go to 30
c
c compute gamma(x) for x .lt. 1.0
c
      n = -n
      if (x.eq.0.d0) call seteru (14hdgamma  x is 0, 14, 4, 2)
      if (x.lt.0.0d0 .and. x+dble(float(n-2)).eq.0.d0) call seteru (
     1  31hdgamma  x is a negative integer, 31, 4, 2)
      if (x.lt.(-0.5d0) .and. dabs((x-dint(x-0.5d0))/x).lt.dxrel) call
     1  seteru (68hdgamma  answer lt half precision because x too near n
     2egative integer, 68, 1, 1)
      if (y.lt.xsml) call seteru (
     1  54hdgamma  x is so close to 0.0 that the result overflows,
     2  54, 5, 2)
c
      do 20 i=1,n
        dgamma = dgamma/(x+dble(float(i-1)) )
 20   continue
      return
c
c gamma(x) for x .ge. 2.0 and x .le. 10.0
c
 30   do 40 i=1,n
        dgamma = (y+dble(float(i))) * dgamma
 40   continue
      return
c
c gamma(x) for dabs(x) .gt. 10.0.  recall y = dabs(x).
c
 50   if (x.gt.xmax) call seteru (32hdgamma  x so big gamma overflows,
     1  32, 3, 2)
c
      dgamma = 0.d0
      if (x.lt.xmin) call seteru (35hdgamma  x so small gamma underflows
     1  , 35, 2, 0)
      if (x.lt.xmin) return
c
      dgamma = dexp ((y-0.5d0)*dlog(y) - y + sq2pil + d9lgmc(y) )
      if (x.gt.0.d0) return
c
      if (dabs((x-dint(x-0.5d0))/x).lt.dxrel) call seteru (
     1  61hdgamma  answer lt half precision, x too near negative integer
     2  , 61, 1, 1)
c
      sinpiy = dsin (pi*y)
      if (sinpiy.eq.0.d0) call seteru (
     1  31hdgamma  x is a negative integer, 31, 4, 2)
c
      dgamma = -pi/(y*sinpiy*dgamma)
c
      return
      end
