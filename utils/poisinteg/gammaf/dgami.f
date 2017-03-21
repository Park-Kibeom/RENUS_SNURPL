      double precision function dgami (a, x)
c july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
c
c evaluate the incomplete gamma function defined by
c
c gami = integral from t = 0 to x of exp(-t) * t**(a-1.0) .
c
c gami is evaluated for positive values of a and non-negative values
c of x.  a slight deterioration of 2 or 3 digits accuracy will occur
c when gami is very large or very small, because logarithmic variables
c are used.
c
      double precision a, x, factor, dlngam, dgamit, dexp, dlog
      external dexp, dgamit, dlngam, dlog
c
      if (a.le.0.d0) call seteru (25hdgami   a must be gt zero, 25, 1,2)
      if (x.lt.0.d0) call seteru (25hdgami   x must be ge zero, 25, 2,2)
c
      dgami = 0.d0
      if (x.eq.0.0d0) return
c
c the only error possible in the expression below is a fatal overflow.
      factor = dexp (dlngam(a) + a*dlog(x))
c
      dgami = factor * dgamit (a, x)
c
      return
      end
