from numpy import array
from sympy import solveset, symbols, linsolve

lam_xe, lam_i = symbols('lam_xe lam_i')
rho_xe, rho_i = symbols('rho_xe rho_i')
yld_xe, yld_i = symbols('yld_xe yld_i')
sig_ff,sig_ft = symbols('sig_ff sig_ft')
sigf_xe, sigt_xe = symbols('sigf_xe sigt_xe')
phi_f, phi_t = symbols('phi_f phi_t')

fxe = - lam_xe * rho_xe - (sigf_xe * phi_f + sigt_xe * phi_t) * rho_xe  + yld_xe * (sig_ff * phi_f + sig_ft * phi_t) + lam_i * rho_i
fi  = - lam_i  * rho_i                                                  + yld_i  * (sig_ff * phi_f + sig_ft * phi_t)

slv = linsolve([fxe,fi],(rho_i, rho_xe)) 

fnorm = 0.11685
mysub={
	lam_i  : 0.287500E-4,
	lam_xe : 0.209167E-4,
	yld_i  : 0.06386,
	yld_xe : 0.00228,
	sigt_xe: 1.45710e6*1.0e-24,
	sigf_xe: 1.05279e2*1.0e-24,
	sig_ff : 0.190224E-02,
	sig_ft : 0.343581E-01,
#	phi_f  : 3.83997188E+14,
#	phi_t  : 1.02478775E+14,
	phi_f  : fnorm * 59.982,
	phi_t  : fnorm * 16.657,
}

print mysub

print slv
for sol in slv:
	aio135 = sol[0].evalf(subs=mysub)
	axe135 = sol[1].evalf(subs=mysub)

# form SKETCH
sxe135 = 1.1061/0.66901E+09 * 1.E+24
sio135 = 1.4626/0.15473E+09 * 1.E+24

print "Analytic: Xe: %6.6E, I: %6.6E"%(axe135, aio135)
print "SKETCH  : Xe: %6.6E, I: %6.6E"%(sxe135, sio135)

print (sigf_xe * phi_f + sigt_xe * phi_t).evalf(subs=mysub)
print (sig_ff * phi_f + sig_ft * phi_t).evalf(subs=mysub)