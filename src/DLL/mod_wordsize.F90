! added in ARTOS ver. 0.2 ( for System-Neutronics Coupling ). 2012_07_03 by SCB
! wordsize - module for determining word sizes in bytes
	module wordsize
	   parameter (N_SIG_DIGITS=6,N_INT_ORDER=8)
	   parameter (NBF=selected_real_kind(N_SIG_DIGITS))
	   parameter (NBI=selected_int_kind(N_INT_ORDER))
	end module