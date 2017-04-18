import numpy as np
import c2raytools as c2t
import cosmological_constants as cosm_cons

def matter_variance(z, M=None, R=None, power_spectra=None):
	"""
	Equation 19
	M: It should have unit of solar mass.
	R: It should have a unit of Mpc.
	"""
	assert M is not None or R is not None
	rho_mean  = cosm_cons.Omega_m*critical_density(z) #kg/m^3
	rho_mean_ = rho_mean*1.47e37
	if R is None: R = (3*M/(4*np.pi*rho_mean_))**0.3333
	if power_spectra is None or len(power_spectra)!=2:
		Pk,k = _power_spectrum(z, [-4,0.5])
	else:
		Pk,k = power_spectra
	integ = k**2*Pk*(3*_j1(k*R)/k/R)**2/2/np.pi**2
	var_M = np.trapz(integ, k)
	return var_M*cosm_cons.sigma_8**2/1.3129049664505633e-06

def critical_density(z):
	H = c2t.hubble_parameter(z)*1000/mpc_to_m(1)
	return 3*H**2/(8*np.pi*cosm_cons.G_grav)

def _j1(x):
	# Used in eq 19
	return (np.sin(x)-x*np.cos(x))/x**2

def _get_overdensity(density):
	# Equation 11
	return density/density.mean(dtype=np.float64) - 1.

def _power_spectrum(z, ks):
	# Used in eq 19
	if type(ks)==np.float: k = ks
	elif len(ks)==2: k = 10**np.linspace(ks[0], ks[1], 1000)
	else: k = ks
	#RH = 3e5/c2t.hubble_parameter(z)  #in Mpc
	#kH = 2*np.pi/RH
	#ns = np.ones(k.size)
	#ns[k>kH] = -3
	Pk_prim  = k#**ns
	Tk = transfer_function(k)
	Dz = growth_factor(z)
	Pk = Pk_prim*Tk*Dz
	return Pk,k

def mpc_to_m(x):
	return 3.086e22*x

def growth_factor(z):
	"""
	Equation 16
	"""
	a = 1./(1.+z)
	Omg_l, Omg_k, Omg_m = cosm_cons.Omega_lam, 0, cosm_cons.Omega_m
	aa = np.linspace(0,a,1000)
	integ = aa**1.5/(Omg_l*aa**3+Omg_k*aa+Omg_m)**1.5
	Dz = np.trapz(integ, aa) * (Omg_l*a**3+Omg_k*a+Omg_m)**0.5/a**1.5
	return Dz

def critical_overdensity_collapse(z):
	"""
	Equation 21
	"""
	return 1.686/growth_factor(z)

def transfer_function(k):
	q  = k/(cosm_cons.Omega_m*cosm_cons.h**2)    ## Mpc^-1
	Lq = np.log(np.e+1.84*q)
	Cq = 14.4 + 325/(1+60.5*q**1.11)
	Tk = Lq/(Lq + Cq*q**2)
	return Tk

def Tvir_to_mass(Tvir, z):
	"""
	Equation 26
	Tvir: in Kelvin.
	Return
	------
	mass: in solar mass.
	"""
	a = 1./(1+z)
	Omega_m_z = cosm_cons.Omega_m*a**3/(cosm_cons.Omega_m*a**3+cosm_cons.Omega_lam)
	d = Omega_m_z-1
	Del_c = 18*np.pi**2+82*d-39*d**2
	mu = 0.6                        #fully ionized case considered
	bla = 1.98e4*(mu/0.6)*(cosm_cons.Omega_m*Del_c/Omega_m_z/(18*np.pi**2))**0.33*(1+z)/10
	mass = (Tvir/bla)**1.5*10**8/cosm_cons.h
	return mass




