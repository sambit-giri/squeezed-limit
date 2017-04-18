import numpy as np
import barkana_loeb_2001
from scipy.special import sici
import functions

def b_mz(m, z, zeta=40):
	Bm = functions.linear_fit_barrier(m, z, zeta=zeta)
	B0 = functions.barrier_density(z, zeta=zeta, var_m=0)
	var_m = barkana_loeb_2001.matter_variance(z, M=m, power_spectra=None)
	return 1 + B0**2/var_m/Bm
	
def profile_function(k, m=0):
	u_k = np.exp(-k**2/2)      # For Gaussian density profile 
	#### NFW profile
	#Si, Ci = sici(k)
	#rho_mean_ = mean_density(z)
	#u_k = 	
	return u_k

def P_1h(k, z, Tvir=1e4, zeta=40):
	m_min = barkana_loeb_2001.Tvir_to_mass(Tvir, z)
	mm    = 10**np.linspace(np.log10(m_min), 20, 100)
	rho_mean_ = functions.mean_density(z)
	u_ks = np.array([profile_function(k, m=m) for m in mm])
	mm_dndm = functions.mass_function(mm, z, zeta=zeta)
	dndm    = mm_dndm/mm
	p1h     = np.trapz(dndm*(mm/rho_mean_)**2*u_ks**2, mm)
	return p1h
	
