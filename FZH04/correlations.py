import numpy as np
import barkana_loeb_2001
import cosmological_constants as cosm_cons
from scipy.special import erfinv
from scipy.misc import derivative
from scipy.integrate import simps
from functions import *

def correlation_xx(r, z, zeta=40):
	#m       = radius_to_mass(r,z)
	m_min   = minimum_source_mass(z, Tvir=1e4)
	#if m_min>m:
	#	print "This radius is less than the minimum source radius."
	#	return 0
	#mm = 10**np.linspace(np.log10(m_min), 20, 100)
	mm = 10**np.linspace(np.log10(m_min), 20, 100)
	Vm = mm/mean_density(z)
	V0 = Vm0(mm, z, r12=r)
	dndm    = mass_function(mm, z, zeta=zeta)/mm
	Q0 = np.trapz(dndm*V0, mm)
	Q1 = np.trapz(dndm*(Vm-V0), mm)
	xi_xx = (1-np.exp(-Q0)) + np.exp(-Q0)*(1-np.exp(-Q1))**2
	return xi_xx

def correlation_xd(r, z, zeta=40):
	return 0

def correlation_dd(r, z):
	assert type(r) == np.ndarray and r.size!=0
	if r.size == 2: rs = 10**np.linspace(np.log10(r[0]), np.log10(r[1]), 100)
	else: rs = r
	ks    = 2*np.pi/rs
	Pk,k  =  barkana_loeb_2001._power_spectrum(z, ks)
	xi_dd = np.fft.ifftshift(Pk)
	return np.real(xi_dd)

	
"""
def correlation_xx_clustering(r, z, zeta=20):
	B0 = barrier_density(z, zeta=zeta, var_m=0)
	m  = radius_to_mass(r)
	if type(m)==np.float:
		Bm = linear_fit_barrier(m, z, zeta=zeta, var_m=None)
		sig2 = barkana_loeb_2001.matter_variance(z, M=m, power_spectra=None) 
	else:
		sig2 = np.array([barkana_loeb_2001.matter_variance(z, M=m1, power_spectra=None) for m1 in m])
		Bm   = np.array([linear_fit_barrier(m1, z, zeta=zeta, var_m=None) for m1 in m])
	bm = 1 + B0**2/sig2/Bm
	def get_b_(m):
		mm = 10**np.linspace(np.log10(m), 20, 100)
		Vm = mm/mean_density(z)
"""
		

