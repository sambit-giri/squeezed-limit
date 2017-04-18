import numpy as np
import barkana_loeb_2001
import cosmological_constants as cosm_cons
from scipy.special import erfinv
from scipy.misc import derivative
from scipy.integrate import simps

def ionized_mass(mass_galaxy, zeta=40):
	"""
	Equation 1
	"""
	return mass_galaxy*zeta

def barrier_density(z, zeta=40, var_m=None, m=None, m_min=None):
	"""
	Equation 4
	"""
	delta_c = barkana_loeb_2001.critical_overdensity_collapse(z)
	K_zeta  = erfinv(1.-1./zeta)
	if var_m is None: var_m   = barkana_loeb_2001.matter_variance(z, M=m, power_spectra=None)
	if m_min is None: m_min   = minimum_source_mass(z, Tvir=1e4)
	var_min = barkana_loeb_2001.matter_variance(z, M=m_min, power_spectra=None)
	return delta_c - np.sqrt(2*(var_min-var_m))*K_zeta

def linear_fit_barrier(m, z, zeta=40, var_m=None):
	B0 = barrier_density(z, zeta=zeta, var_m=0)
	B1 = (barrier_density(z, zeta=zeta, var_m=0.1) - barrier_density(z, zeta=zeta, var_m=0))/0.1
	if var_m is None: var_m   = barkana_loeb_2001.matter_variance(z, M=m, power_spectra=None)
	return B0+B1*var_m

def mass_function(m, z, zeta=40):
	rho_mean_ = mean_density(z)                                        #Solar_mass/Mpc^3
	B0   = barrier_density(z, zeta=zeta, var_m=0)
	def f(x):
		m = np.exp(x)
		y = barkana_loeb_2001.matter_variance(z, M=m)
		return np.sqrt(y)
	
	if type(m)==np.float:
		sig2 = barkana_loeb_2001.matter_variance(z, M=m, power_spectra=None) 
		Bm   = linear_fit_barrier(m, z, zeta=zeta, var_m=None)
		dsigdm = derivative(f, np.log(m), order=15)
	else: 
		sig2 = np.array([barkana_loeb_2001.matter_variance(z, M=m1, power_spectra=None) for m1 in m])
		Bm   = np.array([linear_fit_barrier(m1, z, zeta=zeta, var_m=None) for m1 in m])
		dsigdm = np.array([derivative(f, np.log(m1), order=15) for m1 in m])
	m_dndm = np.sqrt(2/np.pi) * rho_mean_* np.abs(dsigdm) * B0 * np.exp(-Bm**2/(2*sig2))/sig2/m
	return m_dndm

def bubble_distribution(r, z, zeta=40):
	m = radius_to_mass(r,z)
	mass_func = mass_function(m, z, zeta=zeta)
	Vdn_dlogr = 4*np.pi*r**3*mass_func
	return Vdn_dlogr

def filling_factor(z, rs=None, zeta=40):
	#m_max   = radius_to_mass(r_max,z)
	#m_min   = minimum_source_mass(z, Tvir=1e4)
	#assert (m_max>m_min), "The radius is smaller that the minimum source radius."
	#mm = 10**np.linspace( np.log10(m_min), np.log10(m_max), 1000)
	if rs is None: 
		m_min = minimum_source_mass(z, Tvir=1e4)
		mm    = 10**np.linspace(np.log10(m_min), 20, 1000)
	else: mm = radius_to_mass(rs, z)	
	Vm = mm/mean_density(z)
	mm_dndm = mass_function(mm, z, zeta=zeta)
	dndm    = mm_dndm/mm
	Q = np.trapz(dndm*Vm, mm)
	return Q

def mean_density(z):
	rho_mean  = cosm_cons.Omega_m*barkana_loeb_2001.critical_density(z) #kg/m^3
	return rho_mean*1.47e37                                             #Solar_mass/Mpc^3
	
def mass_to_radius(m, z):
	rho_mean_ = mean_density(z)
	return (m*3/(4*np.pi*rho_mean_))**0.33			            # in Mpc

def radius_to_mass(r, z):
	rho_mean_ = mean_density(z)                             
	return 4*np.pi*r**3*rho_mean_/3				            # in Solar_mass

def minimum_source_mass(z, Tvir=1e4):
	return barkana_loeb_2001.Tvir_to_mass(Tvir, z)

def Vm0(m, z, r12=0):
	R  = mass_to_radius(m, z)
	v0 = 4*np.pi*R**3/3 - np.pi*r12*(R**2-r12**2/12)
	v0[r12>2*R] = 0 
	return v0


