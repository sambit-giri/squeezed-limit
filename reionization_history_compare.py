import numpy as np
import c2raytools as c2t
import squeezed_bispectrum
import matplotlib.pyplot as plt
import owntools
import os

from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.family'] = 'small-caps'
plt.rcParams['text.usetex'] = True
rcParams['axes.labelsize']=14
rcParams['font.size']=14 
rcParams['axes.linewidth'] = 1.2

c2t.set_sim_constants(500) 

parent_dir = '/disk/dawn-1/garrelt/Reionization/C2Ray_WMAP7/500Mpc/'
xfrac_dir  = '/disk/dawn-1/garrelt/Reionization/C2Ray_WMAP7/500Mpc/500Mpc_f2_0_300/results/'
dens_dir   = '/disk/dawn-1/garrelt/Reionization/C2Ray_WMAP7/500Mpc/coarser_densities/nc300/'
ph_count_info = owntools.photoncount_info(xfrac_dir)

zs_sim = ph_count_info[:,0]
xv_sim = ph_count_info[:,-2]
xm_sim = ph_count_info[:,-1]

def from_excursion_set(zz, xv):
	command = './Anson/Sambit/src/ESMR.x '+str(zz)+' '+str(xv)+' info.txt'
	os.system(command)
	info = np.loadtxt('./info.txt')
	os.remove('./info.txt')
	return info

xv_exc = np.array([from_excursion_set(zz, 1)[2] for zz in zs_sim])
zetas  = xv_sim/xv_exc


### Idea 16:22 25-04-2017
zs_check = zs_sim[zs_sim>=zz]
f_s = np.array([from_excursion_set(zz, 1)[2] for zz in zs_check])
df_ddelta_s = np.array([from_excursion_set(zz, 1)[0] for zz in zs_check])
