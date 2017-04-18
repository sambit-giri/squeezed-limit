import numpy as np
import c2raytools as c2t
import squeezed_bispectrum
import matplotlib.pyplot as plt
import owntools

from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.family'] = 'small-caps'
plt.rcParams['text.usetex'] = True
rcParams['axes.labelsize']=14
rcParams['font.size']=14 
rcParams['axes.linewidth'] = 1.2

c2t.set_sim_constants(500)

xfrac_dir = '/disk/dawn-1/garrelt/Reionization/C2Ray_WMAP7/500Mpc/500Mpc_f2_0_300/results/'
dens_dir  = '/disk/dawn-1/garrelt/Reionization/C2Ray_WMAP7/500Mpc/coarser_densities/nc300/'
ph_count_info = owntools.photoncount_info(xfrac_dir)

Ncuts = 3

xvs_ = 10**np.linspace(np.log10(3.44e-10),np.log10(0.20),10)
#xvs  = ph_count_info[[np.abs(ph_count_info[:,-2]-x).argmin() for x in xvs_],-2]
zs_  = ph_count_info[[np.abs(ph_count_info[:,-2]-x).argmin() for x in xvs_], 0]

dens_zs = owntools.get_zs_list(dens_dir, file_type='/*n_all.dat')
zs   = dens_zs[[np.abs(dens_zs-i).argmin() for i in zs_]]
xvs  = ph_count_info[[np.abs(ph_count_info[:,0]-i).argmin() for i in zs],-2]

i = 8
z = zs[i]

cube_21 = owntools.coeval_21cm(xfrac_dir, dens_dir, z, mean_subtract=True)
cube_m  = owntools.coeval_overdens(dens_dir, z)
P_dd, ks_m = c2t.power_spectrum_1d(cube_m, kbins=100, box_dims=c2t.conv.LB)
P_21, ks_x = c2t.power_spectrum_1d(cube_21, kbins=100, box_dims=c2t.conv.LB)

f_dd, k_dd  = squeezed_bispectrum._integrated_bispectrum_normalized_cross(cube_m, cube_m, Ncuts=Ncuts)
f_21d, k_xd = squeezed_bispectrum._integrated_bispectrum_normalized_cross(cube_21, cube_m, Ncuts=Ncuts)

### Model
zeta = 2
dx_ddelta = -0.0145128*zeta     # Get from excursion set calculation
xHI = 1. - xvs[i]
A = 2*xHI*dx_ddelta + xHI**2*f_dd
T_b_ = 23.88
dP21_ddelta = T_b_**2*A*P_dd
f_21d_ = dP21_ddelta/P_21

