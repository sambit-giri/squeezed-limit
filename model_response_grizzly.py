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

xfrac_dir  = '/disk/dawn-1/garrelt/Reionization/C2Ray_WMAP7/500Mpc/500Mpc_f2_0_300/results/'
grizzly_dir = '/disk/dawn-1/gragh/sambit_bs/outputs/'
dens_dir    = '/disk/dawn-1/garrelt/Reionization/C2Ray_WMAP7/500Mpc/coarser_densities/nc600/'
ph_count_info = owntools.photoncount_info(xfrac_dir)

Ncuts = 3

xvs_ = 10**np.linspace(np.log10(3.44e-10),np.log10(0.20),10)
#xvs  = ph_count_info[[np.abs(ph_count_info[:,-2]-x).argmin() for x in xvs_],-2]
zs_  = ph_count_info[[np.abs(ph_count_info[:,-2]-x).argmin() for x in xvs_], 0]

dens_zs = owntools.get_zs_list(dens_dir, file_type='/*n_all.dat')
zs   = dens_zs[[np.abs(dens_zs-i).argmin() for i in zs_]]
xvs  = ph_count_info[[np.abs(ph_count_info[:,0]-i).argmin() for i in zs],-2]

z = 9.026 #8.636

gg_xf = c2t.XfracFile(grizzly_dir+str(z)+'xhiifrac.dat').xi
#cc_xf = owntools.coeval_xfrac(xfrac_dir, z)

cube_d  = owntools.coeval_dens(dens_dir, z)
cube_21 = c2t.calc_dt(gg_xf, cube_d, z); cube_21 -= cube_21.mean()
cube_m  = cube_d/cube_d.mean(dtype=np.float64) - 1. #owntools.coeval_overdens(dens_dir, z)

P_dd, ks_m = c2t.power_spectrum_1d(cube_m, kbins=100, box_dims=c2t.conv.LB)
P_21, ks_x = c2t.power_spectrum_1d(cube_21, kbins=100, box_dims=c2t.conv.LB)

f_dd, k_dd  = squeezed_bispectrum._integrated_bispectrum_normalized_cross(cube_m, cube_m, Ncuts=Ncuts)
f_21d, k_xd = squeezed_bispectrum._integrated_bispectrum_normalized_cross(cube_21, cube_m, Ncuts=Ncuts)
f_xx, k_xx  = squeezed_bispectrum._integrated_bispectrum_normalized_cross(cube_x-cube_x.mean(), cube_x-cube_x.mean(), Ncuts=Ncuts)
f_x_ = squeezed_bispectrum._integrated_bispectrum_normalized_cross1(1-cube_x, cube_m, Ncuts=Ncuts)

ks = k_dd.copy()

### Zeta calculation
sources = np.loadtxt(parent_dir+'sources/'+str(z)+'-coarser_sources.dat', skiprows=1)
M_min = sources[np.nonzero(sources[:,-2]),-2].min()*c2t.conv.M_grid*c2t.const.solar_masses_per_gram #solarunit
M_max = sources[np.nonzero(sources[:,-2]),-2].max()*c2t.conv.M_grid*c2t.const.solar_masses_per_gram #solarunit
M_halo_sum = sources[:,-2].sum()*c2t.conv.M_grid*c2t.const.solar_masses_per_gram                    #solarunit
mpc_to_cm  = 3.086e24
den_halo_sum  = M_halo_sum/(c2t.conv.LB*mpc_to_cm)**3  # solarmass/cm^3
den_halo_sum /= c2t.const.solar_masses_per_gram        # gm/cm^3
zeta = (cube_x*cube_d).sum()/cube_x.sum()/den_halo_sum

zs_sim = ph_count_info[:,0]
xm_sim = ph_count_info[:,-1]
xv_sim = ph_count_info[:,-2]

### Model
zz = 9.026 #8.636
xm = xm_sim[zs_sim==zz][0]
xv = gg_xf.mean()
command = './Anson/Sambit/src/ESMR.x '+str(zz)+' '+str(xv)+' info.txt'
os.system(command)
info = np.loadtxt('./info.txt')
os.remove('./info.txt')
plt.clf()
#zeta = 25 
cc = 'g'
zeta = info[1]#*xv/(xv-xv_sim[zs_sim>zz][-1])
df_ddelta = info[0]
dx_ddelta = df_ddelta*zeta     # Get from excursion set calculation
xHI = 1. - xv
A = 2*xHI*dx_ddelta + xHI**2*f_dd
T_b_ = 27.
dP21_ddelta = T_b_**2*A #*P_dd
f_21d_ = A.copy() #+ f_xx #/P_21

plt.semilogx(ks, f_dd, c='k', label='$f_{\delta \delta}$')
plt.semilogx(ks, f_21d, c='b', label='$f_{21 \delta}$')
plt.semilogx(ks, f_21d_, '--', c=cc, label='$f_{21 \delta, model}$($\zeta$='+str(zeta)+')')
plt.legend(loc=0)
plt.xlim(0.02,1.)
plt.ylim(-2.5,6.0)
plt.xlabel('k (Mpc$^{-1}$)')
plt.ylabel('f')
plt.show()






