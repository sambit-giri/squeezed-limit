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

#zs  = [20.134, 17.848, 15.596, 12.603, 10.11, 9.026]
zs = [20.134, 12.603, 9.026, 7.305, 6.549, 6.113]
colors = ['b', 'r', 'g', 'k', 'c', 'm']

xfrac_dir = '/disk/dawn-1/garrelt/Reionization/C2Ray_WMAP7/500Mpc/500Mpc_f2_0_300/results/'
dens_dir  = '/disk/dawn-1/garrelt/Reionization/C2Ray_WMAP7/500Mpc/coarser_densities/nc300/'

Ncuts = 3

plt.figure(1, figsize=(12, 12))

for i in xrange(len(zs)):
	z = zs[i]
	cube_21 = owntools.coeval_21cm(xfrac_dir, dens_dir, z, mean_subtract=True)
	cube_m = owntools.coeval_overdens(dens_dir, z)
	ibk_m_m, k_m_m = squeezed_bispectrum._integrated_bispectrum_normalized_cross(cube_m, cube_m, Ncuts=Ncuts)
	ibk_21_m, k_21_m = squeezed_bispectrum._integrated_bispectrum_normalized_cross(cube_21, cube_m, Ncuts=Ncuts)
	#ibk_21_21, k_21_21 = squeezed_bispectrum._integrated_bispectrum_normalized_cross(cube_21, cube_21, Ncuts=Ncuts)
	### <m m m> - <21 21 m>
	plt.subplot(111)
	plt.plot(k_m_m, ibk_m_m-ibk_21_m, label='z='+str(z), c=colors[i])

plt.ylabel('$\hat{f}_{<\delta,\delta,\delta>}-\hat{f}_{<21,21,\delta>}$')
plt.xlabel('k (Mpc$^{-1}$)')
plt.xlim(0.08,0.6)
plt.ylim(-22,35)

plt.legend(loc=0)
plt.savefig('response_function_difference.png')



