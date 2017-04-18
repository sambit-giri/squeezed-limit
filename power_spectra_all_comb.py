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

plt.figure()
for i in xrange(len(zs)):
	z = zs[i]
	cube_x = owntools.coeval_xfrac(xfrac_dir, z); cube_x = cube_x - cube_x.mean()
	cube_m  = owntools.coeval_overdens(dens_dir, z)
	P_mm, k_mm = c2t.power_spectrum_1d(cube_m, kbins=100)
	P_xx, k_xx = c2t.power_spectrum_1d(cube_x, kbins=100)
	P_mx, k_mx = c2t.cross_power_spectrum_1d(cube_x, cube_m, kbins=100)
	plt.subplot(2,3,i+1)
	plt.title('z='+str(z))
	plt.loglog(k_mm, P_mm, c='b')
	plt.loglog(k_xx, P_xx, c='r')
	plt.loglog(k_mx, P_mx, c='g')
	print 'z=', z, 'done'

plt.suptitle('BLUE: $\delta \delta$, RED: $x_n x_n$, GREEN: $\delta x_n$')
for i in xrange(len(zs)):
	plt.subplot(2,3,i+1)
	plt.ylim(2e-4,2e4)
	plt.xlabel('k')
	plt.ylabel('P(k)')




