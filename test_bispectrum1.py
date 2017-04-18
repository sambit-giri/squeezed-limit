import numpy as np
import c2raytools as c2t
import squeezed_bispectrum
import matplotlib.pyplot as plt
import matplotlib
import owntools

c2t.set_sim_constants(500)

zs  = [8.172, 7.570, 7.305, 7.139, 7.059, 6.905, 6.686]
xvs = [0.09, 0.20, 0.30, 0.35, 0.40, 0.50, 0.70]
zs  = [20.134, 8.515, 7.570, 7.305, 7.059, 6.905, 6.686]
colors = ['b', 'r', 'g', 'k', 'c', 'm', 'y']

xfrac_dir = '/disk/dawn-1/garrelt/Reionization/C2Ray_WMAP7/500Mpc/500Mpc_f2_0_300/results/'
dens_dir  = '/disk/dawn-1/garrelt/Reionization/C2Ray_WMAP7/500Mpc/coarser_densities/nc300/'

Ncuts = 3

plt.figure(1, figsize=(8, 12))

for i in xrange(len(zs)):
	z = zs[i]
	cube_21 = squeezed_bispectrum.make_21cm(xfrac_dir, dens_dir, z)
	cube_m = owntools.coeval_overdens(dens_dir, z)
	ibk_m_m, k_m_m = squeezed_bispectrum._integrated_bispectrum_normalized_cross(cube_m, cube_m, Ncuts=Ncuts)
	ibk_21_m, k_21_m = squeezed_bispectrum._integrated_bispectrum_normalized_cross(cube_21, cube_m, Ncuts=Ncuts)
	ibk_21_21, k_21_21 = squeezed_bispectrum._integrated_bispectrum_normalized_cross(cube_21, cube_21, Ncuts=Ncuts)
	### <m m m>
	plt.subplot(311)
	plt.ylabel('$iB_{\delta,\delta}(k)/P_{\delta}/\sigma_{\delta}^2$')
	plt.xlabel('k (Mpc$^{-1}$)')
	plt.plot(k_m_m, ibk_m_m, label='z='+str(z), c=colors[i])
	plt.xlim(0.08,0.6)
	plt.ylim(0,4)
	### <m 21 21>
	plt.subplot(312)
	plt.ylabel('$iB_{\delta,21}(k)/P_{21}/\sigma_{\delta}^2$')
	plt.xlabel('k (Mpc$^{-1}$)')
	plt.plot(k_21_m, ibk_21_m, c=colors[i])
	plt.xlim(0.08,0.6)
	plt.ylim(-20,10)
	### <21 21 21>
	plt.subplot(313)
	plt.ylabel('$iB_{21,21}(k)/P_{21}/\sigma_{21}^2$')
	plt.xlabel('k (Mpc$^{-1}$)')
	#k_21_21, ibk_21_21 = k_21_21[np.isfinite(ibk_21_21)], ibk_21_21[np.isfinite(ibk_21_21)]
	#mean_amp = ibk_21_21.mean(dtype=np.float64)
	plt.plot(k_21_21, ibk_21_21, c=colors[i])
	plt.xlim(0.08,0.6)
	#plt.ylim(-8e-5, 4e-5)

plt.subplot(311); plt.legend(loc=0)
plt.subplot(313); plt.legend(loc=0)
plt.savefig('all_combination_response_function.png')



