import numpy as np
import c2raytools as c2t
import squeezed_bispectrum
import matplotlib.pyplot as plt
import matplotlib
import owntools

c2t.set_sim_constants(244)

def color_plot_lines(X, Y, parameters, cmap='jet', label='$\delta$', norm=None):
	cmaps = ['jet', 'cool', 'brg', 'hsv']
	cm_func = [matplotlib.cm.jet, matplotlib.cm.cool, matplotlib.cm.brg, matplotlib.cm.hsv]
	if not norm: norm = matplotlib.colors.Normalize(vmin=np.min(parameters), vmax=np.max(parameters))
	# choose a colormap
	c_m = cm_func[cmaps.index(cmap)]
	# create a ScalarMappable and initialize a data structure
	s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
	s_m.set_array([])
	for i in xrange(len(parameters)):
		a,b = X[:,i], Y[:,i]
		x = a[b!=np.nan]
		y = b[b!=np.nan]
		plt.semilogx(x, y, color=s_m.to_rgba(parameters[i]))
	plt.xlim(0.35, X.max())
	plt.ylim(-1, 1)
	plt.colorbar(s_m, label=label)
	return 0

zs  = [8.172, 7.570, 7.305, 7.139, 7.059, 6.905, 6.686]
xvs = [0.09, 0.20, 0.30, 0.35, 0.40, 0.50, 0.70]

xfrac_dir = '/disk/dawn-1/sgiri/simulations/244Mpc/244Mpc_f2_0_250/results'
dens_dir  = '/disk/dawn-1/sgiri/simulations/244Mpc/coarser_densities/nc250/'
zs  = [8.172, 7.570, 6.905, 6.686];xvs = [0.09, 0.20, 0.50, 0.70]
#zs  = [7.391, 7.263, 7.180, 7.059];xvs = [0.25, 0.30, 0.34, 0.40]
#zs  = [7.221, 7.180, 7.139, 7.099];xvs = [0.32, 0.34, 0.36, 0.38]
Ncuts = 5
quantity = 'density'

if quantity == 'density':
	norm = matplotlib.colors.Normalize(vmin=-0.091, vmax=0.114)
	label = '$\delta$'
elif quantity == 'xfrac':
	norm = matplotlib.colors.Normalize(vmin=0.049, vmax=0.956)
	label = '$x_{HII}$'
elif quantity == 'signal':
	norm = matplotlib.colors.Normalize(vmin=-5.1, vmax=5.1)
	label = '$\delta T_b$'

ii = 1
for z in zs:
	field_21 = owntools.coeval_21cm(xfrac_dir, dens_dir, z)
	field_d  = owntools.coeval_dens(dens_dir, z)
	yy, xx, pp = squeezed_bispectrum.position_dependent_cross_correlation_PS(field_21, field_d, xfrac_dir, dens_dir, z, Ncuts=5, quantity=quantity)
	plt.subplot(2,2,ii)
	if quantity == 'signal': pp = pp-pp.mean()
	print pp.min(), pp.max()
	color_plot_lines(xx, yy, pp, cmap='jet', label=label, norm=norm)
	plt.title('$x_v$='+str(xvs[zs.index(z)]))
	plt.ylabel('R$_{\delta,21cm}$')
	plt.xlabel('k (Mpc$^{-1}$)')
	ii += 1


for z in zs:
	cube = owntools.coeval_21cm(xfrac_dir,dens_dir, z)
	cube2 = owntools.coeval_overdens(dens_dir, z)
	ibc, kc = squeezed_bispectrum.integrated_bispectrum_normalized_cross(cube, cube2, Ncuts=Ncuts)
	plt.semilogx(kc, ibc, label='z='+str(z))

plt.legend()

for z in zs:
     cube1 = owntools.coeval_21cm(xfrac_dir,dens_dir, z)
     cube1 = cube1 - cube1.mean(dtype=np.float64)
     cube2 = owntools.coeval_overdens(dens_dir, z)
     ibc0, kc0 = squeezed_bispectrum.integrated_bispectrum_normalized_cross(cube1, cube1, Ncuts=Ncuts)
     ibc1, kc1 = squeezed_bispectrum.integrated_bispectrum_normalized_cross(cube1, cube2, Ncuts=Ncuts)
     ibc2, kc2 = squeezed_bispectrum.integrated_bispectrum_normalized_cross(cube2, cube2, Ncuts=Ncuts)
     plt.subplot(131)
     plt.plot(kc0[np.isfinite(ibc0)], ibc0[np.isfinite(ibc0)], label='$x_v$='+str(xvs[zs.index(z)]))
     plt.subplot(132)
     plt.plot(kc1[np.isfinite(ibc1)], ibc1[np.isfinite(ibc1)])
     plt.subplot(133)
     plt.plot(kc2[np.isfinite(ibc2)], ibc2[np.isfinite(ibc2)], label='z='+str(z))

plt.subplot(131); plt.legend(loc=0)
plt.xlim(0.25,1.); plt.ylim(-0.7,0.2); plt.title('21cm-21cm'); plt.xlabel('k [Mpc$^{-1}$]'); plt.ylabel('$f_{21,21}$')
plt.subplot(132)
plt.xlim(0.25,1.); plt.ylim(-8,10.); plt.title('21cm-$\delta$'); plt.xlabel('k [Mpc$^{-1}$]'); plt.ylabel('$f_{21,\delta}$')
plt.subplot(133); plt.legend(loc=0)
plt.xlim(0.25,1.); plt.ylim(2.,3.2); plt.title('$\delta$-$\delta$'); plt.xlabel('k [Mpc$^{-1}$]'); plt.ylabel('$f_{\delta,\delta}$')


