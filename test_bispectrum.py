import numpy as np
import c2raytools as c2t
import squeezed_bispectrum
import matplotlib.pyplot as plt
import matplotlib

def color_plot_lines(X, Y, parameters, cmap='jet', label='$\delta$', norm=None, ylim_top=None):
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
		plt.loglog(x, y, color=s_m.to_rgba(parameters[i]))
	Y1 = Y.copy()
	Y1[Y1!=Y1] = np.median(Y1)
	plt.xlim(0.7, X.max())
	if not ylim_top: plt.ylim(Y1.min(), Y1.max())
	else: plt.ylim(Y1.min(), ylim_top)

	plt.colorbar(s_m, label=label)
	return 0

xfrac_dir = '/disk/dawn-1/sgiri/simulations/244Mpc/244Mpc_f2_0_250/results'
dens_dir  = '/disk/dawn-1/sgiri/simulations/244Mpc/coarser_densities/nc250/'
#zs  = [8.172, 7.570, 6.905, 6.686];xvs = [0.09, 0.20, 0.50, 0.70]
zs  = [7.391, 7.263, 7.180, 7.059];xvs = [0.25, 0.30, 0.34, 0.40]
#zs  = [7.221, 7.180, 7.139, 7.099];xvs = [0.32, 0.34, 0.36, 0.38]
Ncuts = 5
quantity = 'density'
ps_quantity = 'signal'

if quantity == 'density':
	norm = matplotlib.colors.Normalize(vmin=-0.091, vmax=0.114)
	label = '$\delta$'
elif quantity == 'xfrac':
	norm = matplotlib.colors.Normalize(vmin=0.049, vmax=0.956)
	label = '$x_{HII}$'
elif quantity == 'signal':
	norm = matplotlib.colors.Normalize(vmin=-5.1, vmax=5.1)
	label = '$\delta T_b$'
if ps_quantity == 'density': ylim_top = 5e-60; ps_label='$P_{k,\delta}$'
elif ps_quantity == 'signal': ylim_top = 2e5; ps_label='$P_{k,\delta T}$'
else: ylim_top = None; ps_label='$P_{k,x}$'

ii = 1
for z in zs:
	yy, xx, pp = squeezed_bispectrum.position_dependent_powerspectra(xfrac_dir, dens_dir, z, ps_quantity =ps_quantity ,quantity=quantity, statistic='mean', kbins=100, Ncuts=Ncuts)
	plt.subplot(2,2,ii)
	if quantity == 'signal': pp = pp-pp.mean()
	print pp.min(), pp.max()
	color_plot_lines(xx, yy, pp, cmap='jet', label=label, norm=norm, ylim_top=ylim_top)
	plt.title('$x_v$='+str(xvs[zs.index(z)]))
	plt.ylabel(ps_label)
	plt.xlabel('k (Mpc$^{-1}$)')
	ii += 1
"""
## Plot evolution of variances in xfracs
all_zs  = owntools.get_zs_list(xfrac_dir)
cub_rms = np.zeros(all_zs.size)
sml_rms = []#np.zeros((all_zs.size,125))
for i in xrange(all_zs.size):
	z = all_zs[i]                                    
	xf = owntools.coeval_xfrac(xfrac_dir, z)
	cub_rms[i] = xf.var()
	Lx,Ly,Lz = xf.shape[0]/Ncuts,xf.shape[1]/Ncuts,xf.shape[2]/Ncuts
	#smaller_cubes = [xf[i*Lx:(i+1)*Lx,j*Ly:(j+1)*Ly,k*Lz:(k+1)*Lz] for i in xrange(Ncuts) for j in xrange(Ncuts) for k in xrange(Ncuts)]
	smaller_rms = [xf[i*Lx:(i+1)*Lx,j*Ly:(j+1)*Ly,k*Lz:(k+1)*Lz].var() for i in xrange(Ncuts) for j in xrange(Ncuts) for k in xrange(Ncuts)]
	sml_rms.append(smaller_rms)
    for j in xrange(len(smaller_cubes)):
        sml_rms[i,j] = smaller_cubes[j].var()
"""

