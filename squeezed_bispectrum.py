import numpy as np
import c2raytools as c2t
import owntools
import scipy
import matplotlib
import matplotlib.pyplot as plt

def integrated_bispectrum(cube, Ncuts=4, statistic='mean', kbins=100):
	assert statistic in ['mean']
	assert cube.shape[0]%Ncuts==0 and cube.shape[1]%Ncuts==0 and cube.shape[2]%Ncuts==0
	Lx,Ly,Lz = cube.shape[0]/Ncuts,cube.shape[1]/Ncuts,cube.shape[2]/Ncuts
	smaller_cubes = [cube[i*Lx:(i+1)*Lx,j*Ly:(j+1)*Ly,k*Lz:(k+1)*Lz] for i in xrange(Ncuts) for j in xrange(Ncuts) for k in xrange(Ncuts)]
	B_k = np.zeros(kbins)
	for c in smaller_cubes:
		bk, ks = c2t.power_spectrum_1d(c, kbins=kbins, box_dims=c2t.conv.LB/Ncuts)
		B_k += bk*c.mean(dtype=np.float64)
	B_k = B_k/Ncuts**3
	return B_k, ks

def integrated_bispectrum_normalized(cube, Ncuts=4, statistic='mean', kbins=100):
	assert statistic in ['mean']
	assert cube.shape[0]%Ncuts==0 and cube.shape[1]%Ncuts==0 and cube.shape[2]%Ncuts==0
	Lx,Ly,Lz = cube.shape[0]/Ncuts,cube.shape[1]/Ncuts,cube.shape[2]/Ncuts
	smaller_cubes = [cube[i*Lx:(i+1)*Lx,j*Ly:(j+1)*Ly,k*Lz:(k+1)*Lz] for i in xrange(Ncuts) for j in xrange(Ncuts) for k in xrange(Ncuts)]
	B_k  = np.zeros(kbins, dtype=np.float64)
	P_k  = np.zeros(kbins, dtype=np.float64)
	sig2 = 0
	for c in smaller_cubes:
		pk, ks = c2t.power_spectrum_1d(c, kbins=kbins, box_dims=c2t.conv.LB/Ncuts)
		B_k  += pk*c.mean(dtype=np.float64)
		P_k  += pk
		sig2 += c.mean(dtype=np.float64)**2	#c.var(dtype=np.float64)
	B_k  = B_k/Ncuts**3
	P_k  = P_k/Ncuts**3
	sig2 = sig2/Ncuts**3
	return B_k/P_k/sig2, ks

def integrated_bispectrum_normalized_cross(cube, cube2, Ncuts=4, statistic='mean', kbins=100):
	assert statistic in ['mean']
	assert cube.shape[0]%Ncuts==0 and cube.shape[1]%Ncuts==0 and cube.shape[2]%Ncuts==0
	Lx,Ly,Lz = cube.shape[0]/Ncuts,cube.shape[1]/Ncuts,cube.shape[2]/Ncuts
	smaller_cubes = [cube[i*Lx:(i+1)*Lx,j*Ly:(j+1)*Ly,k*Lz:(k+1)*Lz] for i in xrange(Ncuts) for j in xrange(Ncuts) for k in xrange(Ncuts)]
	smaller_cubes2 = [cube2[i*Lx:(i+1)*Lx,j*Ly:(j+1)*Ly,k*Lz:(k+1)*Lz] for i in xrange(Ncuts) for j in xrange(Ncuts) for k in xrange(Ncuts)]
	B_k  = np.zeros(kbins, dtype=np.float64)
	P_k  = np.zeros(kbins, dtype=np.float64)
	sig2 = 0
	for i in xrange(Ncuts**3):
		c  = smaller_cubes[i]
		c2 = smaller_cubes2[i]
		pk, ks = c2t.power_spectrum_1d(c, kbins=kbins, box_dims=c2t.conv.LB/Ncuts)
		B_k  += pk*c2.mean(dtype=np.float64)
		P_k  += pk
		sig2 += c2.mean(dtype=np.float64)**2   #c2.var(dtype=np.float64)
	B_k  = B_k/Ncuts**3
	P_k  = P_k/Ncuts**3
	sig2 = sig2/Ncuts**3
	return B_k/P_k/sig2, ks

def _integrated_bispectrum_normalized_cross(cube, cube2, Ncuts=4, kbins=20, box_dims=None):
	#assert statistic in ['mean']
	assert cube.shape == cube2.shape
	assert cube.shape[0]%Ncuts==0 and cube.shape[1]%Ncuts==0 and cube.shape[2]%Ncuts==0
	Lx,Ly,Lz = cube.shape[0]/Ncuts,cube.shape[1]/Ncuts,cube.shape[2]/Ncuts
	rLs = [[Lx/2.+i*Lx,Ly/2.+j*Ly,Lz/2.+k*Lz] for i in xrange(Ncuts) for j in xrange(Ncuts) for k in xrange(Ncuts)]
	B_k   = np.zeros(kbins, dtype=np.float64)
	P_k   = np.zeros(kbins, dtype=np.float64)
	sig2  = 0
	n_box = Ncuts**3
	V_L   = (Lx*Ly*Lz)
	for i in xrange(n_box):
		w  = W_L(cube, rLs[i], [Lx,Ly,Lz])
		w2 = W_L(cube2, rLs[i], [Lx,Ly,Lz])
		c  = cube  * w
		c2 = cube2 * w2
		pk, ks = c2t.power_spectrum_1d(c, kbins=kbins, box_dims=box_dims)
		d_mean = c2.sum(dtype=np.float64)/V_L
		B_k   += pk*d_mean
		P_k   += pk
		sig2  += (d_mean)**2   #c2.var(dtype=np.float64)
		print 100*(i+1)/n_box, "%"
	B_k  = B_k/n_box
	P_k  = P_k/n_box
	sig2 = sig2/n_box
	return B_k/P_k/sig2, ks

def _integrated_bispectrum_normalized_cross1(cube, cube2, Ncuts=4):
	#assert statistic in ['mean']
	assert cube.shape == cube2.shape
	assert cube.shape[0]%Ncuts==0 and cube.shape[1]%Ncuts==0 and cube.shape[2]%Ncuts==0
	Lx,Ly,Lz = cube.shape[0]/Ncuts,cube.shape[1]/Ncuts,cube.shape[2]/Ncuts
	rLs = [[Lx/2.+i*Lx,Ly/2.+j*Ly,Lz/2.+k*Lz] for i in xrange(Ncuts) for j in xrange(Ncuts) for k in xrange(Ncuts)]
	B_ = 0
	P_ = 0
	sig2  = 0
	n_box = Ncuts**3
	V_L   = (Lx*Ly*Lz)
	for i in xrange(n_box):
		w  = W_L(cube, rLs[i], [Lx,Ly,Lz])
		w2 = W_L(cube2, rLs[i], [Lx,Ly,Lz])
		c  = cube  * w
		c2 = cube2 * w2
		x_mean = c.sum(dtype=np.float64)/V_L
		d_mean = c2.sum(dtype=np.float64)/V_L
		B_ += x_mean*d_mean
		P_ += x_mean
		sig2  += (d_mean)**2   #c2.var(dtype=np.float64)
		print 100*(i+1)/n_box, "%"
	B_ = B_/n_box
	P_ = P_/n_box
	sig2 = sig2/n_box
	return B_/P_/sig2

def _integrated_bispectrum_cross(cube, cube2, Ncuts=4, kbins=100):
	#assert statistic in ['mean']
	assert cube.shape == cube2.shape
	assert cube.shape[0]%Ncuts==0 and cube.shape[1]%Ncuts==0 and cube.shape[2]%Ncuts==0
	Lx,Ly,Lz = cube.shape[0]/Ncuts,cube.shape[1]/Ncuts,cube.shape[2]/Ncuts
	rLs = [[Lx/2.+i*Lx,Ly/2.+j*Ly,Lz/2.+k*Lz] for i in xrange(Ncuts) for j in xrange(Ncuts) for k in xrange(Ncuts)]
	B_k   = np.zeros(kbins, dtype=np.float64)
	P_k   = np.zeros(kbins, dtype=np.float64)
	sig2  = 0
	n_box = Ncuts**3
	V_L   = (Lx*Ly*Lz)
	for i in xrange(n_box):
		w  = W_L(cube, rLs[i], [Lx,Ly,Lz])
		w2 = W_L(cube2, rLs[i], [Lx,Ly,Lz])
		c  = cube  * w
		c2 = cube2 * w2
		pk, ks = c2t.power_spectrum_1d(c, kbins=kbins)
		d_mean = c2.sum(dtype=np.float64)/V_L
		B_k   += pk*d_mean
		print 100*(i+1)/n_box, "%"
	B_k  = B_k/n_box
	return B_k, ks

def W_L(array, rL, L):
	assert array.ndim == np.array(rL).size
	out = np.zeros(array.shape)
	if np.array(L).size==1:out[(rL[0]-L/2):(rL[0]+L/2),(rL[1]-L/2):(rL[1]+L/2),(rL[2]-L/2):(rL[2]+L/2)] = 1
	else:out[(rL[0]-L[0]/2):(rL[0]+L[0]/2),(rL[1]-L[1]/2):(rL[1]+L[1]/2),(rL[2]-L[2]/2):(rL[2]+L[2]/2)] = 1
	return out

def apply_integrated_bispectrum(xfrac_dir, dens_dir, zs, kbins=100, Ncuts=4, field='21cm'):
	B_k_s = np.zeros((kbins,len(zs)))
	k_s   = np.zeros((kbins,len(zs)))
	for i in xrange(len(zs)):
		if field=='21cm': cube = owntools.coeval_21cm(xfrac_dir, dens_dir, zs[i])
		elif field=='matter' or field=='density': 
			cube = owntools.coeval_dens(dens_dir, zs[i])
			cube = cube/cube.mean(dtype=np.float64) - 1.
		B_k_s[:,i], k_s[:,i] = integrated_bispectrum_normalized(cube, Ncuts=Ncuts, kbins=kbins)
		print "Squeezed-limit bispectrum has been calculated for redshift =", zs[i]
	return B_k_s, k_s
	

def position_dependent_powerspectra(xfrac_dir, dens_dir, z, ps_quantity='signal',quantity='density', statistic='mean', kbins=100, Ncuts=4):
	quantities = ['density', 'xfrac', 'signal']
	stats = ['mean', 'skewness', 'kurtosis']
	assert quantity in quantities
	assert statistic in stats
	if ps_quantity=='signal': cube = owntools.coeval_21cm(xfrac_dir, dens_dir, z)
	elif ps_quantity=='density': cube = owntools.coeval_dens(dens_dir, z)
	elif ps_quantity=='xfrac': cube = owntools.coeval_xfrac(xfrac_dir, z)
	stat_functions  = [np.mean, scipy.stats.skew, scipy.stats.kurtosis]
	P_k_s = np.zeros((kbins,Ncuts**3))
	k_s   = np.zeros((kbins,Ncuts**3))
	stat  = np.zeros(Ncuts**3)
	stat_func = stat_functions[stats.index(statistic)]
	if quantity == 'density':  quant = owntools.coeval_dens(dens_dir, z); quant = quant/quant.mean(dtype=np.float64)-1.
	elif quantity == 'xfrac':  quant = owntools.coeval_xfrac(xfrac_dir, z)
	elif quantity == 'signal': quant = owntools.coeval_21cm(xfrac_dir, dens_dir, z)
	Lx,Ly,Lz = cube.shape[0]/Ncuts,cube.shape[1]/Ncuts,cube.shape[2]/Ncuts
	smaller_cubes = [cube[i*Lx:(i+1)*Lx,j*Ly:(j+1)*Ly,k*Lz:(k+1)*Lz] for i in xrange(Ncuts) for j in xrange(Ncuts) for k in xrange(Ncuts)]
	smaller_quant = [quant[i*Lx:(i+1)*Lx,j*Ly:(j+1)*Ly,k*Lz:(k+1)*Lz] for i in xrange(Ncuts) for j in xrange(Ncuts) for k in xrange(Ncuts)]
	for i in xrange(len(smaller_cubes)):
		pk, ks = c2t.power_spectrum_1d(smaller_cubes[i], kbins=kbins, box_dims=c2t.conv.LB/Ncuts)
		P_k_s[:,i], k_s[:,i] = pk, ks
		stat[i] = stat_func(smaller_quant[i])
	return P_k_s, k_s, stat


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
		plt.loglog(x, y, color=s_m.to_rgba(parameters[i]))
	Y1 = Y.copy()
	Y1[Y1!=Y1] = np.median(Y1)
	plt.xlim(0.2, X.max())
	plt.ylim(Y1.min(), 2e5)

	plt.colorbar(s_m, label=label)
	return 0


def cross_correlation_power_spectra(fielda, fieldb, kbins=100, box_dims=None):
	Pab, kab = c2t.cross_power_spectrum_1d(fielda, fieldb, kbins=kbins, box_dims=box_dims)
	Paa, kaa = c2t.power_spectrum_1d(fielda, kbins=kbins, box_dims=box_dims)
	Pbb, kbb = c2t.power_spectrum_1d(fieldb, kbins=kbins, box_dims=box_dims)
	Rxx = Pab/np.sqrt(Paa*Pbb)
	return Rxx, kab

def position_dependent_cross_correlation_PS(fieldx, fieldy, xfrac_dir, dens_dir, z, quantity='density', kbins=100, Ncuts=4, statistic='mean'):
	assert fieldx.ndim == fieldy.ndim 

	stats = ['mean', 'skewness', 'kurtosis']
	stat_functions  = [np.mean, scipy.stats.skew, scipy.stats.kurtosis]
	stat_func = stat_functions[stats.index(statistic)]
	R_k_s = np.zeros((kbins,Ncuts**3))
	k_s   = np.zeros((kbins,Ncuts**3))
	stat  = np.zeros(Ncuts**3)

	if quantity == 'density':  quant = owntools.coeval_dens(dens_dir, z); quant = quant/quant.mean(dtype=np.float64)-1.
	elif quantity == 'xfrac':  quant = owntools.coeval_xfrac(xfrac_dir, z)
	elif quantity == 'signal': quant = owntools.coeval_21cm(xfrac_dir, dens_dir, z)

	Lx,Ly,Lz = np.array(fieldx.shape)/Ncuts
	smaller_x = [fieldx[i*Lx:(i+1)*Lx,j*Ly:(j+1)*Ly,k*Lz:(k+1)*Lz] for i in xrange(Ncuts) for j in xrange(Ncuts) for k in xrange(Ncuts)]
	smaller_y = [fieldy[i*Lx:(i+1)*Lx,j*Ly:(j+1)*Ly,k*Lz:(k+1)*Lz] for i in xrange(Ncuts) for j in xrange(Ncuts) for k in xrange(Ncuts)]
	smaller_quant = [quant[i*Lx:(i+1)*Lx,j*Ly:(j+1)*Ly,k*Lz:(k+1)*Lz] for i in xrange(Ncuts) for j in xrange(Ncuts) for k in xrange(Ncuts)]

	for i in xrange(len(smaller_x)):
		rk, ks = cross_correlation_power_spectra(smaller_x[i], smaller_y[i], kbins=kbins, box_dims=c2t.conv.LB/Ncuts)
		R_k_s[:,i], k_s[:,i] = rk, ks
		stat[i] = stat_func(smaller_quant[i])
	return R_k_s, k_s, stat

def make_21cm(xfrac_dir, dens_dir, z):
	dt = owntools.coeval_21cm(xfrac_dir, dens_dir, z)
	return dt-dt.mean()



