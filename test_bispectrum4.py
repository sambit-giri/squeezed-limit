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

P_k_m = []
P_k_x = []

## Global Power Spectrum
for i in xrange(len(zs)):
	z = zs[i]
	print z
	cube_x = owntools.coeval_21cm(xfrac_dir, dens_dir, z, mean_subtract=True)
	cube_m = owntools.coeval_overdens(dens_dir, z)
	pk_m, ks_m = c2t.power_spectrum_1d(cube_m, kbins=100, box_dims=c2t.conv.LB)
	pk_x, ks_x = c2t.power_spectrum_1d(cube_x, kbins=100, box_dims=c2t.conv.LB)
	P_k_m.append((pk_m, ks_m))
	P_k_x.append((pk_x, ks_x))

iB_k_mm = []
iB_k_xm = []
### Response Functions
for i in xrange(len(zs)):
	z = zs[i]
	print z
	cube_x = owntools.coeval_21cm(xfrac_dir, dens_dir, z, mean_subtract=True)
	cube_m = owntools.coeval_overdens(dens_dir, z)
	ibk_m_m, k_m_m = squeezed_bispectrum._integrated_bispectrum_normalized_cross(cube_m, cube_m, Ncuts=Ncuts)
	ibk_x_m, k_x_m = squeezed_bispectrum._integrated_bispectrum_normalized_cross(cube_x, cube_m, Ncuts=Ncuts)
	iB_k_mm.append((ibk_m_m, k_m_m))
	iB_k_xm.append((ibk_x_m, k_x_m))
	
### Save data
ii = [0,3,5,6,7,8,9]
for i in ii: np.savetxt('response_function_dd_z_'+str(zs[i])+'.txt', np.array([iB_k_mm[i][1],iB_k_mm[i][0]]).T)
for i in ii: np.savetxt('response_function_xd_z_'+str(zs[i])+'.txt', np.array([iB_k_xm[i][1],iB_k_xm[i][0]]).T)
for i in ii: np.savetxt('global_PS_dd_z_'+str(zs[i])+'.txt', np.array([P_k_m[i][1],P_k_m[i][0]]).T)
for i in ii: np.savetxt('global_PS_xx_z_'+str(zs[i])+'.txt', np.array([P_k_x[i][1],P_k_x[i][0]]).T)

### Make plots
plt.figure(1, figsize=(12, 12))
ii = [0,3,5,6,7,8,9]
plt.subplot(221)
for i in ii:plt.loglog(P_k_m[i][1],P_k_m[i][0], label='z='+str(zs[i])+', $<x_\mathrm{v}>$='+str(xvs[i]))
#plt.legend(loc=0)
plt.xlim(0.02,2.0)
plt.ylim(0.1,1500)
plt.xlabel('k (/Mpc)')
plt.ylabel('$P_{\delta \delta}$ (k)')
#plt.show()

#plt.figure(2, figsize=(12, 12))
plt.subplot(222)
for i in ii:plt.loglog(P_k_x[i][1],P_k_x[i][0], label='z='+str(zs[i])+', $<x_\mathrm{v}>$='+str(xvs[i]))
#plt.legend(loc=0)
plt.xlim(0.02,2.0)
plt.ylim(300,7e5)
plt.xlabel('k (/Mpc)')
plt.ylabel('$P_{x x}$ (k)')
#plt.show()

#plt.figure(3, figsize=(12, 12))
plt.subplot(223)
for i in ii:plt.semilogx(P_k_x[i][1],P_k_x[i][0]/P_k_m[i][0], label='z='+str(zs[i])+', $<x_\mathrm{v}>$='+str(xvs[i]))
#plt.legend(loc=2)
plt.xlim(0.02,2.0)
#plt.ylim(300,7e5)
plt.xlabel('k (/Mpc)')
plt.ylabel('$P_{\delta \delta}$/$P_{x x}$ (k)')

plt.subplot(224)
for i in ii:plt.semilogx([],[], label='z='+str(zs[i])+', $<x_\mathrm{v}>$='+str(xvs[i]))
plt.legend(loc=0)
plt.xticks([])
plt.yticks([])
#plt.xlim(0.02,2.0)
#plt.ylim(300,7e5)
#plt.xlabel('k (/Mpc)')
#plt.ylabel('$P_{\delta \delta}$/$P_{x x}$ (k)')

plt.show()


### Make plots
plt.figure(2, figsize=(16, 8))
ii = [0,3,5,6,7,8,9]
plt.subplot(121)
for i in ii:plt.semilogx(iB_k_mm[i][1],iB_k_mm[i][0], label='z='+str(zs[i]))
plt.legend(loc=0)
plt.xlim(0.1,0.9)
plt.ylim(-4,12)
plt.xlabel('k (/Mpc)')
plt.ylabel('$\mathrm{d} (ln P_{\delta \delta})$/$\mathrm{d} \overline{\delta}$ (k)')
#plt.show()

#plt.figure(2, figsize=(12, 12))
plt.subplot(122)
for i in ii:plt.semilogx(iB_k_xm[i][1],iB_k_xm[i][0], label='$<x_\mathrm{v}>$='+str(xvs[i]))
plt.legend(loc=0)
plt.xlim(0.1,0.9)
plt.ylim(-4,12)
plt.xlabel('k (/Mpc)')
plt.ylabel('$\mathrm{d} (ln P_{x \delta})$/$\mathrm{d} \overline{\delta}$ (k)')

plt.show()



