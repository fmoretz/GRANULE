import numpy as np
import matplotlib.pyplot as plt
from physical_coefficients import *
from influent_input import *
from design_input import *

npt =101
viscosity_H2O = 2.414/100*10**(247.8/(To - 140))   #cP
Dp  = Ro/1000*2
dr = Ro/1000/npt
_D = []; _BC = []; _km = [];
D_pos = []; BC_pos = []; km_pos = []

for k in range(0,6):

    D   = np.linspace(10**(-12+k), 10**(-11+k), 1000)
    Re  = rho_H2O * u/3600 * Dp /(viscosity_H2O/1000)
    Sc  = (viscosity_H2O/1000) / ( rho_H2O * D/3600)
    Sh  = 2 + ( 1.6*Re**(1/3) + 0.6*Re**(0.5) + 5e-3*Re**(0.8) ) * Sc**(1/3)

    km  = Sh * D / Dp

    BC = km - D/dr

    for j in range(0, len(D)-1):
        _D.append(D[j])
        _BC.append(BC[j])
        _km.append(km[j])

for i in range(0, len(_D)):
    if _BC[i] >=0:
        D_pos.append(_D[i])
        BC_pos.append(_BC[i])
        km_pos.append(_km[i])
    else:
        break
 
fig, axs =plt.subplots(2, 2, figsize=(6, 8))
plt.subplots_adjust(
	left   = 0.125,
	bottom = 0.071,
	right  = 0.9,
	top    = 0.971,
	wspace = 0.2,
	hspace = 0.5
	)
axs[0,0].plot(_D, _BC, 'ro-', markevery=1000)
axs[0,0].set_xlabel('Diffusivity - m/h')
axs[0,0].set_ylabel('BC - m2/h')
axs[0,0].set_title('Diffusivity vs BC', fontsize=10)
axs[0,0].grid(True)

axs[0,1].plot(_D, _km, 'bo-', markevery=1000)
axs[0,1].set_xlabel('Diffusivity - m/h')
axs[0,1].set_ylabel('km - m2/h')
axs[0,1].set_title('Diffusivity vs km', fontsize=10)
axs[0,1].grid(True)

axs[1,0].plot(D_pos, BC_pos, 'ro-', markevery=100)
axs[1,0].set_xlabel('Diffusivity - m/h')
axs[1,0].set_ylabel('BC - m2/h')
axs[1,0].set_title('Range of BC positivity - BC vs D', fontsize=10)
axs[1,0].grid(True)

axs[1,1].plot(D_pos, km_pos, 'bo-', markevery=100)
axs[1,1].set_xlabel('Diffusivity - m/h')
axs[1,1].set_ylabel('km - m2/h')
axs[1,1].set_title('Range of BC positivity - BC vs km', fontsize=10)
axs[1,1].grid(True)

plt.show()
