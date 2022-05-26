import numpy as np
import matplotlib.pyplot as plt

from sunpharma_500 import _S as S_500, S_ss as S_ss_500,  _X as X_500, X_ss as X_ss_500,  _E as E_500, E_ss as E_ss_500   
from sunpharma_500 import M_cumulative as M_cumulative_500, _M as M_500, M_ss as M_ss_500,  _M_vol as M_vol_500, M_vol_ss as M_vol_ss_500   
from sunpharma_500 import t as t_500, t_ss as t_ss_500

from sunpharma_1000 import _S as S_1000, S_ss as S_ss_1000,  _X as X_1000, X_ss as X_ss_1000,  _E as E_1000, E_ss as E_ss_1000   
from sunpharma_1000 import M_cumulative as M_cumulative_1000, _M as M_1000, M_ss as M_ss_1000,  _M_vol as M_vol_1000, M_vol_ss as M_vol_ss_1000   
from sunpharma_1000 import t as t_1000, t_ss as t_ss_1000

from sunpharma_1500 import _S as S_1500, S_ss as S_ss_1500,  _X as X_1500, X_ss as X_ss_1500,  _E as E_1500, E_ss as E_ss_1500   
from sunpharma_1500 import M_cumulative as M_cumulative_1500, _M as M_1500, M_ss as M_ss_1500,  _M_vol as M_vol_1500, M_vol_ss as M_vol_ss_1500   
from sunpharma_1500 import t as t_1500, t_ss as t_ss_1500

from sunpharma_2000 import _S as S_2000, S_ss as S_ss_2000,  _X as X_2000, X_ss as X_ss_2000,  _E as E_2000, E_ss as E_ss_2000   
from sunpharma_2000 import M_cumulative as M_cumulative_2000, _M as M_2000, M_ss as M_ss_2000,  _M_vol as M_vol_2000, M_vol_ss as M_vol_ss_2000   
from sunpharma_2000 import t as t_2000, t_ss as t_ss_2000


fig, axs = plt.subplots(4, 2, figsize=(6, 8))
plt.subplots_adjust(
	left   = 0.125,
	bottom = 0.071,
	right  = 0.9,
	top    = 0.971,
	wspace = 0.2,
	hspace = 0.5
	)

axs[1,0].plot(t_500,  S_500,  'b--', linewidth=0.4, label='S 500')
axs[1,0].plot(t_1000, S_1000, 'm--', linewidth=0.4, label='S 1000')
axs[1,0].plot(t_1500, S_1500, 'r--', linewidth=0.4, label='S 1500')
axs[1,0].plot(t_2000, S_2000, 'g--', linewidth=0.4, label='S 2000')
axs[1,0].plot(t_ss_500,  S_ss_500,  'b-o', linewidth=1)
axs[1,0].plot(t_ss_1000, S_ss_1000, 'm-o', linewidth=1)
axs[1,0].plot(t_ss_1500, S_ss_1500, 'r-o', linewidth=1)
axs[1,0].plot(t_ss_2000, S_ss_2000, 'g-o', linewidth=1)
axs[1,0].set_xlabel('t - d')
axs[1,0].set_ylabel('S, - kg/m3')
axs[1,0].set_xlim([0, t_2000[-1]+10])
axs[1,0].set_title('Substrate Concentration Profile', fontsize=10)
axs[1,0].grid(True)

	#axs[1,1].plot(t, _X, 'y--', linewidth=0.4, label='X')
	#axs[1,1].plot(t_ss, X_ss, 'y-o', linewidth=1)
	#axs[1,1].set_xlabel('t - d')
	#axs[1,1].set_ylabel('X - kg/m3')
	#axs[1,1].set_xlim([1000, t[-1]+10])
	#axs[1,1].set_ylim([0, 1.5*X0])
	#axs[1,1].set_title('Active Biomass Concentration Profile', fontsize=10)
	#axs[1,1].grid(True)
#
#
	#axs[2,0].plot(t, _E, 'm--', linewidth=0.4, label='E')
	#axs[2,0].plot(t_ss, E_ss, 'm-o', linewidth=1)
	#axs[2,0].set_xlabel('t - d')
	#axs[2,0].set_ylabel('E - kg/m3')
	#axs[2,0].set_xlim([1000, t[-1]+10])
	#axs[2,0].set_title('Inactive Biomass Concentration Profile', fontsize=10)
	#axs[2,0].grid(True)
#
	#axs[2,1].plot(t_ss, M_cumulative, 'r-o', linewidth=1.5, label='Cumulative Methane')
	#axs[2,1].set_xlabel('t - d')
	#axs[2,1].set_ylabel('Cumulative Methane - m3/d')
	#axs[2,1].set_title('Cumulative Methane Production', fontsize=10)
	#axs[2,1].set_xlim([1000, t_ss[-1]])
	#axs[2,1].grid(True)
#
	#axs[3,0].plot(t, _M, 'g--', linewidth=0.4, label='M')
	#axs[3,0].plot(t_ss, M_ss, 'g-o', linewidth=1)
	#axs[3,0].set_xlabel('t - d')
	#axs[3,0].set_ylabel('M - kg/m3')
	#axs[3,0].set_xlim([1000, t[-1]+10])
	#axs[3,0].set_title('Methane Concentration Profile', fontsize=10)
	#axs[3,0].grid(True)
#
	#axs[3,1].plot(t, _M_vol, 'r--', linewidth=0.4, label='M')
	#axs[3,1].plot(t_ss, M_vol_ss, 'r-o', linewidth=1)
	#axs[3,1].set_xlabel('t - d')
	#axs[3,1].set_ylabel('M - m3/d')
	#axs[3,1].set_xlim([1000, t[-1]+10])
	#axs[3,1].set_title('Methane Volumetric Production', fontsize=10)
	#axs[3,1].grid(True)
#
plt.show()
