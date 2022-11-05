import numpy as np
import matplotlib.pyplot as plt

from main_sunpharma_500 import _S_graph as S_500, S_ss_graph as S_ss_500,  _X_graph as X_500, X_ss_graph as X_ss_500,  _E_graph as E_500, E_ss_graph as E_ss_500   
from main_sunpharma_500 import M_cumulative_graph as M_cumulative_500, _M_graph as M_500, M_ss_graph as M_ss_500,  _M_vol_graph as M_vol_500, M_vol_ss_graph as M_vol_ss_500   
from main_sunpharma_500 import t as t_500, t_ss as t_ss_500

from main_sunpharma_1000 import _S_graph as S_1000, S_ss_graph as S_ss_1000,  _X_graph as X_1000, X_ss_graph as X_ss_1000,  _E_graph as E_1000, E_ss_graph as E_ss_1000   
from main_sunpharma_1000 import M_cumulative_graph as M_cumulative_1000, _M_graph as M_1000, M_ss_graph as M_ss_1000,  _M_vol_graph as M_vol_1000, M_vol_ss_graph as M_vol_ss_1000   
from main_sunpharma_1000 import t as t_1000, t_ss as t_ss_1000

#from sunpharma_1500 import _S as S_1500, S_ss as S_ss_1500,  _X as X_1500, X_ss as X_ss_1500,  _E as E_1500, E_ss as E_ss_1500   
#from sunpharma_1500 import M_cumulative as M_cumulative_1500, _M as M_1500, M_ss as M_ss_1500,  _M_vol as M_vol_1500, M_vol_ss as M_vol_ss_1500   
#from sunpharma_1500 import t as t_1500, t_ss as t_ss_1500
#
#from sunpharma_2000 import _S as S_2000, S_ss as S_ss_2000,  _X as X_2000, X_ss as X_ss_2000,  _E as E_2000, E_ss as E_ss_2000   
#from sunpharma_2000 import M_cumulative as M_cumulative_2000, _M as M_2000, M_ss as M_ss_2000,  _M_vol as M_vol_2000, M_vol_ss as M_vol_ss_2000   
#from sunpharma_2000 import t as t_2000, t_ss as t_ss_2000


fig, axs = plt.subplots(3, 2, figsize=(6, 8))
plt.subplots_adjust(
	left   = 0.125,
	bottom = 0.071,
	right  = 0.9,
	top    = 0.971,
	wspace = 0.2,
	hspace = 0.5
	)

axs[0,0].plot(t_500,  S_500,  'b--', linewidth=0.4, label='S 500')
axs[0,0].plot(t_1000, S_1000, 'b--', linewidth=0.4, label='S 1000')
axs[0,0].plot(t_ss_500,  S_ss_500,  'b-o', linewidth=1)
axs[0,0].plot(t_ss_1000, S_ss_1000, 'b-o', linewidth=1)
axs[0,0].set_xlabel('t - d')
axs[0,0].set_ylabel('S, - kg/m3')
axs[0,0].set_xlim([0, t_1000[-1]+10])
axs[0,0].set_title('Substrate Concentration Profile', fontsize=10)
axs[0,0].grid(True)

axs[0,1].plot(t_500, X_500, 'y--', linewidth=0.4, label='X 500')
axs[0,1].plot(t_1000, X_1000, 'y--', linewidth=0.4, label='X 1000')
axs[0,1].plot(t_ss_500, X_ss_500, 'y-o', linewidth=1)
axs[0,1].plot(t_ss_1000, X_ss_1000, 'y-o', linewidth=1)
axs[0,1].set_xlabel('t - d')
axs[0,1].set_ylabel('X - kg/m3')
axs[0,1].set_xlim([0, t_1000[-1]+10])
axs[0,1].set_ylim([0, 1.5*max(max(X_500), max(X_1000))])
axs[0,1].set_title('Active Biomass Concentration Profile', fontsize=10)
axs[0,1].grid(True)

E_1000[0] = E_500[-1]
E_ss_1000[0] = E_ss_500[-1]


axs[1,0].plot(t_500, E_500, 'm--', linewidth=0.4, label='E 500')
axs[1,0].plot(t_1000, E_1000, 'm--', linewidth=0.4, label='E 500')
axs[1,0].plot(t_ss_500, E_ss_500, 'm-o', linewidth=1)
axs[1,0].plot(t_ss_1000, E_ss_1000, 'm-o', linewidth=1)
axs[1,0].set_xlabel('t - d')
axs[1,0].set_ylabel('E - kg/m3')
axs[1,0].set_xlim([0, t_1000[-1]+10])
axs[1,0].set_title('Inactive Biomass Concentration Profile', fontsize=10)
axs[1,0].grid(True)

axs[1,1].plot(t_500, M_500, 'g--', linewidth=0.4, label='M 500')
axs[1,1].plot(t_1000, M_1000, 'g--', linewidth=0.4, label='M 1000')
axs[1,1].plot(t_ss_500, M_ss_500, 'g-o', linewidth=1)
axs[1,1].plot(t_ss_1000, M_ss_1000, 'g-o', linewidth=1)
axs[1,1].set_xlabel('t - d')
axs[1,1].set_ylabel('M - kg/m3')
axs[1,1].set_xlim([0, t_1000[-1]+10])
axs[1,1].set_title('Methane Concentration Profile', fontsize=10)
axs[1,1].grid(True)

axs[2,0].plot(t_500, M_vol_500, 'r--', linewidth=0.4, label='M 500')
axs[2,0].plot(t_1000, M_vol_1000, 'r--', linewidth=0.4, label='M 1000')
axs[2,0].plot(t_ss_500, M_vol_ss_500, 'r-o', linewidth=1)
axs[2,0].plot(t_ss_1000, M_vol_ss_1000, 'r-o', linewidth=1)
axs[2,0].set_xlabel('t - d')
axs[2,0].set_ylabel('M - m3/d')
axs[2,0].set_xlim([0, t_1000[-1]+10])
axs[2,0].set_title('Methane Volumetric Production', fontsize=10)
axs[2,0].grid(True)

axs[2,1].plot(t_ss_500, M_cumulative_500, 'r-o', linewidth=1)
axs[2,1].plot(t_ss_1000, M_cumulative_1000, 'r-o', linewidth=1)
axs[2,1].set_xlabel('t - d')
axs[2,1].set_ylabel('M - m3/d')
axs[2,1].set_xlim([0, t_ss_1000[-1]])
axs[2,1].set_title('Cumulative Methane - m3/d', fontsize=10)
axs[2,1].grid(True)

plt.show()
