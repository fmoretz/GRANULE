import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import NullFormatter,MultipleLocator, FormatStrFormatter
from matplotlib.ticker import MaxNLocator
from matplotlib.font_manager import FontProperties
import matplotlib.ticker as ticker

from kinetic_coefficients import Y, Ks, B0, Kd
from influent_input import  To
		
        
_m  = (0.013*(To-273.15)-0.129)/24  # h-1
S = np.linspace(0, 50, 1000)

print(_m)
print(Ks)

_mu = []

for i in range(0, len(S)):
    mu = _m * S[i]/ (Ks + S[i])
    _mu.append(mu)

plt.figure(num=None, dpi=80, facecolor='w', edgecolor='k')
matplotlib.rcParams.update({'font.size': 20, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})
matplotlib.rcParams['font.sans-serif'] = "Arial"

figure, axes = plt.subplots(1, 1, figsize=(10, 10))
axes.plot(S, _mu, 'k', linewidth=2.5)
axes.grid(color='k', alpha=0.5, linestyle='dashed', linewidth=0.8)
axes.set_xlabel('S$_P$ - [kg/m$^3$]', fontsize=24)
axes.set_ylabel('$\mu$ - [h$^{-1}$]', fontsize=24)
axes.set_xlim([0, S[-1]])
axes.set_ylim([0, 1.05*max(_mu)])
title = f'Microbial growth rate'
axes.set_title(title)
fig_pathway = '../output/Figures/'
figure.savefig(fig_pathway + 'Microbial growth rate.svg')
