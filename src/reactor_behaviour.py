import numpy as np
import pandas as pd
import time
import openpyxl
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import NullFormatter,MultipleLocator, FormatStrFormatter
from matplotlib.ticker import MaxNLocator
from matplotlib.font_manager import FontProperties
import matplotlib.ticker as ticker

date = time.localtime()
df = pd.read_excel("../input/input_data.xlsx", sheet_name='design input', header=None)

up  = df[4][7]
N   = df[4][10]

N_range = np.linspace(1, 100)
Di  = 1- (1.03*up**(1.11)*0.009**(1/N))/0.009
_Di_curve = []
if Di > 0.75:
    reactor = 'CSTR'
elif Di < 0.25:
    reactor = 'PFR'
else:
    reactor = 'HYBRID APPROACH'

for i in range(0, len(N_range)):
    Di_curve = 1- (1.03*up**(1.11)*0.009**(1/N_range[i]))/(1.03*up**(1.11))
    _Di_curve.append(Di_curve)

plt.figure(num=None, dpi=80, facecolor='w', edgecolor='k')
matplotlib.rcParams.update({'font.size': 24, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})
matplotlib.rcParams['font.sans-serif'] = "Arial"

figure, axes = plt.subplots(1, 1, figsize=(8, 8))
axes.plot(N_range, _Di_curve, 'k', linewidth=2.5)
axes.bar(N, Di, width=2.5, color='red', edgecolor='black')
axes.grid(color='k', alpha=0.5, linestyle='dashed', linewidth=0.8)
axes.set_xlabel('N° CSTR')
axes.set_ylabel('Discretization number - [-]')
axes.set_xlim([0, N_range[-1]])
axes.set_ylim([0, 1])
title = f'Reactor model behavior: {reactor}'
axes.set_title(title)
fig_pathway = '../output/Figures/'
figure.savefig(fig_pathway + 'Reactor model behaviour_'+str(date.tm_mday)+'_'+str(date.tm_mon)+'_'+str(date.tm_year)+'.png')



N = int(N)
