import numpy as np 
import math
from scipy.optimize import fsolve
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import NullFormatter,MultipleLocator, FormatStrFormatter
from matplotlib.ticker import MaxNLocator
from matplotlib.font_manager import FontProperties
import matplotlib.ticker as ticker


def system(x):
    a, b = x

    f = np.empty(2)

    f[0] = 0.0001- a*np.exp(1*b)
    f[1] = 0.005 - a*np.exp(20*b)

    return f 

x0 = [0.0001, 0.1]
solution = fsolve(system, x0, xtol =1e-10)
a, b = solution 
print(a)
print(b)

u = np.linspace(0, 20, 1000)
wup=np.zeros(len(u))

def upflow_washout(u):
    if u < 1:
        wup = 0
    else:
        wup = a*np.exp(b*u)
    return wup

for i in range(0,len(u)):
    wup[i] = upflow_washout(u[i])

plt.figure(num=None, dpi=80, facecolor='w', edgecolor='k', figsize=(10, 10))
matplotlib.rcParams.update({'font.size': 12, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})
matplotlib.rcParams['font.sans-serif'] = "Arial"

plt.figure(num=None, dpi=80, facecolor='w', edgecolor='k')
matplotlib.rcParams.update({'font.size': 20, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})
matplotlib.rcParams['font.sans-serif'] = "Arial"

figure, axes = plt.subplots(1, 1, figsize=(10, 10))
axes.plot(u, wup, 'k', linewidth=2.5)
axes.grid(color='k', alpha=0.5, linestyle='dashed', linewidth=0.8)
axes.set_xlabel('u - [m/h]', fontsize=24)
axes.set_ylabel('w$_{up}$ - [-]', fontsize=24)
plt.xlim([0, u[-1]])
plt.ylim([-10**-4, wup[-1]])
fig_pathway = '../output'
figure.savefig(fig_pathway + 'upflow_washout.svg')
