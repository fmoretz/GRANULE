import numpy as np 
import math
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

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

plt.figure()
plt.plot(u, wup, 'r-o', linewidth=1.5, label='upflow washout', markevery=20)
plt.grid(True)
plt.xlabel('Upflow velocity [m/h]')
plt.ylabel('Upflow washout')
plt.xlim([0, u[-1]])
plt.ylim([-10**-4, wup[-1]])
plt.title('Upflow Washout as a function of the upflow velocity', fontsize=10)
plt.legend(loc='best')
plt.show()