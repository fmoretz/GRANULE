import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

import CSTR_ode as ode

from kinetic_coefficients import *
from physical_coefficients import *
from design_input import *
from influent_input import *

plt.close('all')
N   = 10 		# number of CSTRs

# So =0.22741
# Xo =2
# Eo =0
# Mo =0
# Ro =2.5
# Tw =50
# To =30

# V  = 180
# Dc = 4.05
# H  = 14
# H_ef = 13
# Q = 27
# u = 2.1
# A = 3.1415 * Dc**2 / 4
# Vef = A * H_ef

V    = float(V)
Dc   = float(Dc)
H    = float(H)
H_ef = float(H_ef)
Q    = float(Q)
u    = float(u)

A = float(3.1415 * Dc**2 / 4)

Vef = float(A * H_ef)

m   = rho_biomass * Vef;  			# kg
U   = 400 * 3600 					# J/m2/h/K
A_lateral   = 3.1415*Dc*H/N     	# m2


D   = (7.4 * 10**(-8) *(phi * MM_H2O)**0.5 *To) / (viscosity_H2O * V_H2O**0.6) * 3600 /10**4    #m2/h
_m  = (0.013*(To-273.15)-0.129)/24  # h-1


# sunpharma scenario 1
#day_sim = [10, 20, 50, 100];
day_sim = [365]

plt.figure()
for index in range(0, len(day_sim)):
	day = day_sim[index] 	# days of simulations

	# ∞∞ INIT ∞∞
	npt = 251

	# ∞∞ BCON ∞∞
	Sp  = np.zeros(npt, dtype='float64')

	_S = []; _X = []; _E = []; _T=[]; _Rnew = [];  _eta = []; Nok = []
	_Yv = []; _Vbed = []; _ratio = []; _Ych4 = []; _M = []; _M_vol = []
	_Rnew.append(Ro/1000)

	# ∞∞ CYCL ∞∞
	Np = round((Xo+Eo)/rho_biomass * 3/4 * 1/(3.1415 * (Ro/1000)**3))
	Na = round((Xo)/rho_biomass * 3/4 * 1/(3.1415 * (Ro/1000)**3))

	for k in range(1, N+1):

		dr = _Rnew[-1]/npt
		r = []; _r = 0;

		km  = 2*D/(2*_Rnew[-1])	# m/h - perfect sphere

		for i in range(0, npt):
			_r = _r + dr
			r.append(_r)

		if (Xo + Eo) >= 0.25*rho_biomass:
			print('N:{} passed'.format(k))
			pass
		else:
			print("N:{}, lhs:{}, rhs:{}, r:{}".format(k, Xo+Eo, 0.25*rho_biomass, r[-1]))
			Nok.append(k)
			# Neumann   - Surface of particle
			Sp[-1] = ( km*So - D/dr*Sp[npt-1] )/( km - D/dr )

			for i in range(0, npt):

				if (i == npt-1):
					break
				else:
					etha = (_m/Y/D*Xo) * Sp[i]/(Ks + Sp[i])
					alfa = (2*dr + 2*r[i])/(r[i]*dr**2) * Sp[i]
					beta = -1/dr**2 * Sp[i-1]

					# Solution to Sp
					Sp[i+1] = r[i]*dr**2/(2*dr + r[i]) * (alfa + beta + etha)

					# Dirichlet - Center of particle
					Sp[0]  = Sp[1]

			# efficiency
			eta  = (3 * D * (Sp[-1] - Sp[-2])/dr) / ( r[-1] * (_m*Xo/Y) * Sp[-1]/(Ks + Sp[-1]))
			_eta.append(eta)

			Rnew = ( 3/4 * (Xo+Eo)/(rho_biomass * 3.1415 * Np) )**(1/3)

			_Rnew.append(Rnew)

			Q = km * (So - Sp[0])       #-D * (Sp[0] - Sp[-1])/(2*dr)

			_R = Q * (4*3.1415*r[-1]**2)*Np

			#if Q in globals():
			#	pass
			#else:
			#	Q = u * (3.1415*Dc**2)/4
			#	print(Q)

			HRT = Vef / Q	#h
			B   = B0 *( 1 - K / (_m * HRT/(Kd*HRT + 1) + K-1) ) 	# m3 CH4 / kg COD added
			Yv  = B * So / HRT                                   	# m3 CH4 / kg COD added / h
			Ych4 = Yv * So * Vef									# m3 CH4 / h
			_Ych4.append(Ych4)
			_Yv.append(Yv)

			V_molar = 8.3145 * 10**3 * To / (101325 + rho_biomass*9.8066*H)   				#m3/kmol

			Vload  = 0.25*Vef
			Vbed   = Vload + Yv * MM_CH4 / ( V_molar * Y * _R )
			_Vbed.append(Vbed)
			ratio = Vbed / V
			_ratio.append(ratio)

			t = np.linspace(		# Time span definition
			0,						# Start - h
			day/N*24,				# End   - h
			npt						# Number of iter
			)

			sol = odeint(
				ode.fCSTR,
				[So, Xo, Eo, To, Mo],
				t,
				args=(Q, Vef/N, So, _R, Y, Kd, eta, U, A_lateral, rho_biomass, cp, m, To, Tw, Na),
				atol=1e-7,
				rtol=1e-9,
				mxstep=100000
				)

			S = sol[:,0]
			X = sol[:,1]
			E = sol[:,2]
			T = sol[:,3]
			M = sol[:,4]

			So = S[-1]
			Xo = Xo
			Eo = E[-1]
			To = T[-1]
			Mo = 0

			prod_ch4 = M / rho_ch4 * Vef      #m3 ch4
			M_vol = prod_ch4 / HRT* 24	      #m3/d

			_m  = (0.013*(To-273.15)-0.129)/24  # h-1
			D  = (7.4 * 10**-8 *(phi * MM_H2O)**0.5 *To) / (viscosity_H2O * V_H2O**0.6) * 3600 /10**4    #m2/h


			for j in range(0, len(t)):
				_S.append(S[j])
				_X.append(X[j])
				_E.append(E[j])
				_T.append(T[j])
				_M.append(M[j])
				_M_vol.append(M_vol[j])

	print(_Yv)
	print(_Vbed)
	print(_ratio)

	# ∞∞ SCAL ∞∞
	for i in range(0, npt):
		r[i]  = r[i]  * 1000
		Sp[i] = Sp[i] * 1000

	print(_Rnew)
	for i in range(0, len(Nok)):
		_Rnew[i] = _Rnew[i]*1000

	# ∞∞ VISZ ∞∞
	t = np.linspace(0, day, len(_S))

	#plt.plot(r, Sp, linewidth=2, label='elapsed:{}'.format(day))
	#plt.xlabel('r - mm')
	#plt.ylabel('Sp - g/m3')
	#plt.grid(True)
	#plt.legend()

	fig, axs = plt.subplots(3, 2, figsize=(6, 8))
	plt.subplots_adjust(
		left   = 0.125,
		bottom = 0.071,
		right  = 0.9,
		top    = 0.971,
		wspace = 0.2,
		hspace = 0.287
		)

	axs[0,0].plot(r, Sp, 'k-', linewidth=2)
	axs[0,0].set_xlabel('r - mm')
	axs[0,0].set_ylabel('Sp - g/m3')
	axs[0,0].grid(True)

	axs[0,1].plot(t, _S, 'b-', linewidth=1, label='S')
	axs[0,1].set_xlabel('t - d')
	axs[0,1].set_ylabel('S, - kg/m3')
	axs[0,1].grid(True)

	axs[1,0].plot(t, _X, 'y-', linewidth=2, label='X')
	axs[1,0].set_xlabel('t - d')
	axs[1,0].set_ylabel('X - kg/m3')
	axs[1,0].grid(True)


	axs[1,1].plot(t, _E, 'r-', linewidth=2, label='E')
	axs[1,1].set_xlabel('t - d')
	axs[1,1].set_ylabel('E - kg/m3')
	axs[1,1].grid(True)

	axs[2,0].plot(t, _M, 'g-', linewidth=2, label='M')
	axs[2,0].set_xlabel('t - d')
	axs[2,0].set_ylabel('M - kg/m3')
	axs[2,0].grid(True)

	axs[2,1].plot(t, _M_vol, 'r-', linewidth=2, label='M')
	axs[2,1].set_xlabel('t - d')
	axs[2,1].set_ylabel('M - m3/d')
	axs[2,1].grid(True)




	plt.figure()
	plt.plot(np.linspace(0, Nok[-1], Nok[-1]+1), _Rnew)
	plt.xlabel('N° CSTR')
	plt.ylabel('Particle Radius - mm')
	plt.grid(True)
	plt.title('Increase in Particle Radius along the column')


	plt.figure()
	plt.plot(np.linspace(0, Nok[-1], Nok[-1]-2), _eta[1:-1])
	plt.xlabel('N° CSTR')
	plt.ylabel('Efficiency')
	plt.grid(True)
	plt.title('Internal mass transfer efficiency')

	pic, ass = plt.subplots(2, 1, figsize=(8, 8), gridspec_kw={'height_ratios': [1, 1]})
	plt.subplots_adjust(
		left   = 0.125,
		bottom = 0.071,
		right  = 0.9,
		top    = 0.971,
		wspace = 0.2,
		hspace = 0.287
		)

	ass[0].plot(np.linspace(1, Nok[-1], Nok[-1]-2), _Yv[1:-1],  'k-', linewidth=2)
	ass[0].set_xlabel('N° CSTR ')
	ass[0].set_ylabel('Yv - m3CH4/kgCOD/h')
	ass[0].grid(True)

	ass[1].plot(np.linspace(1, Nok[-1], Nok[-1]-2), _Vbed[1:-1], 'b-', linewidth=2)
	ass[1].set_xlabel('N° CSTR ')
	ass[1].set_ylabel('Vbed - m3')
	ass[1].set_ylim([0, V])
	ass[1].grid(True)

	# Plot Ych4
	plt.figure()
	plt.plot(np.linspace(1, Nok[-1], Nok[-1]-2), _Ych4[1:-1], 'k-', linewidth=2)
	plt.xlabel(xlabel='N° CSTR')
	plt.ylabel(ylabel='Ych4 - m3CH4/kgCOD/h')
	plt.grid(True)

	#Plot methane production M
	plt.figure()
	plt.plot(t, _M, 'k-', linewidth=2)
	plt.xlabel(xlabel='time - d')
	plt.ylabel(ylabel='M - kg/m3')
	plt.grid(True)


# Evaluate area below of the curve _Ych4 with trapezoidal rule
lb = _M[0]		# kg/m3
ub = _M[-1]		# kg/m3
np = len(_M)
step = 0.5
I_trapz = 0.5 * step * (lb + 2*sum(_M[1:np-1]) + ub)
I_trapz_ref = 1
err_traps = 2 - I_trapz / I_trapz_ref
prod_ch4_tot = I_trapz / rho_ch4 * Vef
prod_ch4_daily = prod_ch4_tot / HRT *24

print(f'Total methane produced: {I_trapz} kg/m3')
print(f'volumetric methane daily production: {prod_ch4_daily} m3/d')
print(f'Total volumetric methane production: {prod_ch4_tot} m3')
print(f'Error in trapezoidal rule: {err_traps}')


	# LAST ROW
plt.show()
