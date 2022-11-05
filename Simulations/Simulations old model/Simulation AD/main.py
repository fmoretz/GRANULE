import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, odeint
from CSTR_ode_washout import *
from design_input import *
from influent_input import *
from physical_coefficients import *
from kinetic_coefficients import *

plt.close('all')
N   = 20 		          # number of CSTRs

# ∞∞ DATA ∞∞
m   = rho_biomass * Vef;  		# kg
U   = 400 * 3600 				# J/m2/h/K
A_lateral   = 3.1415*Dc*H/N     # m2

viscosity_H2O = 2.414/100*10**(247.8/(To - 140))    #cP
D   = (7.4 * 10**-8 *(phi * MM_H2O)**0.5 *To) / (viscosity_H2O * V_H2O**0.6) * 3600 /10**4    #m2/h
_m  = (0.013*(To-273.15)-0.129)/24  				# h-1
day_sim = [100]

#whashout weights and limits
alpha = 0.2
beta  = 0.3
gamma = 0.3
delta = 0.2

V_critical = 0.75*Vef   		#m3
OLR_critical = 2  				#m/h
u_critical = 1    				#m/h
M_critical = 5 					#m3/d

#initial volume of sludge bed 


for index in range(0, len(day_sim)):
	day = day_sim[index] 	# days of simulations
	# ∞∞ INIT ∞∞
	npt = 251

	# ∞∞ BCON ∞∞
	Sp  = np.zeros(npt, dtype='float64')


	_S = []; _X = []; _E = []; _T=[]; _Rnew = [];  _eta = []; _Thiele = []
	_Yv = []; _Vbed = []; _ratio = []; _Ych4 = []; _M = []; _M_vol = []
	M_ss = np.zeros(N+1); M_vol_ss = np.zeros(N+1); M_cumulative = np.zeros(N+1)
	S_ss = np.zeros(N+1); X_ss = np.zeros(N+1); E_ss = np.zeros(N+1)
	S_ss[0] = So; X_ss[0] = Xo; E_ss[0] = Eo;
	wM =0;
	_w = [];  _wup = []; _wOLR =[]; _wbed = []; _wM = []
	_Pe = []; _Bi = []
	_omega = []
	_Rnew.append(Ro/1000)

	Vbed = 0.25*Vef
	
	# ∞∞ CYCL ∞∞
	for k in range(1, N+1):

		dr = _Rnew[-1]/npt
		r = []; _r = 0;

		Dp  = _Rnew[-1]*2
		Re  = rho_H2O * u/3600 * Dp /(viscosity_H2O/1000)
		Sc  = (viscosity_H2O/1000) / ( rho_H2O * D/3600)
		Sh  = 2 + ( 1.6*Re**(1/3) + 0.6*Re**(0.5) + 5e-3*Re**(0.8) ) * Sc**(1/3)
		Ar  = (rho_H2O * gravity * (rho_biomass - rho_H2O) * Dp**3) / viscosity_H2O/1000
		
		
		bed_porosity  = ( (18*Re + 2.7*Re**1.69) /Ar)**(1/4.7)
		print(bed_porosity)
		Np = Vbed *(1-bed_porosity)/ (4/3 * 3.14 * _Rnew[-1]**3)
		print(Np)
		
		km  = Sh * D / Dp


		Pe = Re*Sc
		_Pe.append(Pe)
		Bi = km *_Rnew[-1]/D
		_Bi.append(Bi)

		for i in range(0, npt):
			_r = _r + dr
			r.append(_r)

		#Neumann   - Surface of particle
		Sp[-1] = ( km*So - D/dr*Sp[npt-1] )/( km - D/dr )

		for i in range(0, npt):

			if (i == npt-1):
				break
			else:
				g = (_m/Y/D*Xo) * Sp[i]/(Ks + Sp[i])
				a = (2*dr + 2*r[i])/(r[i]*dr**2) * Sp[i]
				b = -1/dr**2 * Sp[i-1]
				# Solution to Sp
				Sp[i+1] = r[i]*dr**2/(2*dr + r[i]) * (a + b + g)
				# Dirichlet - Center of particle
				Sp[0]  = Sp[1]

		# efficiency
		eta  = (3 * D * (Sp[-1] - Sp[-2])/dr) / ( r[-1] * (_m*Xo/Y) * Sp[-1]/(Ks + Sp[-1]))
		_eta.append(eta)
		Rnew = ( 3/4 * (Xo+Eo)/(rho_biomass * 3.1415 * Np) )**(1/3)
		_Rnew.append(Rnew)
		_q = km * (So - Sp[0])       #-D * (Sp[0] - Sp[-1])/(2*dr)

		_R = _q * (4*3.1415*r[-1]**2)*Np
		if Q in locals():
			pass
		else:
			Q = u * (3.1415*Dc**2)/4

		K = 0.6 + 0.0206 * np.exp(0.051*So)
		b = K/Y
		HRT = Vef / Q	#h
		B   = B0 *( 1 - K / (_m * HRT/(Kd*HRT + 1) + K-1) ) 	# m3 CH4 / kg COD added
		Yv  = B * So / HRT                                   	# m3 CH4 / kg COD added / h
		Ych4 = Yv * So * Vef									# m3 CH4 / h
		_Ych4.append(Ych4)
		_Yv.append(Yv)

		V_molar = 8.3145 * 10**3 * To / (101325 + rho_biomass*9.8066*H_ef)   				#m3/kmol
		Vload  = 0.25*Vef
		Vbed   = Vload + Yv * MM_CH4 / ( V_molar * Y * _R )
		_Vbed.append(Vbed)
		ratio = Vbed / V
		_ratio.append(ratio)
	

		if Vbed > V_critical: 
			wbed = 0.001
		else: 
			wbed = 0

		if So*24 > OLR_critical:
			wOLR = 0.001
		else:
			wOLR = 0

		if u >  u_critical: 
			wup = 8*10**(-5) * np.exp(0.205 * u) 
		else:
			wup = 0

		_wOLR.append(wOLR)
		_wM.append(wM)
		_wbed.append(wbed)
		_wup.append(wup)

		w = alpha*wOLR + beta*wup + gamma*wbed + delta*wM
		omega = Q*w
		_w.append(w)
		_omega.append(omega)

		t = np.linspace(	# Time span definition
		0,					# Start - h
		day/N*24,			# End   - h
		npt					# Number of iter
		)

		sol = odeint(
		fCSTR,
		[So, Xo, Eo, To, Mo],
		t,
		args=(Q, Vef/N, So, Xo, Eo, _R, Y, Kd, eta, U, A_lateral, rho_biomass, cp, m, To, Tw, w),
		atol=1e-7,
		rtol=1e-9,
		mxstep=100000
		)

		S = sol[:,0]
		X = sol[:,1]
		E = sol[:,2]
		T = sol[:,3]
		M = sol[:,4]

		S_ss[k] = S[-1]
		X_ss[k] = X[-1]
		E_ss[k] = E[-1]
		M_ss[k] = M[-1]


		So = S[-1]
		Xo = X[-1]
		Eo = E[-1]
		To = T[-1]
		Mo = 0

		prod_ch4 = M / rho_ch4 * Vef     					 #m3 ch4
		M_vol = prod_ch4 / HRT* 24	      					 #m3/d
		M_vol_ss[k] = M_ss[k] / rho_ch4 * Vef / HRT* 24  	 #m3/d
		M_cumulative[k] = M_cumulative[k-1] + M_vol_ss[k]    #m3/d

		for j in range(0, len(M_vol)):
			if M_vol[j] > M_critical:
				wM = 0.001
			else:
				wM = 0
			
	
		_m  = (0.013*(To-273.15)-0.129)/24  # h-1
		viscosity_H2O = 2.414/100*10**(247.8/(To - 140))   #cP
		D  = (7.4 * 10**-8 *(phi * MM_H2O)**0.5 *To) / (viscosity_H2O * V_H2O**0.6) * 3600 /10**4    #m2/h

		Thiele = Rnew / 3 *(_m/Y *Xo /(Ks + Sp[-1])/D )**0.5 
		_Thiele.append(Thiele)
		


		for j in range(0, len(t)):
			_S.append(S[j])
			_X.append(X[j])
			_E.append(E[j])
			_T.append(T[j])
			_M.append(M[j])
			_M_vol.append(M_vol[j])

	r_ext  = np.linspace(r[-1], 1.5*r[-1] , npt)
	Sp_ext = np.zeros(npt)
	Sp_ext = So - (So - Sp[-1])/r_ext * r[-1]

	# ∞∞ SCAL ∞∞
	for i in range(0, npt):
		r[i]  = r[i]  * 1000
		Sp[i] = Sp[i]

	for i in range(0, N +1):
		_Rnew[i] = _Rnew[i]*1000

	# ∞∞ VISZ ∞∞
	t = np.linspace(0, day, len(_S))
	t_ss = np.linspace(0, day, len(M_ss))

	#plt.plot(r, Sp, linewidth=2, label='elapsed:{}'.format(day))
	#plt.xlabel('r - mm')
	#plt.ylabel('Sp - g/m3')
	#plt.grid(True)
	#plt.legend()



	fig, axs = plt.subplots(4, 2, figsize=(6, 8))
	plt.subplots_adjust(
		left   = 0.125,
		bottom = 0.071,
		right  = 0.9,
		top    = 0.971,
		wspace = 0.2,
		hspace = 0.5
		)

	Sp_graph = np.array(Sp) * Xinitial
	Sp_ext_graph = np.array(Sp_ext) * Xinitial 
	_S_graph = np.array(_S) * Xinitial
	S_ss_graph = np.array(S_ss) * Xinitial
	_X_graph = np.array(_X) * Xinitial
	X_ss_graph = np.array(X_ss) * Xinitial
	_E_graph = np.array(_E) * Xinitial
	E_ss_graph = np.array(E_ss) * Xinitial
	M_cumulative_graph = np.array(M_cumulative) * Xinitial
	_M_graph = np.array(_M) * Xinitial	
	M_ss_graph = np.array(M_ss) * Xinitial
	_M_vol_graph = np.array(_M_vol)*Xinitial
	M_vol_ss_graph = np.array(M_vol_ss)*Xinitial


	axs[0,0].plot(r, Sp_graph, 'k-o', linewidth=1.5, label='Sp', markevery=41)
	axs[0,0].plot(r_ext*1000, Sp_ext_graph, 'r-o', linewidth=1.5, label='Sp', markevery=61)
	axs[0,0].set_xlabel('r - mm')
	axs[0,0].set_ylabel('Sp - kg/m3')
	axs[0,0].set_title('Substrate gradient inside the granule', fontsize=10)
	axs[0,0].grid(True)

	axs[0,1].plot(r_ext*1000, Sp_ext_graph, 'k-o', linewidth=1.5, label='Sp', markevery=61)
	axs[0,1].set_xlabel('r - mm')
	axs[0,1].set_ylabel('Sp - kg/m3')
	axs[0,1].grid(True)
	axs[0,1].set_title('Substrate profile outside particle', fontsize=10)

	axs[1,0].plot(t, _S_graph, 'b--', linewidth=0.4, label='S')
	axs[1,0].plot(t_ss, S_ss_graph, 'b-o', linewidth=1)
	axs[1,0].set_xlabel('t - d')
	axs[1,0].set_ylabel('S, - kg/m3')
	axs[1,0].set_xlim([0, t[-1]+10])
	axs[1,0].set_title('Substrate Concentration Profile', fontsize=10)
	axs[1,0].grid(True)

	axs[1,1].plot(t, _X_graph, 'y--', linewidth=0.4, label='X')
	axs[1,1].plot(t_ss, X_ss_graph, 'y-o', linewidth=1)
	axs[1,1].set_xlabel('t - d')
	axs[1,1].set_ylabel('X - kg/m3')
	axs[1,1].set_xlim([0, t[-1]+10])
	axs[1,1].set_ylim([0, 1.5*max(_X_graph)])
	axs[1,1].set_title('Active Biomass Concentration Profile', fontsize=10)
	axs[1,1].grid(True)


	axs[2,0].plot(t, _E_graph, 'm--', linewidth=0.4, label='E')
	axs[2,0].plot(t_ss, E_ss_graph, 'm-o', linewidth=1)
	axs[2,0].set_xlabel('t - d')
	axs[2,0].set_ylabel('E - kg/m3')
	axs[2,0].set_xlim([0, t[-1]+10])
	axs[2,0].set_title('Inactive Biomass Concentration Profile', fontsize=10)
	axs[2,0].grid(True)

	axs[2,1].plot(t_ss, M_cumulative_graph, 'r-o', linewidth=1.5, label='Cumulative Methane')
	axs[2,1].set_xlabel('t - d')
	axs[2,1].set_ylabel('Cumulative Methane - m3/d')
	axs[2,1].set_title('Cumulative Methane Production', fontsize=10)
	axs[2,1].set_xlim([0, t_ss[-1]])
	axs[2,1].grid(True)

	axs[3,0].plot(t, _M_graph, 'g--', linewidth=0.4, label='M')
	axs[3,0].plot(t_ss, M_ss_graph, 'g-o', linewidth=1)
	axs[3,0].set_xlabel('t - d')
	axs[3,0].set_ylabel('M - kg/m3')
	axs[3,0].set_xlim([0, t[-1]+10])
	axs[3,0].set_title('Methane Concentration Profile', fontsize=10)
	axs[3,0].grid(True)

	axs[3,1].plot(t, _M_vol_graph, 'r--', linewidth=0.4, label='M')
	axs[3,1].plot(t_ss, M_vol_ss_graph, 'r-o', linewidth=1)
	axs[3,1].set_xlabel('t - d')
	axs[3,1].set_ylabel('M - m3/d')
	axs[3,1].set_xlim([0, t[-1]+10])
	axs[3,1].set_title('Methane Volumetric Production', fontsize=10)
	axs[3,1].grid(True)


	#plt.plot(t_ss, M_cumulative, 'r-o', linewidth=1.5, label='Cumulative Methane')
	#plt.xlabel('t - d')
	#plt.ylabel('Cumulative Methane - m3/d')
	#plt.xlim([0, t_ss[-1]])
	#plt.grid(True)
	#plt.title('Cumulative Methane Production - m3/d', fontsize=10)

	fig2, axs2 =plt.subplots(3, 2, figsize=(6, 8))
	plt.subplots_adjust(
		left   = 0.125,
	    bottom = 0.071,
	    right  = 0.9,
	    top    = 0.971,
	    wspace = 0.2,
	    hspace = 0.5
		)

	axs2[0,0].plot(np.linspace(0, N, N+1), _Rnew, 'k-o', linewidth=1.5)
	axs2[0,0].set_xlabel('N° CSTR')
	axs2[0,0].set_ylabel('Particle Radius - mm')
	axs2[0,0].set_xlim([0,N])
	axs2[0,0].set_title('Increase of the Particle Radius', fontsize=10)
	axs2[0,0].grid(True)

	axs2[0,1].plot(np.linspace(1, N, N-2), _eta[1:-1], 'b-o', linewidth=1, label='S')
	axs2[0,1].set_xlabel('N° CSTR')
	axs2[0,1].set_ylabel('Efficiency')
	axs2[0,1].set_xlim([0,N])
	axs2[0,1].set_title('Internal Mass Transfer Efficiency', fontsize=10)
	axs2[0,1].grid(True)


	axs2[1,0].plot(np.linspace(1, N, N-2), _Yv[1:-1],  'k-o', linewidth=2)
	axs2[1,0].set_xlabel('N° CSTR ')
	axs2[1,0].set_ylabel('Ych4 - m3CH4/kgCOD/h')
	axs2[1,0].set_xlim([0,N])
	axs2[1,0].set_title('Methane Production Rate', fontsize=10)
	axs2[1,0].grid(True)

	axs2[2,1].plot(np.linspace(1, N, N-2), _Vbed[1:-1], 'b-o', linewidth=2)
	axs2[2,1].set_xlabel('N° CSTR ')
	axs2[2,1].set_ylabel('Vbed - m3')
	axs2[2,1].set_xlim([0,N])
	axs2[2,1].set_ylim([0, V])
	axs2[2,1].set_title('Bed Volume Expansion', fontsize=10)
	axs2[2,1].grid(True)

	axs2[2,0].plot(np.linspace(1, N, N-2), _Ych4[1:-1], 'k-o', linewidth=2)
	axs2[2,0].set_xlabel('N° CSTR')
	axs2[2,0].set_ylabel('Ych4 - m3CH4/kgCOD/h')
	axs2[2,0].set_xlim([0,N])
	axs2[2,0].set_title('Methane Production Rate', fontsize=10) #to change
	axs2[2,0].grid(True)

	axs2[1,1].plot(np.linspace(1, N, N-2), _Thiele[1:-1], 'b-o', linewidth=1, label='S')
	axs2[1,1].set_xlabel('N° CSTR')
	axs2[1,1].set_ylabel('Thiele Modulus')
	axs2[1,1].set_xlim([0,N])
	axs2[1,1].set_title('Thiele Modulus', fontsize=10)
	axs2[1,1].grid(True)

	fig3, axs3 =plt.subplots(3, 2, figsize=(6, 8))
	plt.subplots_adjust(
		left   = 0.125,
	    bottom = 0.071,
	    right  = 0.9,
	    top    = 0.971,
	    wspace = 0.2,
	    hspace = 0.5
		)

	axs3[0,0].plot(np.linspace(0, N, N), _wOLR, 'b-o', linewidth=1, label='washout OLR')
	axs3[0,0].set_xlabel('N° CSTR')
	axs3[0,0].set_ylabel('OLR washout function')
	axs3[0,0].set_xlim([0,N])
	axs3[0,0].set_title('OLR washout function', fontsize=10)
	axs3[0,0].grid(True)

	axs3[0,1].plot(np.linspace(0, N, N), _wbed, 'b-o', linewidth=1, label='Washout bed')
	axs3[0,1].set_xlabel('N° CSTR')
	axs3[0,1].set_ylabel('Bed washout function')
	axs3[0,1].set_xlim([0,N])
	axs3[0,1].set_title('Bed washout function', fontsize=10)
	axs3[0,1].grid(True)

	axs3[1,0].plot(np.linspace(0, N, N), _wup, 'b-o', linewidth=1, label='Washout upflow velocity')
	axs3[1,0].set_xlabel('N° CSTR')
	axs3[1,0].set_ylabel('Upflow velocity washout')
	axs3[1,0].set_xlim([0,N])
	axs3[1,0].set_title('Upflow velocity washout', fontsize=10)
	axs3[1,0].grid(True)
 
	axs3[1,1].plot(np.linspace(0, N, N), _wM, 'b-o', linewidth=1, label='Washout methane')
	axs3[1,1].set_xlabel('N° CSTR')
	axs3[1,1].set_ylabel('Methane washout')
	axs3[1,1].set_xlim([0,N])
	axs3[1,1].set_title('Methane washout', fontsize=10)
	axs3[1,1].grid(True)

	axs3[2,0].plot(np.linspace(0, N, N), _w, 'b-o', linewidth=1, label='Washout')
	axs3[2,0].set_xlabel('N° CSTR')
	axs3[2,0].set_ylabel('Washout function')
	axs3[2,0].set_xlim([0,N])
	axs3[2,0].set_title('Washout function', fontsize=10)
	axs3[2,0].grid(True)

	axs3[2,1].plot(np.linspace(0, N, N), _omega, 'b-o', linewidth=1, label='omega')
	axs3[2,1].set_xlabel('N° CSTR')
	axs3[2,1].set_ylabel('Omega [m3/h]')
	axs3[2,1].set_xlim([0,N])
	axs3[2,1].set_title('Washout flow [m3/h]', fontsize=10)
	axs3[2,1].grid(True)

	fig4, axs4 =plt.subplots(2, 2, figsize=(6, 8))
	plt.subplots_adjust(
		left   = 0.125,
	    bottom = 0.071,
	    right  = 0.9,
	    top    = 0.971,
	    wspace = 0.2,
	    hspace = 0.5
		)

	axs4[0,0].plot(np.linspace(1, N, N-2), _eta[1:-1], 'b-o', linewidth=1, label='eta')
	axs4[0,0].set_xlabel('N° CSTR')
	axs4[0,0].set_ylabel('Efficiency')
	axs4[0,0].set_xlim([0,N])
	axs4[0,0].set_title('Internal Mass Transfer Efficiency', fontsize=10)
	axs4[0,0].grid(True)

	axs4[1,0].plot(np.linspace(1, N, N-2), _Thiele[1:-1], 'b-o', linewidth=1, label='Thiele')
	axs4[1,0].set_xlabel('N° CSTR')
	axs4[1,0].set_ylabel('Thiele Modulus')
	axs4[1,0].set_xlim([0,N])
	axs4[1,0].set_title('Thiele Modulus', fontsize=10)
	axs4[1,0].grid(True)

	axs4[0,1].plot(np.linspace(0, N, N), _Bi, 'b-o', linewidth=1, label='Biot')
	axs4[0,1].set_xlabel('N° CSTR')
	axs4[0,1].set_ylabel('Biot number')
	axs4[0,1].set_xlim([0,N])
	axs4[0,1].set_title('Biot number', fontsize=10)
	axs4[0,1].grid(True)

	axs4[1,1].plot(np.linspace(0, N, N), _Pe, 'b-o', linewidth=1, label='Peclet')
	axs4[1,1].set_xlabel('N° CSTR')
	axs4[1,1].set_ylabel('Peclet number')
	axs4[1,1].set_xlim([0,N])
	axs4[1,1].set_title('Peclet number', fontsize=10)
	axs4[1,1].grid(True)

	x = np.linspace(0, r[-1], 100)
	y = (r[-1]**2 - x**2)**0.5

	color = 'white'
	x_pos1 = 1.0
	x_pos2 = 3
	y_pos = 0.99*Sp_ext_graph[-1]
	circle = plt.Circle((0,0), r[-1]+0.05, color = 'orange', alpha=0.2)
	fig5, axs5 = plt.subplots()
	axs5.plot(r, Sp_graph, 'k', alpha=0.6, linewidth=2, label='Sp', markevery=41)
	axs5.plot(r_ext*1000, Sp_ext_graph, 'k', alpha=0.6, linewidth=2, label='Sp', markevery=61)
	axs5.set_xlabel('x - mm')
	axs5.set_ylabel('Sp - kg/m3')
	axs5.tick_params(axis='y')
	axs6 = axs5.twinx()
	axs6.add_patch(circle)
	axs6.tick_params(axis='y', labelcolor=color, color=color)
	axs6.set_ylim([0, 2.7])
	axs5.set_ylim([Sp_graph[0], Sp_ext_graph[-1]])
	axs5.set_xlim([0, r_ext[-1]*1000])
	axs5.text(x_pos1, y_pos, "Inner Particle")
	axs5.text(x_pos2, y_pos, "Boundary Layer")
	fig5.tight_layout()





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

print(f'volumetric methane daily production: {prod_ch4_daily} m3/d')
print(f'Total volumetric methane production: {prod_ch4_tot} m3')
print(f'Error in trapezoidal rule: {err_traps}')

	# LAST ROW
plt.show()
