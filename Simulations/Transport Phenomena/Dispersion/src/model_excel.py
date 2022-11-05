from CSTR_ode_washout import* 
from reactor_behaviour import*
def model():
	
	import numpy as np
	import time 
	start = time.time()
	import matplotlib
	import xlwt
	from xlwt import Workbook
	import xlsxwriter
	import matplotlib.pyplot as plt
	from matplotlib.colors import Normalize
	from matplotlib.ticker import NullFormatter,MultipleLocator, FormatStrFormatter
	from matplotlib.ticker import MaxNLocator
	from matplotlib.font_manager import FontProperties
	import matplotlib.ticker as ticker
	from scipy.integrate import odeint
	from design_input import V, Dc, H, H_ef, Q, u, A, Vef
	from influent_input import Sinitial, Xinitial, Eo, Mo, Ro, Tw, To, day, So, Xo
	from physical_coefficients import rho_biomass, cp, phi, rho_H2O, MM_H2O, V_H2O, MM_CH4, rho_ch4 
	from kinetic_coefficients import Y, Ks, B0, Kd

 		          # number of CSTRs

	# ∞∞ DATA ∞∞
	m   = rho_biomass * Vef;  		# kg
	U   = 400 * 3600 				# J/m2/h/K
	A_lateral   = 3.1415*Dc*H/N     # m2

	_D = []; _km = [] 

	viscosity_H2O = 2.414/100*10**(247.8/(To - 140))    #cP
	D   = (7.4 * 10**-8 *(phi * MM_H2O)**0.5 *To) / (viscosity_H2O * V_H2O**0.6) * 3600 /10**4    #m2/h
	_D.append(D)
	_m  = (0.013*(To-273.15)-0.129)/24  				# h-1

	#whashout weights and limits
	alpha = 0.2
	beta  = 0.3
	gamma = 0.3
	delta = 0.2

	V_critical = 0.75*Vef   		#m3
	OLR_critical = 2  				#m/h
	u_critical = 1    				#m/h
	M_critical = 5 					#m3/d

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
	_D_TIS = []; _Pe_TIS = []; _Dispersion_coefficient = []
	_Rnew.append(Ro/1000)

	# ∞∞ CYCL ∞∞
	Np = round((Xo+Eo)/rho_biomass * 3/4 * 1/(3.1415 * (Ro/1000)**3))
	Na = round((Xo)/rho_biomass * 3/4 * 1/(3.1415 * (Ro/1000)**3))

	for k in range(1, N+1):

		dr = _Rnew[-1]/npt
		r = []; _r = 0;

		Dp  = _Rnew[-1]*2
		Re  = rho_H2O * u/3600 * Dp /(viscosity_H2O/1000)
		Sc  = (viscosity_H2O/1000) / ( rho_H2O * D/3600)
		Sh  = 2 + ( 1.6*Re**(1/3) + 0.6*Re**(0.5) + 5e-3*Re**(0.8) ) * Sc**(1/3)

		km  = Sh * D / Dp
		_km.append(km)

		Pe = Re*Sc
		_Pe.append(Pe)
		Bi = km *_Rnew[-1]/D
		_Bi.append(Bi)

		Dispersion_coefficient  = 1.03 *u**(1.11)*0.009**(k/N)
		_Dispersion_coefficient.append(Dispersion_coefficient)
	
		Pe_TIS = 2*(N-1)
		_Pe_TIS.append(Pe_TIS)
		D_TIS  = u*H_ef/Pe_TIS
		_D_TIS.append(D_TIS)


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
		mxstep=100000,
		hmax=0.1
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
		_D.append(D)

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


	# Set up the features of plot 
	plt.figure(num=None, dpi=800, facecolor='w', edgecolor='k')
	matplotlib.rcParams.update({'font.size': 18, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})
	matplotlib.rcParams['font.sans-serif'] = "Arial"

	fig, axs = plt.subplots(2, 4, figsize=(25, 10))
	plt.subplots_adjust(
		left   = 0.125,
		bottom = 0.071,
		right  = 0.9,
		top    = 0.95, 
		wspace = 0.45,
		hspace = 0.75
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
	M_production_average = np.mean(_M_vol_graph) 

	axs[0,0].plot(t_ss, S_ss_graph, 'k', linewidth=2, marker='.', markersize=10, linestyle='solid', label='S')
	axs[0,0].set_xlabel('t - [d]')
	axs[0,0].set_ylabel('S - [kg/m$^3$]')
	axs[0,0].grid(color='k', alpha=0.5, linestyle='dashed', linewidth=0.8)
	axs[0,0].set_xlim([0, t[-1]])
	axs[0,0].set_title('Substrate Concentration Profile', y=1.025, fontsize=18)

	axs[0,1].plot(t, _X_graph, 'k', linewidth=2, marker='.', markersize=10, markevery=round(len(t)/20), linestyle='solid', label='X')
	axs[0,1].set_xlabel('t - [d]')
	axs[0,1].set_ylabel('X - [kg/m$^3$]')
	axs[0,1].set_xlim([0, t[-1]])
	axs[0,1].set_ylim([0.8*min(_X_graph), 1.2*max(_X_graph)])
	axs[0,1].grid(color='k', alpha=0.5, linestyle='dashed', linewidth=0.8)
	axs[0,1].set_title('Active Biomass Concentration Profile',  y=1.025, fontsize=18)

	axs[0,2].plot(t, _E_graph, 'k', linewidth=2, marker='.', markersize=10, markevery=round(len(t)/20), linestyle='solid', label='E')
	axs[0,2].set_xlabel('t - [d]')
	axs[0,2].set_ylabel('E - [kg/m$^3$]')
	axs[0,2].set_xlim([0, t[-1]])
	axs[0,2].grid(color='k', alpha=0.5, linestyle='dashed', linewidth=0.8)
	axs[0,2].set_title('Inactive Biomass Concentration Profile',  y=1.025, fontsize=18)

	axs[0,3].plot(t_ss, M_vol_ss_graph, 'k', linewidth=2, marker='.', markersize=10, linestyle='solid', label='M')
	axs[0,3].set_xlabel('t - [d]')
	axs[0,3].set_ylabel('M - [m$^3$/d]')
	axs[0,3].set_ylim([0, 1.05*max(_M_vol_graph)])
	axs[0,3].set_xlim([0, t[-1]])
	axs[0,3].grid(color='k', alpha=0.5, linestyle='dashed', linewidth=0.8)
	axs[0,3].set_title('Methane Volumetric Production',  y=1.025, fontsize=18)

	axs[1,0].plot(t_ss, M_cumulative_graph, 'k', linewidth=2, marker='.', markersize=10, linestyle='solid', label='Cumulative Methane')
	axs[1,0].set_xlabel('t - [d]')
	axs[1,0].set_ylabel('Cumulative Methane - [m$^3$/d]')
	axs[1,0].set_title('Cumulative Methane Production', y=1.025, fontsize=18)
	axs[1,0].set_xlim([0, t_ss[-1]])
	axs[1,0].grid(color='k', alpha=0.5, linestyle='dashed', linewidth=0.8)

	axs[1,1].plot(np.linspace(1, N, N), _eta, 'k', linewidth=2, marker='.', markersize=10, linestyle='solid', label='eta')
	axs[1,1].set_xlabel('N° CSTR')
	axs[1,1].set_ylabel('Efficiency - [-]')
	axs[1,1].set_xlim([1,N])
	axs[1,1].set_title('Internal Mass Transfer Efficiency', y=1.025, fontsize=18)
	axs[1,1].grid(color='k', alpha=0.5, linestyle='dashed', linewidth=0.8)

	axs[1,2].plot(np.linspace(1, N, N), _km, 'k', linewidth=2, marker='.', markersize=10, linestyle='solid', label='eta')
	axs[1,2].set_xlabel('N° CSTR')
	axs[1,2].set_ylabel('km - [m/h]')
	axs[1,2].set_xlim([1,N-1])
	axs[1,2].set_title('External Mass Transfer Coefficient', y=1.025, fontsize=18)
	axs[1,2].grid(color='k', alpha=0.5, linestyle='dashed', linewidth=0.8)

	axs[1,3].plot(np.linspace(1, N, N), _D[1:],  'k', linewidth=2, marker='.', markersize=10, linestyle='solid', label='eta')
	axs[1,3].set_xlabel('N° CSTR')
	axs[1,3].set_ylabel('D - [m$^2$/h]')
	axs[1,3].set_xlim([1,N])
	axs[1,3].set_title('Diffusivity coefficient', y=1.025, fontsize=18)
	axs[1,3].grid(color='k', alpha=0.5, linestyle='dashed', linewidth=0.8)

	fig_pathway = '../output/Figures/'
	fig.savefig(fig_pathway + 'Concentration profiles_'+str(date.tm_mday)+'_'+str(date.tm_mon)+'_'+str(date.tm_year)+'.png')	

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

	removal_efficiency =(Sinitial - _S_graph[-1])/Sinitial*100
	radius_growth = (_Rnew[-1]- Ro)/Ro*100


	end = time.time()
	wb = xlsxwriter.Workbook('../output/simulation_'+str(date.tm_mday)+'_'+str(date.tm_mon)+'_'+str(date.tm_year)+'.xlsx')
	bold = wb.add_format({'bold': True})
	number_format = wb.add_format({'num_format':'0.00'})

	sheet1 = wb.add_worksheet('Numerical results')
	sheet1.write(0,0, f'RESULTS AFTER {round(day)} DAYS OF SIMULATION:', bold)
	sheet1.write(2,0, f'Efficiency of substrate removal:')
	sheet1.write_number(2,4, round(removal_efficiency, 2))
	sheet1.write(2,5, f'%')
	sheet1.write(3,0, f'Final Substrate concentration:')
	sheet1.write_number(3,4, round(_S_graph[-1], 2))  
	sheet1.write(3,5, f'[kg/m3]')          
	sheet1.write(4,0, f'Final Biomass concentration:') 
	sheet1.write_number(4,4, round(_X_graph[-1], 2))	
	sheet1.write(4,5, f'[kg/m3]')            
	sheet1.write(5,0, f'Final Methane volumetric concentration:')   
	sheet1.write_number(5,4, round(_M_vol_graph[-1], 2))
	sheet1.write(5,5, f'[m3/d]')  
	sheet1.write(6,0, f'Daily average Methane production:')  
	sheet1.write_number(6,4, round(M_production_average, 2))
	sheet1.write(6,5, f'[m3/d]')        
	sheet1.write(7,0, f'Radius particle growth:') 
	sheet1.write_number(7,4, round(radius_growth, 2))  
	sheet1.write(7,5, f'%')               
	sheet1.write(8,0, f'Final particle radius:')                   
	sheet1.write_number(8,4, round(_Rnew[-1], 2))
	sheet1.write(8,5, f'[mm]')
	end = time.time()
	sheet1.write(10, 0, f'Time elapsed for the simulation:')
	sheet1.write_number(10, 4, round(end-start, 2))
	sheet1.write(10, 5, f'seconds')
	sheet2 = wb.add_worksheet('Concentration profiles')
	sheet2.insert_image(1,1, '../output/Figures/Concentration profiles_'+str(date.tm_mday)+'_'+str(date.tm_mon)+'_'+str(date.tm_year)+'.png')

	sheet3 = wb.add_worksheet('Reactor model behaviour')
	sheet3.insert_image(1,1, '../output/Figures/Reactor model behaviour_'+str(date.tm_mday)+'_'+str(date.tm_mon)+'_'+str(date.tm_year)+'.png')
	wb.close()

	print('\n=======================================')
	print('Simulation finished --> Go to output folder')
	print('=======================================\n')
