import numpy as np

def fCSTR(y, t, Q, V, Sin, Xin, Ein, R, Y, Kd, eta, U, A, rb, cp, m, Tin, Tw, w, Ro, D, rho_biomass):
	S, X, E, T, M = y
	omega = Q * w 

	dy = np.empty(5)
	dy[0] = (Q-omega)/V * (Sin - S) - R*eta
	dy[1] = omega/V*(Xin - X) + R*eta*Y - Kd*X
	dy[2] = omega/V*(Ein - E) + Kd*X
	dy[3] = ( U*A*(Tw - T) + Q*rb*cp*(Tin - T) ) / (m * cp)
	dy[4] = -Q/V*M + R*eta
	#dy[5] = D*Y/rho_biomass*(Sin - S)/(Ro + r) - r*Kd 
 
	if abs(Tw - T) < 1e-4:
		dy[3] = 0
	
	if S < 0:
		S = 0


	return dy
