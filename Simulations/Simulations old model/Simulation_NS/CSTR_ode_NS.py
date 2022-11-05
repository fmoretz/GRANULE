import numpy as np

def fCSTR(t, y, Q, V, Sin, R, Y, Kd, eta, U, A, rb, cp, m, Tin, Tw, Na, Ks, _m, Xo, wi):

	S, X, E, T, M = y


	dy = np.empty(5)
	dy[0] = Q/V * (Sin - S) - (_m*X/Y)*(S/(Ks + S))*eta
	dy[1] = Q/V * (Xo - X) + (_m*S/(Ks + S)-Kd)*X
	dy[2] = Kd*X
	dy[3] = ( U*A*(Tw - T) + Q*rb*cp*(Tin - T) ) / (m * cp)
	dy[4] = -Q/V*M + R*eta
 
	if abs(Tw - T) < 1e-4:
		dy[3] = 0
	
	if S < 0:
		S = 0


	return dy
