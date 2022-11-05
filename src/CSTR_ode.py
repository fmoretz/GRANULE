import numpy as np

def fCSTR(y, t, Q, V, Sin, Xin, Ein, R, Y, Kd, eta, w):
	S, X, E, M = y
	omega = Q * w 

	dy = np.empty(4)
	dy[0] = Q/V * (Sin - S) - R*eta
	dy[1] = omega/V*(Xin - X) + R*eta*Y - Kd*X
	dy[2] = omega/V*(Ein - E) + Kd*X
	dy[3] = -Q/V*M + R*eta

	if S < 0:
		S = 0


	return dy
