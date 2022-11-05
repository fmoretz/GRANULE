import numpy as np

def fCSTR(y, t, Q, V, Sin, Xin, Ein, R, Y, Kd, eta, w):
	S, X, E, M = y
	omega = Q * w 

	dy = np.empty(4)
	dy[0] = - R*eta
	dy[1] = + R*eta*Y - Kd*X
	dy[2] =  + Kd*X
	dy[3] = -Q/V*M + R*eta

	if S < 0:
		S = 0


	return dy
