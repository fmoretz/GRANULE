import numpy as np

def fTemperature(y, t, Q, U, A, rb, cp, m, Tin, Tw):
	T = y

	dy = ( U*A*(Tw - T) + Q*rb*cp*(Tin - T) ) / (m * cp)
 
	if abs(Tw - T) < 1e-4:
		dy = 0

	return dy
