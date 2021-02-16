#!python3

import numpy as np


"""
	This code implements the algorithm as described in the McNamara and Wiesenfeld paper : https://doi.org/10.1103/PhysRevA.39.4854
"""

def randnum(mode='uni'):
	# Returns a random number of a specific distribution

	# Parameters :
	#	mode (str):The distribution of choice (currently only uniform and normal)

	# Returns:
	#	random_number from 0 to 1

	if mode=='uni':
		return np.random.uniform()
	elif mode=='norm':
		return np.random.normal()
	else:
		return randnum()

def transition_rate(alpha_0, alpha_1, eta_0, omega_s, t, state=+1):
	# Returns the transition rate for a given system

	# Parameters:
	#	alpha_0 : Ratio of potential barrier to noise strength
	# 	alpha_1 and eta_0 : Parameters that represent modulation amplitude
	# 	omega_s : Modulation frequency
	#	state : Specifies which state the system is in
	return (1/2)*(alpha_0 + state*(alpha_1*eta_0*np.cos(omega_s*t)))

# time interval
dt = 0.005
# modulating frequency
f_s = 0.001

# specifying the parameters
alpha_0 = 10
alpha_1 = 1
eta_0 = 5

x = np.zeros(np.arange(0,100,dt).shape)


# starts at the x=-c state
curr_state = -1	

for i,t in enumerate(np.arange(0,100,dt)):
	xi = randnum() # produces random number

	p = dt*transition_rate(alpha_0,alpha_1,eta_0,2*np.pi*f_s,t,state=curr_state) # calculates transition probability
	
	x[i] = curr_state # records current state

	if xi < p:
		# checks for transition
		curr_state = curr_state*-1 # switching states

# saving the file
np.save("twostate.npy",x)
