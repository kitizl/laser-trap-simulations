#!python3
import numpy as np
import matplotlib.pyplot as plt
import time

# here we are going to simulate a single well oscillator

def potential(x, params, grad=False):
	"""
	A function that returns the value of the potential for
	a given position x, modulated by a list of parameters params.
	"""
	k = params # for now the only parameter is stiffness

	if grad: # if the gradient is required
		return k*x
	else: # if we need the potential directly
		return (k*x**2)/2

def var_stiffness(t,t0=2):
	"""
	A function where we can vary the stiffness of the trap
	as a function of time some given function
	"""
	if t<t0: # after some specific t0, stiffness shoots up
		return 0.5
	else:
		return 1.0


def trapSolver(params,save_frequency=2):
	"""
	Returns the times, positions and velocities of
	the nanosphere in the trap for a given set of 
	parameters (damping and so on) and a rate of saving data

	params : stiffness, max_time, damping constant gamma, k_B T
	save_frequency : number of steps between saving data
	"""

	def update_x(x,v,dt):
		# 'macro' to update position
		return x + v*dt/2
	def update_v(v,F,dt):
		# 'macro' to update velocity
		return v+F*dt/2

	def random_update(v,gamma,kBT, dt):
		# 'macro' to update velocity with the random noise
		R = np.random.normal()
		damping = np.exp(-gamma*dt) # this is following the BAOAB scheme
		random_kick = np.sqrt(1-damping*damping)*np.sqrt(kBT)
		return damping*v + R*random_kick

	stiffness, max_time, gamma, kBT = params # setting parameters
	# all initial conditions set to 0
	x = 0
	v = 0
	t = 0
	dt = 1e-3
	step_number = 0
	positions = []
	velocities= []
	total_energies = []
	save_times = []

	while(t<max_time):
		# B

		v = update_v(v,-potential(x,stiffness(t),True),dt)

		# A

		x = update_x(x,v,dt)

		# O

		v = random_update(v,gamma,kBT,dt)

		# A

		x = update_x(x,v,dt)

		# B

		v = update_v(v,-potential(x,stiffness(t),True),dt)

		if abs(step_number%save_frequency) <1e-6 and step_number>0:
			e_total = 0.5*v*v + potential(x,stiffness(0))

			positions.append(x)
			velocities.append(v)
			total_energies.append(e_total)
			save_times.append(t)
		t = t+dt
		step_number += 1
	return save_times, positions, velocities, total_energies

def save_data(data_series, file_index=0):
	"""
	Saves data in a .npy format for quick saving and loading
	"""
	np.save(f"./data/experiment-{file_index}.npy",np.vstack(data_series))

# stiffness, max time, gamma, kBT

