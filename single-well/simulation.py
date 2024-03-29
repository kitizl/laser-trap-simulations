#!python3
import numpy as np
import matplotlib.pyplot as plt
import time
from helper import make_dir
from sys import exit

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

def electric_field_amplitude(power, waist):
	"""
		Returns electric field amplitude given the power of the laser
		and the waist at the focus.

		Since electric field amplitude is not directly accessible t
	"""
	permittivity = 8.8541878128e-12 # SI units, from CODATA
	c = 299792458 # SI units, from CODATA
	return (4*power)/(permittivity*c*np.pi*waist**2)

def gieseler_stiffness():
	"""
		A physically realistic function that determines stiffness
		of the trap at time t (constant for now) using constants
		determined by Jan Gieseler in their PhD thesis.
	"""
	laser_power = 300e-3 # 300 mW is the power of the laser
	waist = 687e-9 # meters
	perm =  8.8541878128e-12 # farad/meter
	wave = 2*np.pi/(1064e-9) # wavenumber of the laser
	a = 71.5e-9  # particle radius [m]
	V = 4/3*np.pi*a**3 # volume of particle (assuming full spherical symmetry)
	# in air, dielectric constant of medium is basically the same as permittivity
	em = perm
	# but that of the particle is
	ep = 11.9*perm #:shrug:
	alpha = 3*V*perm*(ep-em)/(ep+2*em) # Clausius-Mossotti relation
	pol_eff = alpha/(1 + (((wave**3)/(6*np.pi*perm))**2) * (alpha**2)) # calculating (real) effective polarizability

	k = pol_eff*electric_field_amplitude(laser_power,waist)/waist**2 # calculating stiffness from polarizability, laser intensity and waist width
	return k

def square_stiffness(t, omega):
	from scipy.signal import square

	k_min = giessen_stiffness()
	k_max = 10*k_min
	return (1/2)*((k_max + k_min)*square(omega*t) + (k_max - k_min))

def var_stiffness(t,t0=0.5):
	"""
	A function where we can vary the stiffness of the trap
	as a function of time some given function
	"""
	# using a logistic (sigmoid) function
	# k_max is maximum value
	# steep is the steepness
	# t0 is when the stiffness is halfway to the max value
	k_min = gieseler_stiffness()
	k_max = 10*k_min
	steep = 10.0
	return k_min + (k_max-k_min)/(1+np.exp(-steep*(t-t0)))

def langevinEuler(params, save_frequency=2):
	"""
	DO NOT USE THIS FUNCTION : Use trapSolver() INSTEAD!

	Uses Euler-Maruyama scheme to integrate
	Returns the times, positions and velocities of
	the nanosphere in the trap for a given set of 
	parameters (damping and so on) and a rate of saving data

	params : stiffness, max_time, damping constant gamma, k_B T
	save_frequency : number of steps between saving data
	"""
	stiffness, mass, max_time, gamma, kBT = params # setting parameters
	# all initial conditions set to 0
	x = 1e-8 # initial position for time being
	v = 0
	t = 0
	dt = 1e-6
	step_number = 0
	positions = []
	velocities= []
	total_energies = []
	save_times = []

	while t < max_time:

		a = (-potential(x, stiffness(t), True))/mass # acceleration
		x = x + v*dt # updating position
		v = v + a*dt # updating velocity

		if abs(step_number%save_frequency) <1e-6 and step_number>0:
			e_total = 0.5*v*v + potential(x,stiffness(t))

			positions.append(x)
			velocities.append(v)
			total_energies.append(e_total)
			save_times.append(t)
		t = t+dt
		step_number += 1
	return save_times, positions, velocities, total_energies

def trapSolver(params,timestep,save_frequency=2):
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

	def random_update(v,gamma,mass,kBT, dt):
		# 'macro' to update velocity with the random noise
		R = np.random.normal()
		damping = np.exp(-gamma*dt) # through solving stochastic eq.
		random_kick = np.sqrt(1-damping*damping)*np.sqrt(kBT/mass)
		return damping*v + R*random_kick

	stiffness, mass, max_time, gamma, kBT = params # setting parameters
	# all initial conditions set to 0
	x = 0 # initial position for time being
	v = 0
	t = 0
	dt = timestep # setting timestep from user param
	step_number = 0
	positions = []
	velocities= []
	kinetic_energy = []
	potential_energy = []
	save_times = []

	while(t<max_time):
		# B

		v = update_v(v,-potential(x,stiffness(t),True)/mass,dt)

		# A

		x = update_x(x,v,dt)

		# O

		v = random_update(v,gamma,mass,kBT,dt)

		# A

		x = update_x(x,v,dt)

		# B

		v = update_v(v,-potential(x,stiffness(t),True)/mass,dt)

		v = update_v(v,0,dt)
		if abs(step_number%save_frequency) <1e-6 and step_number>0:
			kin = 0.5*mass*v*v # calculating kinetic energy
			pot = potential(x,stiffness(t)) # calculating potential energy

			positions.append(x)
			velocities.append(v)
			kinetic_energy.append(kin)
			potential_energy.append(pot)
			save_times.append(t)
		t = t+dt
		step_number += 1

		if np.abs(x) > 1e100:
			# if the "distance" of the particle is too large
			print(x)
			print("Nanosphere lost.")
			exit(-1) # exiting now.
	return save_times, positions, velocities, kinetic_energy, potential_energy

def save_data(data_series, DIR_NAME="data", file_index=0):
	"""
	Saves data in a .npy format for quick saving and loading
		data_series : the dataset that needs to be saved
		DIR_NAME : the directory where the data needs to be saved
		file_index : an index to label individual runs
	"""
	# saves each trial run in the specific directory
	np.save(f"{DIR_NAME}/experiment-{file_index}.npy",np.vstack(data_series))

def duffingSolver(params,timestep,save_frequency=2):
	"""
	Returns the times, positions and velocities of
	the nanosphere in the trap for a given set of 
	parameters (damping and so on) and a rate of saving data,
	modelling the system as a stochastic Duffing Oscillator

	params : stiffness, max_time, damping constant gamma, k_B T
	save_frequency : number of steps between saving data
	"""

	def update_x(x,v,dt):
		# 'macro' to update position
		return x + v*dt/2
	def update_v(v,F,dt):
		# 'macro' to update velocity
		return v+F*dt/2

	def random_update(v,gamma,mass,kBT, dt):
		# 'macro' to update velocity with the random noise
		R = np.random.normal()
		damping = np.exp(-gamma*dt) # through solving stochastic eq.
		random_kick = np.sqrt(1-damping*damping)*np.sqrt(kBT/mass)
		return damping*v + R*random_kick

	stiffness, mass, max_time, gamma, kBT = params # setting parameters
	# all initial conditions set to 0
	x = 0 # initial position for time being
	v = 0
	t = 0
	dt = timestep # setting timestep from user param
	step_number = 0
	positions = []
	velocities= []
	kinetic_energy = []
	potential_energy = []
	save_times = []

	while(t<max_time):
		# B

		v = update_v(v,-potential(x,stiffness(t),True)/mass,dt)

		# A

		x = update_x(x,v,dt)

		# O

		v = random_update(v,gamma,mass,kBT,dt)

		# A

		x = update_x(x,v,dt)

		# B

		v = update_v(v,-potential(x,stiffness(t),True)/mass,dt)

		v = update_v(v,0,dt)
		if abs(step_number%save_frequency) <1e-6 and step_number>0:
			kin = 0.5*mass*v*v # calculating kinetic energy
			pot = potential(x,stiffness(t)) # calculating potential energy

			positions.append(x)
			velocities.append(v)
			kinetic_energy.append(kin)
			potential_energy.append(pot)
			save_times.append(t)
		t = t+dt
		step_number += 1

		if np.abs(x) > 1e100:
			# if the "distance" of the particle is too large
			print(x)
			print("Nanosphere lost.")
			exit(-1) # exiting now.
	return save_times, positions, velocities, kinetic_energy, potential_energy
