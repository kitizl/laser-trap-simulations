#!python3

"""
This is script is to run a set of tests :
	1. Brownian Motion 
	2. Constant Stiffness
	3. Variable Stiffness
It is identical to driver.py in every other way.
"""
import os
import argparse
import sys
import datetime
import matplotlib.pyplot as plt
from math import pi, sqrt

def make_dir(dir_name):
	"""
	A function that creates a directory dir_name, if it doesn't exist.
	"""
	path = dir_name
	try:
		os.makedirs(path)
	except FileExistsError:
		print(f"Directory {path} already exists")
	except OSError:
		print (f"Creation of the directory %s failed" % path)
	else:
		print ("Successfully created the directory %s " % path)

def parse_arguments():
	# a helper function that parses the arguments passed
	# in the terminal

	# note : usage and help is included in the argparse package

	# constructing arugment parser

	all_args = argparse.ArgumentParser()

	# Add arguments to the parser
	all_args.add_argument("-t","--max_time",required=True,help="maximum time for simulation (in s)")
	all_args.add_argument("-p","--pressure",required=True,help="pressure (in mbar)")
	all_args.add_argument("-T","--temperature",required=True,help="temperature (in K)")
	all_args.add_argument("-s","--saving_freq",required=True,help="saving frequency (number of steps per save)")
	all_args.add_argument("-N","--numTrials", required=True,help="number of trials that need to be run")
	all_args.add_argument("-r","--resolution",required=True,help="resolution for the histogram")
	all_args.add_argument("-l","--label",required=False, help="label for the experiment")
	all_args.add_argument("-d", "--timestep",required=False, help="time step in simulation")

	args = vars(all_args.parse_args())

	# checking if the values are registered as intended
	print(f"Just to confirm, your values are : \n\
	Maximum Runtime : {args['max_time']} s\n\
	Pressure : {args['pressure']} mbar \n\
	Temperature : {args['temperature']} K \n\
	Saving Frequency : {args['saving_freq']} steps per save\n\
	Number of trials : {args['numTrials']} trials\n\
	Resolution : {args['resolution']}\n\
	Timestep: {args['timestep']} s\n")
	choice = input("Continue or abort? [y/n]")
	if choice.upper() == 'Y':
		# order is pressure, temperature, saving_freq, NUM_TRIALS, resolution
		return args
	else:
		print("Please run the script again.")
		sys.exit(0)

def physicalize(args):
	"""
		A function that takes in real world parameters such as 
		pressure, temperature, etc. and turns them into the
		numerical values required to solve the Langevin equation
	"""
	p = float(args['pressure']) * 100 # converting from mbar to Pa
	eta = 18.27e-6 # viscosity of air [Pa*s] at ambient temperature
	a = 71.5e-9  # particle radius [m]
	m_a = 28.97e-3 # molecular mass of air [kg/mol]
	N_A = 6.02214086e23  # Avogadro constant [mol-1]
	rho = 1850 # particle density kg*m^-3
	T = float(args['temperature']) # temperature [K]
	kB = 1.38064852e-23 # Boltzmann constant [m^2 kg s^-2 K^ -1]
	m = rho*(4/3)*pi*a**3 # mass of particle [kg]
	l = (eta*(sqrt((kB*T* N_A*pi) /(2*m_a))))/p # mean free path of air molecules
	Kn = l/a #knudsen number

	# Define the following constants for simplicity:
	A = 6*pi*eta*a/(m)
	ck = 0.31*Kn/(0.785 + 1.152*Kn + Kn**2)
	B = (1+ck)

	# damping rate
	Gamma = A*(0.619/(0.619+Kn))*B
	Gamma2 = Gamma / (2*pi)

	return [float(args['max_time']), m,Gamma2, kB*T, int(args['saving_freq']), int(args['numTrials']), int(args['resolution']),args['label'],float(args['timestep'])]

if __name__ == "__main__":

	import simulation
	import statistics
	params = parse_arguments() # importing experimental parameters
	max_time, mass, gamma, kBT, saving_freq, NUM_TRIALS, resolution, label, timestep = physicalize(params)


	# directory format is YYYYMMDD_HHMM

	# the datetime object is formmated to provide the the string
	# we desire for our directories
	# crucial assumption here being two RUNS will not complete
	# within the minute
	timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M")
	DATA_DIR = "data"
	PLOT_DIR = "plot"

	if label is None:
		# if a label is not provided, then use default naming conventions
		DATA_DIR += f"_n{NUM_TRIALS}_{timestamp}"
		PLOT_DIR += f"_n{NUM_TRIALS}_{timestamp}"
	else:
		# if a label is provided, then use label for directory naming
		DATA_DIR += f"_{label}"
		PLOT_DIR += f"_{label}"
	
	test_labels= ["brownian","const-stiffness","variable-stiffness"]
	for label in test_labels:
		make_dir(DATA_DIR+f"/{label}")
		make_dir(PLOT_DIR+f"/{label}")

	print("Running Simulations!")
	# first, running simulations
	print(f"Maximum simulation time\t:{max_time}\nDamping (gamma)\t:{gamma}\nTemperature (kB T)\t:{kBT}\nSaving frequency\t:{saving_freq} steps per save")
	# doing this multiple times so as to generate an average

	print("====== BROWNIAN TEST ======")
	# BROWNIAN TEST
	for trial_num in range(NUM_TRIALS):
		print(f"\r{trial_num}/{NUM_TRIALS}",end="")
		ts,xs,vs,ks,ps = simulation.trapSolver([lambda t : 0,mass,max_time, gamma, kBT], timestep, saving_freq)
		# lambda function to return 0 for all values of t
		# save this data
		simulation.save_data([ts,xs,vs,ks,ps],DIR_NAME=DATA_DIR+"/brownian",file_index=trial_num)
	print("\n")
	print("Producing signal distribution")
	statistics.signal_ensemble(DATA_DIR+"/brownian",resolution,PLOT_DIR+"/brownian")
	# Providing directory where the data is, resolution (ie number of bins), directory where the plots end up
	print("Producing signal distribution")
	statistics.signal_variance(DATA_DIR+"/brownian",PLOT_DIR+"/brownian")

	print("Producing velocity distribution")
	statistics.velocity_distribution(mass, kBT, DATA_DIR+"/brownian", resolution, PLOT_DIR+"/brownian")

	# Providing directory where the data is, resolution (ie number of bins), directory where the plots end up
	print("Producing energy evolution plots")
	statistics.energy_evolution(DATA_DIR+"/brownian",PLOT_DIR+"/brownian")
	# once this is done, you should have a folder of plot data and a movie
	print("Analysis complete!")
	
	print("===== CONSTANT STIFFNESS TEST =====")
	# CONSTANT STIFFNESS
	for trial_num in range(NUM_TRIALS):
		print(f"\r{trial_num}/{NUM_TRIALS}",end="")
		ts,xs,vs,ks,ps = simulation.trapSolver([lambda t : simulation.gieseler_stiffness(),mass,max_time, gamma, kBT], timestep, saving_freq)
		# lambda function to return a constant k for all values of t
		# save this data
		simulation.save_data([ts,xs,vs,ks,ps],DIR_NAME=DATA_DIR+"/const-stiffness",file_index=trial_num)
	print("\n")
	print("Producing signal distribution")
	statistics.signal_ensemble(DATA_DIR+"/const-stiffness",resolution,PLOT_DIR+"/const-stiffness")
	# Providing directory where the data is, resolution (ie number of bins), directory where the plots end up
	print("Producing signal distribution")
	statistics.signal_variance(DATA_DIR+"/const-stiffness",PLOT_DIR+"/const-stiffness")

	print("Producing velocity distribution")
	statistics.velocity_distribution(mass, kBT, DATA_DIR+"/const-stiffness", resolution, PLOT_DIR+"/const-stiffness")


	# Providing directory where the data is, resolution (ie number of bins), directory where the plots end up
	print("Producing energy evolution plots")
	statistics.energy_evolution(DATA_DIR+"/const-stiffness",PLOT_DIR+"/const-stiffness")
	# once this is done, you should have a folder of plot data and a movie
	print("Analysis complete!")


	print("===== VARIABLE STIFFNESS TEST =====")
	# VARIABLE STIFFNESS
	for trial_num in range(NUM_TRIALS):
		print(f"\r{trial_num}/{NUM_TRIALS}",end="")
		ts,xs,vs,ks,ps = simulation.trapSolver([simulation.var_stiffness,mass,max_time, gamma, kBT], timestep, saving_freq)
		# save this data
		simulation.save_data([ts,xs,vs,ks,ps],DIR_NAME=DATA_DIR+"/variable-stiffness",file_index=trial_num)
	print("\n")
	print("Producing signal distribution")
	statistics.signal_ensemble(DATA_DIR+"/variable-stiffness",resolution,PLOT_DIR+"/variable-stiffness")
	# Providing directory where the data is, resolution (ie number of bins), directory where the plots end up
	print("Producing signal distribution")
	statistics.signal_variance(DATA_DIR+"/variable-stiffness",PLOT_DIR+"/variable-stiffness")
	# Providing directory where the data is, resolution (ie number of bins), directory where the plots end up

	print("Producing velocity distribution")
	statistics.velocity_distribution(mass, kBT, DATA_DIR+"/variable-stiffness", resolution, PLOT_DIR+"/variable-stiffness")

	print("Producing energy evolution plots")
	statistics.energy_evolution(DATA_DIR+"/variable-stiffness",PLOT_DIR+"/variable-stiffness")
	# once this is done, you should have a folder of plot data and a movie
	print("Analysis complete!")
	
	print("All tests are complete...")

	# once all of this is complete, we then add a METADATA file
	# that contains the parameters that went into producing these
	# numerical results
	print("Saving ")
	with open(f"{DATA_DIR}/METADATA-experiment-{timestamp}.txt",'w') as fhand:
		fhand.write(f"# Metadata\n\
## Experimental Parameters\n\
	Maximum runtime : {params['max_time']}\n\
	Pressure : {params['pressure']} mbar \n\
	Temperature : {params['temperature']} K \n\
	Saving Frequency : {params['saving_freq']} steps per save\n\
	Number of trials : {params['numTrials']} trials\n\
	Resolution : {params['resolution']}\n\
	---\n\
## Directories\n\
	Data stored in : {DATA_DIR}\n\
	Plots stored in : {PLOT_DIR}\n\
	---\n\
## Additional Comments\n\
	<add comments if you want>\n")
		fhand.close()

