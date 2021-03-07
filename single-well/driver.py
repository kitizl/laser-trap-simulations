#!python3

"""
This is the main driver code. Run this to get data.
"""
import matplotlib.pyplot as plt
import os
import argparse
import sys
import datetime

def make_dir(dir_name):
	"""
	A function that creates a directory dir_name, if it doesn't exist.
	"""
	path = dir_name
	try:
		os.mkdir(path)
	except FileExistsError:
		print(f"Director {path} already exists")
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

	# TODO : add optional argument to turn verbose mode off
	#		 which is currently the default
	# TODO: add optional argument for labelling the directories
	#		default being : `dir_n{NUM_TRIALS}_t{max_time}_datetime`

	args = vars(all_args.parse_args())

	# checking if the values are registered as intended
	print(f"Just to confirm, your values are : \n\
	Maximum Runtime : {args['max_time']} s\n\
	Pressure : {args['pressure']} mbar \n\
	Temperature : {args['temperature']} K \n\
	Saving Frequency : {args['saving_freq']} steps per save\n\
	Number of trials : {args['numTrials']} trials\n\
	Resolution : {args['resolution']}\n")
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
	# defining Boltzmann constant

	kB = 1.381e-23 # J/K ; from CODATA 2019
	return [float(args['max_time']), 10, kB*float(args['temperature']), int(args['saving_freq']), int(args['numTrials']), int(args['resolution'])]

if __name__ == "__main__":

	import simulation
	params = parse_arguments() # importing experimental parameters
	max_time, gamma, kBT, saving_freq, NUM_TRIALS, resolution = physicalize(params)

	# directory format is YYYYMMDD_HHMM

	# the datetime object is formmated to provide the the string
	# we desire for our directories
	# crucial assumption here being two RUNS will not complete
	# within the minute
	timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M")\

	DATA_DIR = f"data_n{NUM_TRIALS}_{timestamp}"
	PLOT_DIR = f"plot_n{NUM_TRIALS}_{timestamp}"

	make_dir(DATA_DIR) # making the directory where the data will end up

	print("Running Simulations!")
	# first, running simulations
	print(f"Maximum simulation time\t:{max_time}\nDamping (gamma)\t:{gamma}\nTemperature (kB T)\t:{kBT}\nSaving frequency\t:{saving_freq} steps per save")
	# doing this multiple times so as to generate an average
	for trial_num in range(NUM_TRIALS):
		print(f"\r{trial_num}/{NUM_TRIALS}",end="")
		ts,xs,vs,es = simulation.trapSolver([simulation.var_stiffness,max_time, gamma, kBT], saving_freq)
		# save this data
		simulation.save_data([ts,xs,vs,es],DIR_NAME=DATA_DIR,file_index=trial_num)
	
	print("\nSimulations are complete")

	# running statistics, particularly the histogram
	import statistics
	print("Analyzing simulations...")

	print("Producing signal distribution")
	statistics.signal_ensemble(DATA_DIR,resolution,PLOT_DIR)
	# Providing directory where the data is, resolution (ie number of bins), directory where the plots end up
	print("Producing energy evolution plots")
	statistics.energy_evolution(DATA_DIR,PLOT_DIR)
	# once this is done, you should have a folder of plot data and a movie

	# once all of this is complete, we then add a METADATA file
	# that contains the parameters that went into producing these
	# numerical results

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

