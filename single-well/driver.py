#!python3

"""
This is the main driver code. Run this to get data.
"""
import matplotlib.pyplot as plt
import os
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
	# if we ever come to passing arguments to the driver
	# this is where you should put them
	return [10,0.1,0.1,2]

if __name__ == "__main__":
	print("Running Simulations!")
	# first, running simulations
	import simulation
	params = parse_arguments() # importing experimental parameters
	max_time, gamma, kBT, saving_freq = params
	DATA_DIR = "data_n10000_short_time"
	PLOT_DIR = "plot_n10000_short_time"
	
	make_dir(DATA_DIR) # making the directory where the data will end up

	NUM_TRIALS = 10000 # number of trials
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
	statistics.signal_ensemble(DATA_DIR,100,PLOT_DIR)
	# Providing directory where the data is, resolution (ie number of bins), directory where the plots end up

	# once this is done, you should have a folder of plot data and a movie

