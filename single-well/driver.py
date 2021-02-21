#!python3

"""
This is the main driver code. Run this to get data.
"""
import matplotlib.pyplot as plt

def parse_arguments():
	# if we ever come to passing arguments to the driver
	# this is where you should put them
	return [10,0.1,0.1,10]

if __name__ == "__main__":
	print("Running Simulations!")
	# first, running simulations
	import simulation
	params = parse_arguments() # importing experimental parameters
	max_time, gamma, kBT, saving_freq = params

	print(f"Maximum simulation time\t:{max_time}\nDamping (gamma)\t:{gamma}\nTemperature (kB T)\t:{kBT}\nSaving frequency\t:{saving_freq} steps per save")
	# doing this multiple times so as to generate an average
	for trial_num in range(100):
		print(f"Running {trial_num}/100...")
		ts,xs,vs,es = simulation.trapSolver([simulation.var_stiffness,max_time, gamma, kBT], saving_freq)
		# save this data
		simulation.save_data([ts,xs,vs,es],DIR_NAME="data",file_index=trial_num)
	
	print("Simulations are complete")

	# running statistics, particularly the histogram
	import statistics
	statistics.signal_ensemble("data",20,"plots") # providing DIR for accessing data

	# once this is done, you should have a folder of plot data and a movie

