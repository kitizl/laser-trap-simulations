#!python3


import time
import numpy as np
import glob
import matplotlib.pyplot as plt

def fourier_spectrum(X, sample_freq=1e-3):
	"""
	Returns the power spectrum (in arbitary units) for a given
	signal X and a sampling frequency sample_freq
	"""
	# find FFT of signal and its amplitude
	ps = np.abs(np.fft.fft(X))**2 
	# store frequencies and indices corresponding to the above FFT
	freqs = np.fft.fftfreq(X.size, sample_freq)
	idx = np.argsort(freqs)

	return (freqs[idx], ps[idx])

def signaltonoise(fs, ps, f_d):
	"""
	Returns the signal to noise ratio in dB for a given range of frequencies
	fs for a powerspectrum ps, specifically for a given frequency
	"""

	# find the index of the specified frequency
	idx = np.where(fs==f_d) 
	# find the value of the spectrum at that frequency
	signal = ps[idx]

	# removing edge cases of when the indices are at the edge
	# of the spectrum
	if idx==0:
		return 0
	elif idx==len(fs):
		return 0

	# calculating noise by averaging/interpolating
	# the neighbouring spectral values from the 
	# specified frequency
	noise = (1/2)*(ps[idx-1] + ps[idx+1])

	return 10*np.log10(signal/noise)


def residence_time(x,flags=1.0):
	"""
	Returns the residence time distribution for a given signal x
	and for a given set of flags (where the timekeeping is triggered)
	"""

	# measure the time between it going between +flag and -flag
	# and add it to an array
	# return the array and/or produce a historgram so it also
	# gives you the distribution
	steps = len(x)

	# keeping track of the crossing indices, * being the null label
	crossing_indices = np.array([(0,"*")])

	tol = 1e-3 # tolerance for when signal crosses the flag

	dt = 1e-3 # took this value from the simulator

	for i in range(steps):
		if x[i] == 0:
			continue
		if (np.abs(x[i]-flags) < tol) and (x[i]-x[i-1] > 0):
			# this is when the signal is breaching the positive crossing
			if len(crossing_indices) > 1:
				# if this isn't the first time a crossing occured
				if crossing_indices[-1][-1] == "-":
					# and if the last crossing was at the other
					# level crossing
					crossing_indices = np.vstack((crossing_indices,[(i,"+")]))
					# append this to the new crossing
			else:
				crossing_indices = np.vstack((crossing_indices,[(i,"+")]))
			# otherwise, ignore completely
		elif (np.abs(x[i]+flags) < tol) and (x[i]-x[i-1] < 0):
			# this is the negative crossing
			if len(crossing_indices) > 1:
				# if this isn't the first time a crossing occured
				if crossing_indices[-1][-1] == "+":
					# and if the last crossing was at the other
					# level crossing
					crossing_indices = np.vstack((crossing_indices,[(i,"-")]))
					# append this to the new crossing
			else:
				crossing_indices = np.vstack((crossing_indices,[(i,"-")]))
			# otherwise, ignore completely


	crossing_indices = crossing_indices[1:] # ignoring the first crossing index since it's (0,*)
	rtseries = np.array([dt*int(i) for i,_ in crossing_indices])

	# setting the first crossing to t=0
	rtseries = rtseries-rtseries[0]

	# finding the distribution by subtracting from an offset timing series
	rtdistribution = rtseries[1:] - rtseries[:-1]

	# returns the indices when the flags were triggered, the timing series and the distribution overall
	return (crossing_indices, rtseries, rtdistribution)


def signal_ensemble(data_loc):
	"""
	A function that returns heatmap of the position distributions
	"""
	# making a list of all the files
	file_list = glob.glob(f"{data_loc}/experiment*")
	# importing all of the data from the experiments
	print("Importing data...")
	all_data = np.array([np.load(file) for file in file_list])
	# extracting time series (assumes common time scaling across exps)
	ts = all_data[0][0]
	# extracting all position data
	pos_data = np.array([all_data[i][1] for i in range(len(all_data))])

	print("Producing plots...")
	# making histogram plots
	for i in range(len(ts)):
		print(f"\r{i}/{len(ts)}",end="")
		plt.clf() # clear figure
		plt.xlim(-1,1) # setting common x axis
		# we are taking the histogram across experiments
		# for each timestep, hence the transposing
		plt.hist(pos_data.T[i]) # plotting histogram
		plt.title(f"Time : {ts[i]} units") # keeping track of time
		plt.savefig(f"./plots/step-{i:05n}.png")
	print("\nPlot production complete!")
