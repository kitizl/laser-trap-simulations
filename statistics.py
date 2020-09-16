#!python3


import time
import numpy as np

def fourier_spectrum(X, sample_freq=1e-3):

	# takes in the signal and spits out the power spectrum
	# in arbitrary units

	ps = np.abs(np.fft.fft(X))**2
	freqs = np.fft.fftfreq(X.size, sample_freq)
	idx = np.argsort(freqs)

	return (freqs[idx], ps[idx])

def signaltonoise(fs, ps, f_d, dB=True):
	
	# returns the signal to noise ratio in decibels

	idx = np.where(fs==f_d)
	signal = ps[idx]

	# removing the edge cases because I do not 
	# want to deal with them
	if idx==0:
		return 0
	elif idx==len(fs):
		return 0

	noise = (1/2)*(ps[idx-1] + ps[idx+1])

	return 10*np.log10(signal/noise)

# writing a function to return residence time distribution
def residence_time(x,flags=1.0):
	# measure the time between it going between +flag and -flag
	# and add it to an array
	# return the array and/or produce a historgram so it also
	# gives you the distribution
	steps = len(x)
	crossing_indices = np.array([(0,"*")])

	tol = 1e-3

	dt = 1e-3 # took this value from the simulator

	for i in range(steps):
		if x[i] == 0:
			continue
		if (np.abs(x[i]-flags) < tol) and (x[i]-x[i-1] > 0):
			# this is when I'm breaching the positive crossing
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
	# this thing is O(n)

	crossing_indices = crossing_indices[1:]
	rtseries = np.array([dt*int(i) for i,_ in crossing_indices])

	# setting the first crossing to t=0

	rtseries = rtseries-rtseries[0]

	# finding the distribution by shifting the boys

	rtdistribution = rtseries[1:] - rtseries[:-1]

	return (crossing_indices, rtseries, rtdistribution)



def juliaStats():
	# performs statistics on data simulated in Julia
	# can be changed depending on the kind of outputs the simulator provides
	# works with .jld files
	import h5py
	
	for i in range(1,9):
		# importing the file handle
		start = time.time()

		print(f"Round {i} of 9")

		f = h5py.File(f"data-mod-{i}.jld","r")
		x = f["x"]

		print("\tData loaded.")

		# finding the power spectrum
		fs, ps = fourier_spectrum(x)
		print("\tPower spectrum found")
		# saving the power spectrum
		np.save(f"freq-{i}.npy",fs)
		np.save(f"power-spec-{i}.npy",ps)
		print("\tFiles saved.")
		end = time.time()

		print(f"Run {i} of 9. Time taken : {end-start} seconds.")

def pyStats():
	# performs statistics on data simulated in python
	# primarily works with .npy files
