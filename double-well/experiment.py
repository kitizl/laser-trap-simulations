#! python3

import numpy as np
import matplotlib.pyplot as plt
import time


def trapSolver(modAmplitude=0.0,f=0.0):
	# takes in some data and returns (t,x)
	T = 10000

	# a function that returns potentials
	def U(x,grad=False):
		a = 1
		b = 16
		if grad:
			return -a*x + b*x**3
		else:
			return -(a/2)*x**2 + (b/4)*x**4

	# measuring runtime
	start = time.time()

	numOfSteps = int(1e3)*T
	dt = T/numOfSteps

	x = np.zeros(numOfSteps)
	v = np.zeros(numOfSteps)

	# running the Euler steps + random noise + modulation
	for i in range(numOfSteps-1):
	    x[i+1] = x[i] + dt*(v[i])
	    v[i+1] = v[i] + dt*(-U(x[i],True) - 1*v[i] + 10*np.random.normal(0,1) + modAmplitude*np.cos(2*np.pi*f*dt*i)) 

	end = time.time()

	ts = np.linspace(0,T,numOfSteps)

	plt.plot(ts,x)
	plt.xlabel("Time")
	plt.ylabel("Position")

	print(f"Time taken to complete computation : {end-start} seconds.")

	return ts, x

# a series of "experiments" with different modulation amplitudes
t, unmodulated_x = trapSolver()
_, mod_x_1 = trapSolver(modAmplitude=0.01, f = 1/1000)
_, mod_x_2 = trapSolver(modAmplitude=0.01, f = 1/100)
_, mod_x_3 = trapSolver(modAmplitude=0.01, f = 1/10)
_, mod_x_4 = trapSolver(modAmplitude=0.01, f = 1)
_, mod_x_5 = trapSolver(modAmplitude=0.01, f = 10)



# saving the data

print("All position data collected.")

for i,pos in enumerate([unmodulated_x, mod_x_1, mod_x_2, mod_x_3, mod_x_4, mod_x_5]):
	np.save(f"exp2/positions-{i}.npy", pos)
np.save("exp2/timeseries.npy",t)

# delete these arrays later if the memory doesn't exist?

print("Position data saved.")
