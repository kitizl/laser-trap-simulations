#include <stdio.h>
#include <iostream>

int main(){
	// xdot = v
	// vdot = -gamma_0 * v - omega_0^2 * x + (1/m) (F_th + F_fb)
	// each variable is progressed one half-timestep out of sync
	// x_{n+1} = x_n + v_{n + 0.5} dt
	// v_{n + 0.5} = v_{n-0.5} + dot v_N dt
	// F_th is a Gaussian random number with 0 mean and variance of (2mkT gamma/dt)
	// Measurement noise is simulated by S_nn/dt which is uncorrelated white noise
	// gamma_0 = (1 + pi/8) * 4pi/3 * MNR^2 v_T/m
	// vT = sqrt(2kT/pi M)
	return 1;
}