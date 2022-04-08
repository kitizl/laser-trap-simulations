#include <stdio.h>
#include <iostream>

#include <numeric> // to get std::iota from
#include <vector>
#include <cmath>
#include <random>

#include <fstream>
// Implementing the BAOAB algorithm for a stochastic modulated Duffing oscillator

const double PI = 3.14159265358979323;


// stolen trying to implement arange
template<typename T>
std::vector<T> arange(T start, T stop, T step = 1) {
    std::vector<T> values;
    for (T value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}

float pressure_to_gamma(float pressure, float kBT){
	// input in mbar
	float p = pressure * 100;
	float pi = PI;
	float eta = 18.27e-6; // viscosity of air [Pa*s] at ambient temperature
	float a = 71.5e-9;  // particle radius [m]
	float m_a = 28.97e-3; // molecular mass of air [kg/mol]
	float N_A = 6.02214086e23;  // Avogadro constant [mol-1]
	float rho = 1850; // particle density kg*m^-3
	float m = rho*(4/3)*pi*a*a*a; // mass of particle [kg]
	float l = (eta*(sqrt((kBT* N_A*pi) /(2*m_a))))/p; // mean free path of air molecules
	float Kn = l/a; //knudsen number

	// Define the following constants for simplicity:
	float A = 6*pi*eta*a/(m);
	float ck = 0.31*Kn/(0.785 + 1.152*Kn + Kn*Kn);
	float B = (1+ck);

	// damping rate
	float Gamma = A*(0.619/(0.619+Kn))*B;
	float Gamma2 = Gamma / (2*pi);
	return Gamma2;
}

int main(){
	std::ofstream myfile;
	myfile.open("electric-driving-simulations/qE_sine.csv");
	myfile << "t,x,v\n";

	std::cout << "Setting up simulations..." << std::endl;

	double x = 0.0;
	double v = 0.0;
	double t = 0.0;
	long long i = 0;
	double max_time = 0.5; // seconds
	double dt = 1e-6; // 1 microsecond

	// experimental parameters
	float density = 1850;
	float radius = (143/2)*1e-9;
	float mass = (4*PI/3)*radius*radius*radius*density;
	float omega0 = 2*PI*153.646e3; // setting nat. frequency to 100 kHz
	float kBT = 1.38e-23 * 300; // 300 K temperature
	float gamma = 15602.0;// gamma from the calibrator
	float q = -7*1.6e-19;
	float E = +70000; // V/m

	std::default_random_engine generator;
  	std::normal_distribution<double> R(0.0,1.0);

  	std::cout << "Starting simulations" <<std::endl;
	while (t < max_time) {

		// B
		v = v - (dt/2)*(omega0*omega0*x + q*E*std::sin((6.28*200e3)*t));
		// A
		x = x + (dt/2)*v;
		// O
		
		double damping = std::exp(-gamma*dt);// through solving stochastic eq.
		double random_kick = std::sqrt(1-damping*damping)*std::sqrt(kBT/mass);
		v = damping*v + R(generator)*random_kick; // main O-step

		// A
		x = x + (dt/2)*v;
		// B
		v = v - (dt/2)*(omega0*omega0*x + q*E*std::sin((6.28*200e3)*t));

		// save t,x,v
		if (x > 1e100) {
			std::cout << "Particle lost" <<std::endl;
		}
		std::cout<< t << "," << x << "," << v << std::endl;
		myfile <<  t << "," << x << "," << v << std::endl;
		t = t + dt; // updating time
		i++; // incrementing counter

	}
	myfile.close();
	return 0;
}
