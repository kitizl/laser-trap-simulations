#include <stdio.h>
#include <iostream>

#include <numeric> // to get std::iota from
#include <vector>
#include <cmath>
#include <random>

#include <string> // for string processing

#include <fstream>

const double PI = 3.14159265359;

float pToGamma(float pressure, float kBT, float a, float rho, float m){
    float p = pressure * 100;
	float pi = PI;
	float eta = 18.27e-6; // viscosity of air [Pa*s] at ambient temperature
	float m_a = 28.97e-3; // molecular mass of air [kg/mol]
	double N_A = 6.02214086e23;  // Avogadro constant [mol-1]
	float l = (eta*(std::sqrt((kBT* N_A*pi) /(2*m_a))))/p; // mean free path of air molecules
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

int main(int argc, char *argv[]){
    // A C++ Program specifically dedicated to simulating Linear Feedback Cooling
    std::cout << "Hello!" <<std::endl;
    if (argc < 3) {
        std::cout << "Not enough parameters specified, halting." << std::endl;
        return -1;
    }
    /**
        Parameters:
        1. Filename
        2. Charge in `e`
        3. Electric field in V/m
        4. Radius of particle in meters
        5. Pressure in mbar
        6. Frequency (f)
        7. Number of steps before saving (stick to 1.)
    */
    for (int i = 0 ; i < argc ; i++) {
        std::cout << i << "\t" << argv[i] << std::endl;
    }
    
    std::string filename = argv[1]; // extracting the filename
    
    float q = std::stof(argv[2])*1.6e-19; // extracting charge
    
    float E = std::stof(argv[3]); // extracting electric field
    // getting experimental parameters
    float density = 1860; // kg/m^3
    float kBT = 1.38e-23 * (273+22); // joules
    float radius = std::stof(argv[4]); // radius of particle
    std::cout << "Safe till now." << std::endl;
    float pressure = std::stof(argv[5]); // pressure
    
    float mass = (4/3)*PI*radius*radius*radius*density; // mass of the particle calculated from these params
    
    float gamma = pToGamma(pressure, kBT, radius, density, mass); // extracting gamma value through our fit.
    
    float omega = 2*PI*std::stof(argv[6]); // getting natural frequency of particle
    
    int steps_per_save = std::stoi(argv[7]); // number of steps per save

    std::ofstream myfile;
	myfile.open(filename); // create a file here
	myfile << "t,x,v" << std::endl;
    
    float max_time = 1.0; // 1 second
    double dt = 1e-6; // 1 microsecond sampling time
    long Nsteps = 1e6; // one million samples

    // initializing the position and velocity vectors
    float x = 100e-9;
    float v = 0;
    float a = 0;
    double t = 0;

    std::default_random_engine generator;
    double mean_force = 0.0; // average force has to be 0
    double std_force = std::sqrt(2*kBT*gamma/(dt*mass)); // standard deviation is propto the size of the autocorrelation function
  	std::normal_distribution<double> R(mean_force,std_force);

    float electric_force = q*E/mass; // initializing the force
    float noise_force = 0.0;

    for(int i = 1; i < Nsteps-1; i++){
        noise_force = R(generator);
        electric_force = q*E*v/mass;

        a = -omega*omega*x - gamma*v + noise_force + electric_force;
        v = v + 0.5*dt*a;
        x = x + dt*v;
        a = -omega*omega*x - gamma*v + noise_force + electric_force;
        v = v + 0.5*dt*a;

        if (std::abs(x) > 1.0){
            std::cout << "Particle lost" << std::endl;
            break;
        }
        if ( i % steps_per_save == 0 ) {
            //std::cout<< i << "," << t << "," << x << "," << v << std::endl;
            myfile << t << "," << x << "," << v << std::endl;
        }

        t = t+dt;
    }
    myfile.close();
    return 0;
}