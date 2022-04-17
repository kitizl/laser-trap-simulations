#include <iostream>

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
        2. Radius of particle in meters
        3. Pressure in mbar
        4. Frequency (f)
        5. Number of steps before saving (stick to 1.)
		6. Modulation Depth (Epsilon)
		7. Modulation Phase Shift (Phi)
    */
    for (int i = 0 ; i < argc ; i++) {
        std::cout << i << "\t" << argv[i] << std::endl;
    }
    
    std::string filename = argv[1]; // extracting the filename

    // getting experimental parameters
    float density = 1860; // kg/m^3
    float kBT = 1.38e-23 * (273+22); // joules
    float radius = std::stof(argv[2]); // radius of particle
    std::cout << "Safe till now." << std::endl;
    float pressure = std::stof(argv[3]); // pressure
    
    float mass = (4/3)*PI*radius*radius*radius*density; // mass of the particle calculated from these params
    
    float gamma = pToGamma(pressure, kBT, radius, density, mass); // extracting gamma value through our fit.
    float omega = 2*PI*std::stof(argv[4]); // getting natural frequency of particle
    
    int steps_per_save = std::stoi(argv[5]); // number of steps per save

    std::ofstream myfile;
	myfile.open(filename); // create a file here
	myfile << "t,x,v" << std::endl;
    
    float max_time = 1.0; // 1 second
    double dt = 1e-6; // 1 nanosecond sampling time
    long Nsteps = 1e6; // one billion samples

    // initializing the position and velocity vectors
    float x = 100e-9;
    float v = 0;
    float a = 0;
    double t = 0;

    std::default_random_engine generator;
  	std::normal_distribution<double> R(0.0,1.0);

    float noise_force = 0.0;

	// getting modulation depth

	float mod_depth = std::stof(argv[6])
	float param_phase_shift = std::stof(argv[7])

    myfile << t << "," << x << "," << v << std::endl;
    for(int i = 1; i < Nsteps-1; i++){
        noise_force = R(generator)*std::sqrt(2*kBT*gamma/(mass));

        a = -omega*omega*(1 + mod_depth*std::cos(2*omega*t + param_phase_shift))*x - gamma*v + noise_force;
    
        v = v + 0.5*dt*a;
        x = x + dt*v;
        a = -omega*omega*(1 + mod_depth*std::cos(2*omega*t + param_phase_shift))*x - gamma*v + noise_force;
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