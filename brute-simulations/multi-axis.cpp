#include <iostream>

#include <cmath>
#include <random>

#include <string> // for string processing

#include <fstream> // for file IO

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
        4. x-axis frequency (fx)
        5. y-axis frequency (fy)
        6. Number of steps before saving (stick to 1.)
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
    float omega_x = 2*PI*std::stof(argv[4]); // getting x-frequency of particle
    float omega_y = 2*PI*std::stof(argv[5]); // getting y-frequency of particle
    
    int steps_per_save = std::stoi(argv[6]); // number of steps per save

    std::ofstream myfile;
	myfile.open(filename); // create a file here
	myfile << "t,x,y,vx,vy" << std::endl;
    
    float max_time = 1.0; // 1 second
    double dt = 1e-6; // 1 microsecond sampling time
    long Nsteps = max_time/dt;

    // initializing the position and velocity vectors
    float x = 0.0;
    float y = 0.0;
    float vx = 0;
    float vy = 0;
    float ax = 0;
    float ay = 0;
    double t = 0;

    std::default_random_engine generator;
  	std::normal_distribution<double> R(0.0,1.0);

    float noise_force_x = 0.0;
    float noise_force_y = 0.0; // noise has to be uncorrelated between axes

    myfile << t << "," << x << "," << y << "," << vx << "," << vy << std::endl;
    // updating the first line first
    for(int i = 1; i < Nsteps-1; i++){
        noise_force_x = R(generator)*std::sqrt(2*kBT*gamma/(mass));
        noise_force_y = R(generator)*std::sqrt(2*kBT*gamma/(mass));
        // find acceleration first        
        ax = -omega_x*omega_x*x - gamma*vx + noise_force_x;
        ay = -omega_y*omega_y*y - gamma*vy + noise_force_y;
        // update velocity 
        vx = vx + 0.5*dt*ax;
        vy = vy + 0.5*dt*ay;

       // update position
        x = x + dt*vx;
        y = y + dt*vy;

        // update acceleration

        ax = -omega_x*omega_x*x - gamma*vx + noise_force_x;
        ay = -omega_y*omega_y*y - gamma*vy + noise_force_y;
        // update velocity 

        vx = vx + 0.5*dt*ax;
        vy = vy + 0.5*dt*ay;


        if (std::abs(x) > 1.0){
            std::cout << "Particle lost" << std::endl;
            break;
        }
        if ( i % steps_per_save == 0 ) {
            myfile << t << "," << x << "," << y << "," << vx << "," << vy << std::endl;
        }

        t = t+dt;
    }
    myfile.close();
    return 0;
}