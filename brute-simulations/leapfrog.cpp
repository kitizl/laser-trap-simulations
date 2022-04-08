#include <stdio.h>
#include <iostream>

#include <numeric> // to get std::iota from
#include <vector>
#include <cmath>
#include <random>

#include <fstream>

int main(){
    // We're just gonna make a leapfrog algorithm, fuck the BAOAB scheme
    std::ofstream myfile;
	myfile.open("electric-driving-simulations/pfc-attempt-exp.csv");
	myfile << "t,x,v\n";
    const float PI = 3.14159265358;
    float max_time = 1.0; // 1 second
    double dt = 1e-6; // 1 nanosecond sampling time
    long Nsteps = 1e6; // one billion samples

    // initializing the position and velocity arrays
    float x = 0;
    float v = 0;
    float a = 0;
    double t = 0;

    int steps_per_save = 1;
    // other constants

    float q = -7*1.6e-19; // -7e charges
    float E = 10000; // 5000 V/m

    // parameters
    float density = 1860; // kg/m^3
    float radius = 71.5e-9;
    float kBT = 1.38e-23 * (273+22); // joules
    float gamma = 32602.0; // Hz
    float omega = 2*PI*(40e3);
    float mass = (4/3)*PI*radius*radius*radius*density;

    std::default_random_engine generator;
  	std::normal_distribution<double> R(0.0,1.0);

    float electric_force = q*E/mass;
    float noise_force = 0.0;

    float mod_depth = 0.04;
    float mod_freq = 2 * omega;

    for(int i = 1; i < Nsteps-1; i++){
        noise_force = R(generator)*std::sqrt(2*kBT*gamma/mass);//
        electric_force = 0;//q*E*v/mass;
        a = -omega*omega*(1+mod_depth*x*v)*x - gamma*v + noise_force + electric_force;
        v = v + 0.5*dt*a;
        x = x + dt*v;
        a = -omega*omega*(mod_depth*x*v)*x - gamma*v + noise_force + electric_force;
        v = v + 0.5*dt*a;

        if (std::abs(x) > 1.0){
            std::cout << "Particle lost" << std::endl;
            break;
        }
        if ( i % steps_per_save == 0 ) {
            std::cout<< i << "," << t << "," << x << "," << v << std::endl;
            myfile << t << "," << x << "," << v << std::endl;
        }

        t = t+dt;
    }
    myfile.close();
    return 0;
}