#include <gsl/gsl_rng.h>
#include <string>
#include <fstream>
#include "Psi.hpp"
#include "Potential.hpp"

#ifndef _Walker_HEADER_  
#define _Walker_HEADER_

/* Walker header file
 */
using namespace std;

/*
 * Class definition for Particle
 * Holds data related to a particle
 */
class Particle {
	public:
		void equal(Particle const &other, int dim);
		double *pos;
};

//class Psi;
//class Potential;

/*
 * Class definiton for Walker
 * Holds varius data for a system of particles
 * Have functions for input output and manipulation and MPI
 */
class Walker {
	public:
		Walker (); 
		Walker(int nElectron, gsl_rng *random, double spaceSize, int dim, 
			Psi *psiTrial, Potential *potential);
		~Walker();
		void customize(int nElectron, gsl_rng *random, double spaceSize, int dim,
                        Psi *psiTrial, Potential *potential);
		void equal(Walker const &other);
		void zero();
		double distance(int const p1, int const p2);
		double distance(int const p1);
		double distance(Walker &w);
		double distance();
		void vDrift(double *array,double step, double psiG);
		double localEnergy(double step);
		void diffusion(double *sqrtDT);
		void metropolis(double *sqrtDT);
		void guide(double dt);
		void branching(double referenceE, double *dt);	
	
		Particle *par; // particle array
		int size; // number of particles
		double energy; // energy of the walker
		double oldEnergy; // energy of the previous step
		int fitness; // energy related fitness of walker
		int index; // index of walker in the array
		int dim; //dimensionality of space
		int boyut;
		Psi *psiTrial;
		gsl_rng *random;
		Potential *potential;
		Walker *oldWalker;
	
	private:	
		void particleInitialize(gsl_rng *random,double spaceSize);
			
};

//int Walker::dim;
//int Walker::boyut;
//Psi Walker::psi;

#endif

