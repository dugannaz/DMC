#include <iostream>
#include <fstream> 
#include <math.h>
#include "/home/nazim/program/gsl-1.16/gsl/gsl_rng.h"
#include "/home/nazim/program/gsl-1.16/gsl/gsl_randist.h"
#include "CoulombPotential.hpp"
#include "HarmonicPotential.hpp"
#include <string>
#include <time.h>
#include "HeliumPsi.hpp"
#include "HarmonicPsi.hpp"

/*
 * Importance Sampling Diffusion Monte Carlo
 * reads parameters from `dmc.input` file
 */

/* application class for DMC  
 */
class DMC {
	public:
	
	DMC(int spaceDimension, int seed, int nElectron, int nWalker, 
		double dt, int thermalTimeSteps, int timeSteps, double alpha, int outputStep, int nBin, 
		double min, double max, Potential *potential, Psi *psiTrial);

	void dmc();
	void refreshWalkerArray();
    
	Psi *psiTrial;
	Potential *potential; // potential to be used
	Walker *walker; // walker array
	gsl_rng *random; // random number generator
	int nElectron;
	int seed; // seed of rng
	int nWalkerMax; // maximum number of walkers
	int nWalker; // number of walkers
	double dt; // time step
	int thermalTimeSteps; // number of timesteps for thermalization 
	int timeSteps; // number of timesteps for data collection
	double alpha; // parameter used in reference energy adjustment
	double referenceE; // reference energy used in DMC
	double timeSumRefE;
	double varRefE;
	double varE;
	int outputStep; // number of time steps at which output is written
	int nBin; // number of bins in the walker histogram
	double min; // coordinate minimum for wavefunction
	double max; // coordinate maximum for wavefunction
	double dx; // space partition
	double timeSumE; // 
	double sqrtDT;
	double averageE;
	int boyut;
	int acceptance;
	double normalization;
	double startTime;
	double stopTime;
	int spaceDimension;
	bool IS;
	bool AS;
};

using namespace std;

/* constructor of application class DMC
 */
DMC::DMC(int spaceDimension, int seed, int nElectron, int nWalker, 
		double dt, int thermalTimeSteps, int timeSteps, double alpha, int outputStep, int nBin, 
		double min, double max, Potential *potential, Psi *psiTrial) {

	this->spaceDimension = spaceDimension;
	this->seed = seed;
	this->dt = dt;
	this->thermalTimeSteps = thermalTimeSteps;
	this->timeSteps = timeSteps;
	this->alpha = alpha;
	this->nWalkerMax = nWalker* 4;
	this->nWalker = nWalker;
	this->outputStep = outputStep;
	this->min = min;
	this->max = max;
	this->normalization = nWalker;
	this->potential = potential;
	this->psiTrial = psiTrial;
	this->AS = false;

	boyut= nElectron*spaceDimension;

	// initialize random number generator

	random = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(random, this->seed);

	/** initialize walkers **/

	walker = new Walker[nWalkerMax+1];

	for (int w = 0; w < nWalkerMax +1 ; w++) {
		walker[w].customize(nElectron, random, max - min, spaceDimension, psiTrial, potential);
		walker[w].oldWalker = new Walker(nElectron, random, max - min, spaceDimension, 
					psiTrial, potential);
		if (psiTrial != NULL)
			psiTrial->symmetry(walker[w]);

	}


	for (int w=0; w < nWalkerMax+1; w++) {
		potential->V(walker[w]);
		walker[w].oldEnergy = walker[w].energy;
		walker[w].fitness = 0;
	}

	dx = (max - min) / double(nBin);

	sqrtDT = sqrt(dt / 1.0);
	
}

/* main calculation loop of DMC implementation
 */
void DMC::dmc() {

	cout << endl << "Diffusion Monte Carlo" << endl << endl;

	/** initialize variables **/
	
	referenceE = 1.0;
	timeSumE = 0.0;
	timeSumRefE = 0.0;
	varRefE = 0.0;
	varE = 0.0;
	int o = 0;

	double timeSumNWalker = 0.0;

	// VMC
	if (true) {
		cout << endl << "VMC initialization ..." << endl << endl;
	
		#pragma omp parallel default(shared)
        	{
			#pragma omp for
			for (int w=0; w < nWalker; w++) 
				for (int t=0; t < 2000; t++) 
					walker[w].metropolis(&sqrtDT);

			if (psiTrial != NULL) {
				#pragma omp for
				for (int w=0; w < nWalker; w++) 
					psiTrial->symmetry(walker[w]);
			}
		}
	}


	cout << "starting DMC ..." << endl;
	
	/** start calculation **/

	for (int t=0; t < timeSteps + thermalTimeSteps; t++) {

		o++; // screen output counter

		/** energy calculation **/
	
		double walkerEnergy = 0.0;
		#pragma omp parallel default(shared)
                {
			double walkerEnergyPart = 0.0;	
			#pragma omp for
			for (int w=0; w < nWalker; w++) {
				walkerEnergyPart += walker[w].energy;
			}

			#pragma omp critic
			walkerEnergy += walkerEnergyPart;
		}

		averageE = walkerEnergy / double(nWalker);		

		/** diffusion step **/

		bool refresh = false;
                #pragma omp parallel default(shared)
		{
			if (IS) {
				#pragma omp for
				for (int w=0; w < nWalker; w++) { 
					walker[w].guide(dt);
				}
			} else {
				#pragma omp for
                                for (int w=0; w < nWalker; w++)
                                        walker[w].diffusion(&sqrtDT);

				if (AS) {
					refresh = true;
					#pragma omp for
					for (int w=0; w < nWalker; w++)
                                        	psiTrial->antisymmetry(walker[w]);
				}
			}
	
		}

		if (refresh)
			refreshWalkerArray();

		/** adjust reference energy **/

		referenceE = averageE + alpha * (1.0 - double(nWalker)  / normalization ) / dt;

                #pragma omp parallel default(shared)
		{
			#pragma omp for
			for (int w=0; w < nWalker; w++) {
				walker[w].branching(referenceE, &dt); 	
			}
		}

		refreshWalkerArray();	

		/** thermalization, data collection and screen output **/

		if (t >= thermalTimeSteps) {

			if (t == thermalTimeSteps)
				startTime = double(time(NULL));

			timeSumE += averageE;
			timeSumRefE += referenceE;
			varRefE += referenceE*referenceE;
			varE += averageE*averageE;
			timeSumNWalker += double(nWalker);

			if (o == outputStep) {
	
				cout << "t = " << t << " N = " << nWalker << " --- <ER>t = " << 
				     timeSumRefE/double(t - thermalTimeSteps +1) 
				     <<  " <E>t = " << timeSumE / double(t - thermalTimeSteps +1);

				cout << endl;		
	
				o = 0;
			}
		} else {
			if (o == outputStep) {
				cout << "t = " << t << " N = " << nWalker << " ER = " << referenceE <<  
				" E = " << averageE << endl;
				o = 0;

			}
		}
	}
	
	/** final output **/
	
	stopTime = double(time(NULL));

	cout << "Computation time : " << double(stopTime - startTime) / 60.0 << " minutes" << endl;
}

/* Applies changes to the walker array according to branching
 */
void DMC::refreshWalkerArray() {

	// death
	int nWalkerNew = 0;
   	int death = 0;
    	int dummy = nWalker-1;
   
    	for (int w = 0; w < nWalker; w++) {
        	if (walker[w].fitness > -1)            
            		nWalkerNew++;
         	else
            		death++;
    	}

    	for (int w = 0; w < nWalkerNew; w++) {
        	if (walker[w].fitness <= -1) {            
        		while(walker[dummy].fitness <= -1)
        			dummy--;
           		walker[w].equal(walker[dummy]); 
		 	dummy--;	
		 }
     	}
     
    	nWalker = nWalkerNew;

	// birth
	for (int w = 0; w < nWalker; w++) {
		
		while(walker[w].fitness >0)	{	
			if (nWalkerNew < nWalkerMax) {
				walker[nWalkerNew].equal(walker[w]);		
				walker[nWalkerNew].fitness=0;	
				nWalkerNew++;			
			}

			walker[w].fitness --;
		}
		
	}
	nWalker = nWalkerNew;
}

/* main function of c++ application.
 * constructs an instance of DMC and starts calculation.
 */
int main() {

	//double startTime = double(time(NULL));

	int spaceDimension, nNucleus, nElectron, nWalker, seed, thermalTimeSteps, 
	     timeSteps, outputStep, nBin, IS, AS;
	double dt, alpha, min, max;
	char *bekle;
	string name;

	ifstream fd("dmc.input");
	
	fd >> name;fd >> spaceDimension;
	fd >> name;fd >> nNucleus;
	fd >> name;fd >> nElectron;
	fd >> name;fd >> nWalker;
	fd >> name;fd >> seed;
	fd >> name;fd >> dt;
	fd >> name;fd >> thermalTimeSteps;
	fd >> name;fd >> timeSteps;
	fd >> name;fd >> alpha;
	fd >> name;fd >> outputStep;
	fd >> name;fd >> nBin;
	fd >> name;fd >> min;
	fd >> name;fd >> max;
	fd >> name;fd >> IS;
	fd >> name;fd >> AS;

	fd.close();
	
	// rng seed from clock
	if (seed == 0)
		seed = time(NULL);

	/*** set potential function and trial wavefunction ***/

	Potential *potential = new CoulombPotential(spaceDimension, nNucleus);
	//Potential *potential = new HarmonicPotential(spaceDimension);

	Psi *psiTrial = new HeliumPsi();
	//Psi *psiTrial = new HarmonicPsi();
	//Psi *psiTrial = NULL;

	/*****************************************************/

	DMC &dmc = *(new DMC(spaceDimension, seed, nElectron, nWalker, dt, 
			thermalTimeSteps, timeSteps, alpha, outputStep, nBin, min, max,
			potential, psiTrial));

	dmc.IS = IS;

	if (dmc.IS == false ) 
		dmc.AS = AS;

	if ((dmc.IS || dmc.AS) && psiTrial==NULL) {
		cout << "Trial Wavefunction needed !!!" << endl;
		return 0;
	}

	dmc.dmc();

	return 0;
}



