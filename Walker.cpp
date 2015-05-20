#include "Walker.hpp"
#include <iostream>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

/* equalize two particles
 */
void Particle::equal(Particle const &other, int dim) {

	if (this != &other) {
		for (int d = 0; d < dim; d++) {
			this->pos[d] = other.pos[d];
		}
	}
}

/* equalize two walkers
 */
void Walker::equal(Walker const &other) {

	if (this != &other) {
		for (int p = 0; p < other.size; p++)
			this->par[p].equal(other.par[p], dim);
		this->size = other.size;
		this->energy = other.energy;
		this->oldEnergy = other.oldEnergy;
		this->fitness = other.fitness;	
		this->index = other.index;
		this->dim = other.dim;
		this->boyut = other.boyut;
		this->potential = other.potential;
		this->psiTrial = other.psiTrial;
		this->random = other.random;
	}	
}

/* Constructors
 */
Walker::Walker () {
	
}

Walker::Walker (int nElectron, gsl_rng *random, double spaceSize, int dim, 
			Psi *psiTrial, Potential *potential) {

	this->customize(nElectron, random, spaceSize, dim, psiTrial, potential);
}

/* Destructor
 */
Walker::~Walker() {
	
	delete [] par;
}

/* Customize walker according to given parameters
 */
void Walker::customize(int nElectron, gsl_rng *random, double spaceSize, int dim,
                        Psi *psiTrial, Potential *potential) {	

	size = nElectron;

	par = new Particle[size];
	energy = 0.0;
	oldEnergy = 0.0;
	fitness = 0;
	this->dim = dim;
	this->boyut = dim * nElectron;
	this->random = random;
	this->psiTrial = psiTrial;
	this->potential = potential;
	particleInitialize(random, spaceSize);
}

/* initialize particles according to given parameters
 */
void Walker::particleInitialize(gsl_rng *random, double spaceSize) {

	for (int p = 0; p < size; p++) {

		par[p].pos = new double[dim];
		for (int d = 0; d < dim; d++) 
			par[p].pos[d] = 0.0;

		for (int d = 0; d < dim; d++)  {
			if (random != NULL)
				par[p].pos[d] = spaceSize * gsl_rng_uniform(random) 
								   - spaceSize / 2.0;
			else
				par[p].pos[d] = 0.0;
		}	
	}
}

/*
 * Make all the coordinates of all the particles equal to zero
 */
void Walker::zero() {

	for (int p = 0; p < size; p++) 
		for (int d = 0; d < dim; d++) 
			par[p].pos[d] = 0.0;
}

/*
 * Find distance between two particles
 */
double Walker::distance(int const p1, int const p2) {

	double square = 0.0;
	double distance;
	for (int d = 0; d < dim; d++) {
		distance = par[p1].pos[d]-par[p2].pos[d];
		square += distance * distance;
	}
	return sqrt(square);
}

double Walker::distance(int const p1) {

        double square = 0.0;
        double distance;
        for (int d = 0; d < dim; d++) {
                distance = par[p1].pos[d];
                square += distance * distance;
        }
        return sqrt(square);
}


/*
 * Find distance in phase space of current walker with w 
 */
double Walker::distance(Walker &w) {

	double square = 0.0;
	double distance;
	for (int p = 0; p < size; p++) 
		for (int d = 0; d < dim; d++) {
			distance = par[p].pos[d] - w.par[p].pos[d];
			square += distance * distance;
		}

	return sqrt(square);
}

double Walker::distance() {

        double square = 0.0;
        double distance;
        for (int p = 0; p < size; p++)
                for (int d = 0; d < dim; d++) {
                        distance = par[p].pos[d];
                        square += distance * distance;
                }

        return sqrt(square);
}

void Walker::vDrift(double *array,double step, double psiG){

        double ara[boyut+1];
        int parcacik,dimension;

        for (int k=0;k<boyut;k++){
                parcacik=int(k/dim);
                dimension=k-parcacik*dim;
                par[parcacik].pos[dimension]+=step;
                ara[k]= psiTrial->psiA(*this);
                par[parcacik].pos[dimension]-=2.0*step;
                ara[k]-= psiTrial->psiA(*this);
                par[parcacik].pos[dimension]+=step;
                array[k]=ara[k]/(2.0*step*psiG);
        }
}

double Walker::localEnergy(double step){
        double ara;
        double temp, delkare;
        int parcacik,dimension;

        delkare = 0.0;
        temp= psiTrial->psiA(*this);
        for (int k=0;k<boyut;k++){
            parcacik=(k/dim);
            dimension=k-parcacik*dim;
            ara= -2.0 * temp;
            par[parcacik].pos[dimension]+=step;
            ara+= psiTrial->psiA(*this);
            par[parcacik].pos[dimension]-=step;
            par[parcacik].pos[dimension]-=step;
            ara+= psiTrial->psiA(*this);
            par[parcacik].pos[dimension]+=step;
            delkare+=ara;

        }
        delkare/=(step*step); 
        ara=potential->V(*this)-(delkare/2.0)/temp;
        return ara;
}

/* applies diffusion step to the particles of walker
 */
void Walker::diffusion(double *sqrtDT) {

        for (int p = 0; p  < size; p++) {
                for (int d = 0; d < dim; d++) {
                        par[p].pos[d] += gsl_ran_gaussian (random, 1.0) * (*sqrtDT);

                }
        }
        oldEnergy = energy;
        energy=potential->V(*this);
	//cout << energy << endl;
}


void Walker::metropolis(double *sqrtDT) {

	int acceptance = 0;

        for (int n=0; n<1; n++) {
                oldWalker->equal(*this);

                bool accept = false;
                while (!accept) {

                        acceptance++;

                        for (int p = 0; p  < size; p++)
                                for (int d = 0; d < dim; d++)
                                        par[p].pos[d] +=gsl_ran_gaussian (random, 1.0) * (*sqrtDT) * 1.0;

                        double acceptanceProbability = pow(psiTrial->psiA(*this),2.0) 
						     / pow(psiTrial->psiA(*oldWalker),2.0);
                        if (acceptanceProbability > gsl_rng_uniform(random))
                                accept = true;
                        else
                                this->equal(*oldWalker);
                }
        }
}

void Walker::guide(double dt) {

        double temp[boyut+1];
        double temp2[boyut+1];
        double sum,ara;

        bool confirmed=false;

        double psiG = psiTrial->psiA(*this);

        vDrift(temp,0.0001, psiG);
        oldWalker->equal(*this);
        while(!confirmed) {

                confirmed = true;

                for (int p = 0; p  < size; p++)
                        for (int d = 0; d < dim; d++)
                                par[p].pos[d] = oldWalker->par[p].pos[d]+gsl_ran_gaussian(random, 1.0)
                                                * sqrt(dt) + temp[p*dim+d]*dt;

                double psiG1 = psiTrial->psiA(*this);

                if (psiG1 < 0.0) {
                        this->equal(*oldWalker);
                        confirmed = false;
                }

                if (confirmed) {
                        vDrift(temp2,0.0001, psiG1);

                        double vectorSquare = 0.0;
                        for (int p = 0; p  < size; p++)
                                for (int d = 0; d < dim; d++)
                                        vectorSquare += pow(oldWalker->par[p].pos[d] - par[p].pos[d]
                                                        - temp2[p*dim+d]*dt,2.0);

			double vectorSquare1 = 0.0;
                        for (int p = 0; p  < size; p++)
                                for (int d = 0; d < dim; d++)
                                        vectorSquare1 += pow(par[p].pos[d] - oldWalker->par[p].pos[d]
                                                        - temp[p*dim+d]*dt,2.0);

                        double acceptanceProbability = exp(-vectorSquare/2.0/dt) * pow(psiG1,2.0)
                                                      / exp(-vectorSquare1/2.0/dt) / pow(psiG,2.0);

                        if (acceptanceProbability < gsl_rng_uniform(random)) {
                                this->equal(*oldWalker);
                                confirmed = false;
                        }

                }
        }

        oldEnergy = energy;
        energy=localEnergy(0.0001);

}

/* applies branching to the walker
 */
void Walker::branching(double referenceE, double *dt) {

        double weight = exp(-(*dt) *( (oldEnergy + energy) / 2.0 -  referenceE));
        int fitness = int(weight + gsl_rng_uniform(random));

        if (fitness > 2) {
                //cout << "uyari!!! branching " << fitness << endl;
                fitness = 2;
        }
        else if (fitness < 0)
                fitness = 0;

        this->fitness = fitness-1;
}


