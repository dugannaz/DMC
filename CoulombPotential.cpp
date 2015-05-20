#include "CoulombPotential.hpp"
#include "math.h"

using namespace std;

CoulombPotential::CoulombPotential(const int dim, int nNucleus) {
	
	this->dim = dim;
	this->nNucleus = nNucleus;
	scale = 1.0;
}

double CoulombPotential::V(Walker &w) {
	
	double v = 0.0;
	
	for (int p = 0; p < w.size-1; p++) {
		for (int p1 = p+1; p1 < w.size; p1++) {
			
			double distance = w.distance(p,p1);
				
			if (distance < 0.0001)
				distance = 0.0001;

			v += 1.0/distance;
		}
	}

        for (int p = 0; p < w.size; p++) {

		double distance = w.distance(p);

                if (distance < 0.0001)
                	distance = 0.0001;

		v -= 2.0/distance;
	}

	w.energy = v;
	return v; 
}
