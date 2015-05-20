#include "HarmonicPotential.hpp"

using namespace std;

HarmonicPotential::HarmonicPotential(const int dim) {
	
	this->dim = dim;
}

double HarmonicPotential::V(Walker &w) {
	
	double v = 0.0;
	
	for (int p = 0; p < w.size; p++)
		for (int d=0; d < dim; d++)
			v+= 0.015 * w.par[p].pos[d] * w.par[p].pos[d];

	w.energy = v;
	return v; 
}
