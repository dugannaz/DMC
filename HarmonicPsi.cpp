#include "HarmonicPsi.hpp"
#include <math.h>

//using namespace NEWMAT;

HarmonicPsi::HarmonicPsi() {

}

double HarmonicPsi::psiA(Walker &w) {

	double psiValue;

	double r2=0.0;
     
        for (int p=0; p < w.size; p++)      
           	for (int d=0;d<w.dim;d++)
               		r2+= w.par[p].pos[d]*w.par[p].pos[d];

	// slater determinant

	int n = w.size;

	Matrix &slater = *(new Matrix(n,n));

	// orbital definitions

	for (int p=0; p < n; p++) {
	
		slater(p+1,1) = w.par[p].pos[0];
		slater(p+1,2) = w.par[p].pos[1];
		slater(p+1,3) = w.par[p].pos[2];
		slater(p+1,4) = w.par[p].pos[0]*w.par[p].pos[1];
		slater(p+1,5) = w.par[p].pos[0]*w.par[p].pos[2];
		slater(p+1,6) = w.par[p].pos[1]*w.par[p].pos[2];
		slater(p+1,7) = w.par[p].pos[0]*w.par[p].pos[1]*w.par[p].pos[2];
		slater(p+1,8) = 1.0;			
	}

	// total fermionic function
	
	// -0.0866
        psiValue = exp( (-0.07)*r2) * slater.Determinant();

	delete &slater;

	return psiValue;
}

double HarmonicPsi::psiS(Walker &w) {
       
	double p_a = 0.5;
      	double p_c = 0.0;
	
     	double psiT= psiA(w);
    	double r2=0.0;
     
      	for (int p=0;p< w.size;p++)      
           for (int d=0;d<w.dim;d++)
               r2+= w.par[p].pos[d]*w.par[p].pos[d];
              
	return  sqrt(p_c * exp(-p_a*r2) + psiT * psiT);
}

/* Applies fermionic antisymmetry condition
 */
void HarmonicPsi::antisymmetry(Walker &w) {

        if (psiA(w) < 0.0) {
                w.fitness = -1;
        }
}

void HarmonicPsi::symmetry(Walker &w) {

        if (psiA(w) < 0.0) {
                for (int d=0; d < w.dim; d++) {
                        double temp = w.par[0].pos[d];
                        w.par[0].pos[d] = w.par[1].pos[d];
                        w.par[1].pos[d] = temp;
                }
        }
}

