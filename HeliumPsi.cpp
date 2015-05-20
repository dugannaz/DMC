#include "HeliumPsi.hpp"

using namespace std;

HeliumPsi::HeliumPsi() {

}

double HeliumPsi::psiA(Walker &w) {
	double ara;
	double a0 = 0.5;

	double r1 = 0.0;
	for (int d=0; d < w.dim; d++)
		r1 += w.par[0].pos[d] * w.par[0].pos[d];
	r1 = sqrt(r1);

	double r2 = 0.0;
	for (int d=0; d < w.dim; d++)
		r2 += w.par[1].pos[d] * w.par[1].pos[d];
	r2 = sqrt(r2);

	/*** exact node r1=r2 ***/
	ara = -(exp(-r1/a0) * (2.0-r2/a0) * exp(-r2/2.0/a0)-exp(-r2/a0) * (2.0-r1/a0) * exp(-r1/2.0/a0));
	/***********************/

	/*** wrong node ***/
	//double a=0.65;
	//double a=0.35;
	//ara = -(exp(-2.0*r1-a*r2)*(1.0-a*r2)-exp(-2.0*r2-a*r1)*(1.0-a*r1));
	/******************/	

	return ara;
}

double HeliumPsi::psiS(Walker &w) {
       
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
void HeliumPsi::antisymmetry(Walker &w) {

        if (psiA(w) < 0.0) {
                w.fitness = -1;
        }
}

void HeliumPsi::symmetry(Walker &w) {

        if (psiA(w) < 0.0) {
                for (int d=0; d < w.dim; d++) {
                        double temp = w.par[0].pos[d];
                        w.par[0].pos[d] = w.par[1].pos[d];
                        w.par[1].pos[d] = temp;
                }
        }
}

/*void HeliumPsi::symmetry(Walker &w) {

         double r1=0.0;
         for (int d=0;d<w.dim;d++)
               r1+= w.par[0].pos[d]*w.par[0].pos[d];
         double r2=0.0;
         for (int d=0;d<w.dim;d++)
               r2+= w.par[1].pos[d]*w.par[1].pos[d];

        if (r1 > r2) {
                for (int d=0; d < w.dim; d++) {
                        double temp = w.par[0].pos[d];
                        w.par[0].pos[d] = w.par[1].pos[d];
                        w.par[1].pos[d] = temp;
                }

        }
}*/


