#include <math.h>

#ifndef _Psi_HEADER_  
#define _Psi_HEADER_

class Walker;

class Psi {
	public:
		virtual double psiA(Walker &w) = 0;
    		virtual double psiS(Walker &w) = 0;
		virtual void antisymmetry(Walker &w) = 0;
		virtual void symmetry(Walker &w) = 0;
   	
};

#endif
