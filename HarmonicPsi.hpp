#include "Psi.hpp"
#include "Walker.hpp"
#include "newmat.h"

class HarmonicPsi : public Psi {
	public:
		HarmonicPsi();

		virtual double psiA(Walker &w);
    		virtual double psiS(Walker &w);
		virtual void antisymmetry(Walker &w);
		virtual void symmetry(Walker &w);
   	
};

