#include "Psi.hpp"
#include "Walker.hpp"

class HeliumPsi : public Psi {
	public:
		HeliumPsi();

		virtual double psiA(Walker &w);
    		virtual double psiS(Walker &w);
		virtual void antisymmetry(Walker &w);
		virtual void symmetry(Walker &w);
   	
};

