#include "Potential.hpp"
#include "Walker.hpp"

class CoulombPotential : public Potential {
	public:
		CoulombPotential(const int dim, int nNucleus);
		
		virtual double V(Walker &walker);
		
	private:	
		int dim; //dimensionality of space 
		int nNucleus; //number of nucleus
};
