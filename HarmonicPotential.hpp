#include "Potential.hpp"
#include "Walker.hpp"

class HarmonicPotential : public Potential {
	public:
		HarmonicPotential(const int dim);
		
		virtual double V(Walker &walker);
		
	private:	
		int dim; //dimensionality of space 
};
