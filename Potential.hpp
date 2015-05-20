#ifndef _Potential_HEADER_  
#define _Potential_HEADER_

class Walker;

using namespace std;

class Potential {
	public:
		virtual double V(Walker &walker) = 0;
		double scale;
};

#endif
