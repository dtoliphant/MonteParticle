#ifndef INTERACTION_H_
#define INTERACTION_H_

#include "Particle.h"
#include "RandomNumbers.h"
#include "interpolation.h" // see http://www.alglib.net/
#include "stdafx.h"

class interaction{
	protected:
		unsigned long int current,size; // the most recently used index for E, SIGMA, etc likely close to the next use so use it to start next search
		double LAMBDA,SIGMA,Energy; // current lambda, SIGMA and Energy

		RandNum RandGenerator;  //argument is seed
		alglib::spline1dinterpolant sMax1d;   	// spine interpolant for sigma the max value of Sigma(theta,Energ)  to use with RMS to box in function
		alglib::spline1dinterpolant s1d;   	// spine interpolant for sigma total for a given energy
		alglib::spline2dinterpolant s2d;	// spine interpolant for sigma  for a given energy and Angle

		string particle_type; // electron for most things
		string material; // Gold (Au) for our case
		string interaction_type; // describes the type of interaction: elastic, inelastic, plasmon, Auger, etc.
		double n_density; // the number density used to calculate mean free path from sigma in units of number per cubic Angstroms.

	    double str2double(char*);
	    void linepars(string ,std::vector<double>&,char delim =','); // function to turn a line of data from a file to a array of doubles char is deliminator default is ,
	    
	public:
	    // constructors
		interaction(char *file_name, unsigned long seed); //for when the interaction data is in a file
		interaction(char *file_name, int numAngles, unsigned long seed); //for when the interaction data is in a file with angular info also
		interaction(); // for when SIGMA is constant

		double sigma(double energy,double angle); // Gives the value of sigma for a given energy and angle.
		double sigma(double energy); // returns the total cross section for given energy, also interpolates between values of energy.
		double sigma(); // returns the latest value/current value of sigma
		double MaxSig(double energy); // returns the max value of sigma(theta). used for getting random number representing sigma(theta) distribution.
		double lambda(double energy);   // returns the mean free path for given energy, also interpolates between values of energy
		//double lambda(particle&);
		double lambda();
		static double H(double, double, double);
		void process(std::vector<particle>&); // for inelastic interactions changes the properties of the particle according to the interaction
		void process(particle&);  // for elastic interactions changes the properties of the particle according to the interaction
};
#endif /*INTERACTION_H_*/
