#ifndef RANDOMNUMBERS_H_
#define RANDOMNUMBERS_H_
#include <boost/random.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "interpolation.h" // see http://www.alglib.net/
#include "stdafx.h"

typedef boost::mt19937 base_generator_type;  // Could use boost::rand48 instead it is about 2x faster but only 2^48 before it repeats, see http://www.boost.org/doc/libs/1_43_0/doc/html/boost_random/reference.html#boost_random.reference.concepts.uniform_random_number_generator

//base_generator_type gen; // seed
//boost::uniform_real<> R_dist(0,1); // range
//boost::variate_generator<base_generator_type&, boost::uniform_real<> > R(gen, R_dist);


class RandNum{
 private:
	unsigned long seed;
	base_generator_type gen; //

	boost::uniform_real<> R_dist(double,double); // range
    boost::normal_distribution<> G_dist(double,double); // first is average then sigma

	boost::variate_generator<base_generator_type&, boost::uniform_real<> > R(base_generator_type gen, boost::uniform_real<> R_dist);
	boost::variate_generator<base_generator_type&, boost::normal_distribution<> > G(base_generator_type, boost::normal_distribution<> );

 public:
	// constructors
	RandNum(unsigned long); // constructor
	RandNum(); // constructor
	// Using default copy and destructor functions
	int RandNum::HistoGram(double, double,unsigned int, char* ,std::vector<double>& Data); // a simple histogram function outputs histogram to "DataFile.dat"


	double Triangular(double,double, double); //Triangular Distribution Function
	double Logistic(double,double, double);  // Logistic Distribution Function
	double Lorentz(double, double); // Lorentz or Cauchy Distribution Function.  See page A-19 and A-20 in Regress+ Appendix A A Compendium of Common Probability Distributions Version 2.3
	double Weibull(double,double, double);  // Weibull Distribution Function. (0,1,2) is the Rayleigh distribution
	double R(double xmin=0,double xmax=1); // a uniform distribution between xmin and xmax.
	double Normal(double,double); // a Normal or Gaussian distribution with average and sigma.


	// RSM uses the Rejection Sampling Method to return a random number that will satisfy the probability distribution function (pdf).
	// no effort is made to make this particular efficient i.e. a simple box (y=0..Ymax x=xmin .. xmax) is used for the area around the pdf function.
	// Also no effort is made to check for errors, such as the possibility that the box doesn't encompass the function --
	// !!!!!!!!!!!!!  care should be taken to avoid such issues when using RSM    !!!!!!!!!!!
	// Parameters are (A,B,C,function(A,B,C,x),Xmin,Xmax,Ymax)
	// The function RSM uses Polymorphism so that it will work for pdf functions with different number of parameters: from 0 up to 3.
	// Also the last one uses the same method but a cubic spline from data as its function
	double RSM(double,double,double, 	    double (*pdf)(double,double,double,double),	double,double,double  ); // (A,B,C,function(A,B,C,x),Xmin,Xmax,Ymax)
	double RSM(double,double, 			    double (*pdf)(double,double,double),		double,double,double  ); // (A,B,  function(A,B,  x),Xmin,Xmax,Ymax)
	double RSM(double, 					    double (*pdf)(double,double),				double,double,double  ); // (A,,   function(A,    x),Xmin,Xmax,Ymax)
	double RSM(		 					    double (*pdf)(double),						double,double,double  ); // (      function(      x),Xmin,Xmax,Ymax)
	double RSM(const alglib::spline1dinterpolant& S, double,double, double); 									 // (S						,Xmin,Xmax,Ymax)
	double RSMx(const alglib::spline2dinterpolant& S, double,double, double, double); 							//  (S						,Xmin,Xmax,Ymax, other parameter in 2d besides x)
	double RSMy(const alglib::spline2dinterpolant& S, double,double, double, double); 							//  (S						,Xmin,Xmax,Ymax, other parameter in 2d besides y)
	alglib::real_1d_array RSM(const alglib::spline2dinterpolant& S, double Xmin,double Xmax,double Ymin,double Ymax, double Fmax);

	// These are functions for stitching together data with a cubic spline to use with the last of the RSM functions
	alglib::real_1d_array stdvec2array(std::vector<double>& x);  // Converts a std:vector<double> into an array type used by ALGLIB for interpolation
	alglib::spline1dinterpolant spline(alglib::real_1d_array& x,alglib::real_1d_array& y); // returns the interpolant when given ALGLIB data array.  Used for interpolation
	alglib::spline1dinterpolant spline(std::vector<double>& x,std::vector<double>& y);  // returns the interpolant when given std::vector<double> data type.  Used for interpolation
	alglib::spline2dinterpolant spline2d(alglib::real_1d_array& x,alglib::real_1d_array& y, alglib::real_1d_array& f);
	alglib::spline2dinterpolant spline2d(std::vector<double>& x,std::vector<double>& y, std::vector<double>& f);

};

#endif /* RANDOMNUMBERS_H_ */
