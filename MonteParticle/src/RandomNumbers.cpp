#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <math.h>
/*
#include <boost/random.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
*/
#include "RandomNumbers.h"
#include "interpolation.h" // see http://www.alglib.net/
#include "stdafx.h"

#define Pi 3.14159265358979323846
// This is a typedef for a random number generator.
typedef boost::mt19937 base_generator_type;  // Could use boost::rand48 instead it is about 2x faster but only 2^48 before it repeats, see http://www.boost.org/doc/libs/1_43_0/doc/html/boost_random/reference.html#boost_random.reference.concepts.uniform_random_number_generator

//$$$ need to check for divide by zeros and sqrt (-1)

//constructor
RandNum::RandNum(unsigned long s=33u) : seed(s){
	gen.seed(s);
};

double RandNum::Triangular(double A,double B, double C){ //Triangular Distribution Function
	if (A>B){double temp=B; B=A;A=temp;}
	if (C>B){double temp=B; B=C;C=temp;}
	if (A>C){double temp=C; C=A;A=temp;}

	//boost::uniform_real<> R_dist(0,1); // range
	//boost::variate_generator<base_generator_type&, boost::uniform_real<> > R(gen, R_dist);
	double qMode=(C-A)/(B-A);
	double r=R();
    if (qMode>=r){
    	return A+sqrt(r*(B-A)*(C-A));
    }else{
    	return B-sqrt((1-r)*(B-A)*(B-C));
  	}
};
double RandNum::Logistic(double A,double B, double C){ // Logistic Distribution Function
	double r=R();
	return A+B*log(r/(1-r));
};
double RandNum::Weibull(double A,double B, double C){ // Weibull Distribution Function. (0,1,2) is the Rayleigh distribution
	double r=R();
	return A+B*pow(-log(r),1/C);
};
double RandNum::Lorentz(double A,double B){ // Lorentz or Cauchy Distribution Function.  See page A-19 and A-20 in Regress+ Appendix A A Compendium of Common Probability Distributions Version 2.3
	double r=R();
	return A+B*tan(Pi*(r-0.5));
};
double RandNum::R(double xmin,double xmax){ // a uniform distribution between xmin and xmax.
	boost::uniform_real<> R_dist(xmin,xmax); // range
	boost::variate_generator<base_generator_type&, boost::uniform_real<> > rnum(gen, R_dist);
	return rnum();
};
double RandNum::Normal(double avg,double sigma){ // a Normal or Gaussian distribution with average and sigma.
   boost::normal_distribution<> G_dist(avg,sigma); // first is average then sigma
   boost::variate_generator<base_generator_type&, boost::normal_distribution<> > G(gen, G_dist);
   return G();
};

// RSM uses the Rejection Sampling Method to return a random number that will satisfy the probability distribution function (pdf).
// no effort is made to make this particular efficient i.e. a simple box (y=0..Ymax x=xmin .. xmax) is used for the area around the pdf function.
// Also no effort is made to check for errors, such as the possibility that the box doesn't encompass the function --
// !!!!!!!!!!!!!  care should be taken to avoid such issues when using RSM    !!!!!!!!!!!
// Parameters are (A,B,C,function(A,B,C,x),Xmin,Xmax,Ymax)
// The function RSM uses Polymorphism so that it will work for pdf functions with different number of parameters: from 0 up to 3.
// Also the last one uses the same method but a cubic spline from data as its function
double RandNum::RSM(double A,double B,double C, double (*pdf)(double, double, double, double),double Xmin,double Xmax, double Ymax  ){
	double ry=R(0,Ymax);
	double rx=R(Xmin,Xmax);
	while (ry>pdf(A,B,C,rx)){ // this could be very slow -- even be an infinite loop.
		ry=R(0,Ymax);
		rx=R(Xmin,Xmax);
	}
	return rx;
}
double RandNum::RSM(double A,double B, double (*pdf)(double, double, double),double Xmin,double Xmax, double Ymax  ){
	double ry=R(0,Ymax);
	double rx=R(Xmin,Xmax);
	while (ry>pdf(A,B,rx)){ // this could be very slow -- even be an infinite loop.
		ry=R(0,Ymax);
		rx=R(Xmin,Xmax);
	}
	return rx;
}
double RandNum::RSM(double A, double (*pdf)(double, double),double Xmin,double Xmax, double Ymax  ){
	double ry=R(0,Ymax);
	double rx=R(Xmin,Xmax);
	while (ry>pdf(A,rx)){ // this could be very slow -- even be an infinite loop.
		ry=R(0,Ymax);
		rx=R(Xmin,Xmax);
	}
	return rx;
}
double RandNum::RSM( double (*pdf)(double),double Xmin,double Xmax, double Ymax  ){
	double ry=R(0,Ymax);
	double rx=R(Xmin,Xmax);
	while (ry>pdf(rx)){ // this could be very slow -- even be an infinite loop.
		ry=R(0,Ymax);
		rx=R(Xmin,Xmax);
	}
	return rx;
}
double RandNum::RSM(const alglib::spline1dinterpolant& S, double Xmin,double Xmax, double Ymax  ){
	double ry=R(0,Ymax);
	double rx=R(Xmin,Xmax);
	while (ry>spline1dcalc(S,rx)){ // this could be very slow -- even be an infinite loop.
		ry=R(0,Ymax);
		rx=R(Xmin,Xmax);
	}
	return rx;
}
double RandNum::RSMx(const alglib::spline2dinterpolant& S, double Xmin,double Xmax, double Fmax  , double y){
	double rf=R(0,Fmax);
	double rx=R(Xmin,Xmax);
	double f=spline2dcalc(S,y,rx);
	//std::cout<<"RSMx:  rx= "<<rx<<"  y= "<<y<<"  f= "<<f<<" fmax= "<<Fmax<<std::endl;
	while (rf>f){ // this could be very slow -- even be an infinite loop.
		rf=R(0,Fmax);
		rx=R(Xmin,Xmax);
		f=spline2dcalc(S,y,rx);
		//std::cout<<"RSMx:  rx= "<<rx<<"  y= "<<y<<"  f= "<<f<<" fmax= "<<Fmax<<std::endl;
	}
	return rx;
}
double RandNum::RSMy(const alglib::spline2dinterpolant& S, double Xmin,double Xmax, double Fmax  , double y){
	double rf=R(0,Fmax);
	double rx=R(Xmin,Xmax);
	double f=spline2dcalc(S,rx,y);
	//std::cout<<"RSMy: rx= "<<rx<<"  y= "<<y<<"  f= "<<f<<" fmax= "<<Fmax<<std::endl;
	while (rf>f){ 
		rf=R(0,Fmax);
		rx=R(Xmin,Xmax); 
		f=spline2dcalc(S,rx,y);
		//std::cout<<"RSMy: rx= "<<rx<<"  y= "<<y<<"  f= "<<f<<" fmax= "<<Fmax<<std::endl;
	}
	return rx;
}
alglib::real_1d_array RandNum::RSM(const alglib::spline2dinterpolant& S, double Xmin,double Xmax,double Ymin,double Ymax, double Fmax){
	double rf=R(0,Fmax);
	double rx=R(Xmin,Xmax);
	double ry=R(Ymin,Ymax);
	double f=spline2dcalc(S,rx,ry);
	while (rf>f){ // this could be very slow -- even be an infinite loop.
		rf=R(0,Fmax);
		rx=R(Xmin,Xmax);
		ry=R(Ymin,Ymax);
	}
	double _r[]={rx,ry};
	alglib::real_1d_array r;
	r.setcontent(2,_r);
	return r;
}

// Below are the functions that make it possible to have a pdf from data
alglib::real_1d_array RandNum::stdvec2array(std::vector<double>& x){
	alglib::real_1d_array x_a;
	x_a.setlength(x.size());
	for (unsigned int i=0;i<x.size();i++){
		x_a(i)=x[i];
	}
	return x_a;
}
alglib::spline1dinterpolant RandNum::spline(alglib::real_1d_array& x,alglib::real_1d_array& y){
    alglib::spline1dinterpolant s;
    alglib::spline1dbuildcubic(x, y, s);
    return s;
}
alglib::spline1dinterpolant RandNum::spline(std::vector<double>& x,std::vector<double>& y){
	alglib::spline1dinterpolant s;
	alglib::real_1d_array x_a;
	alglib::real_1d_array y_a;
	x_a=stdvec2array(x);
	y_a=stdvec2array(y);
    alglib::spline1dbuildcubic(x_a,y_a, s);
    return s;
}
alglib::spline2dinterpolant RandNum::spline2d(alglib::real_1d_array& x,alglib::real_1d_array& y, alglib::real_1d_array& f){
    alglib::spline2dinterpolant s;
    int n=x.length();
    int m=y.length();
    int p=f.length();

    if (n*m==p){
    	alglib::spline2dbuildbicubicv(x,n,y,m,f,1,s);
    }else {std::cout<<"spline2d didn't work"<<std::endl;  return s;}
    return s;
}
alglib::spline2dinterpolant RandNum::spline2d(std::vector<double>& x,std::vector<double>& y, std::vector<double>& f){
    alglib::spline2dinterpolant s;
	alglib::real_1d_array x_a="[0.0, 100 , 1000.0]";
	alglib::real_1d_array y_a= "[0.0, 180.0]";
	alglib::real_1d_array f_a ="[303.00,302.00,301.00,300.30,300.20,300.10]"; // [f(x[0],y[0])  ,  f(x[1],y[0])  , f(x[2],y[0])  , f(x[0],y[1])  , .....        ]
	/*for (int i=0;i<5;i++){
		//std::cout<<i<<" x= "<<x[i]<<std::endl;
		std::cout<<i<<" y= "<<y[i]<<std::endl;
	}*/

	x_a=stdvec2array(x);
	y_a=stdvec2array(y);
	f_a=stdvec2array(f);



    int n=x_a.length();
    int m=y_a.length();
    int p=f_a.length();
    //std::cout<<"size of E="<<n<<"  size of theta="<<m<<" size of f(x,y)="<<p<<std::endl;
    if (n*m==p){
    	//std::cout<<"before "<<std::endl;
    	alglib::spline2dbuildbicubicv(x_a,n,y_a,m,f_a,1,s);
    	//std::cout<<"after "<<std::endl;
    }else {
    	std::cout<<"sline2d didn't work"<<std::endl;
    	std::cout<<"size of E="<<n<<"  size of theta="<<m<<" size of f(x,y)="<<p<<std::endl;
    	return s;}
    return s;
}
// end function for pdf from data

int RandNum::HistoGram(double minX, double maxX,unsigned int numbins, char *file_name, std::vector<double>& Data){
    std::ofstream datafile;
    datafile.open (file_name);
    datafile<<"##,##"<<std::endl;
    datafile<<"#@Histogram of probability distribution function (pdf)"<<std::endl;
    datafile<<"#bins counts"<<std::endl;

    numbins=numbins++; // this is so there is one more bin than needed to put the rest of the data in.  It is then skipped at the end.
    if (maxX<minX){double temp=maxX;maxX=minX;minX=temp;}

    std::vector<double> bins;
    std::vector<int> counts;

    bins.reserve(numbins+4);
    counts.reserve(numbins+4);

    double stepsize=1.;
    stepsize=(maxX-minX)/numbins;

    // first set the values for each bin
    for (unsigned int i=0;i<=numbins;i++){
    	bins.push_back(minX+i*stepsize);
    	counts.push_back(0);
    }

    for (unsigned int i=0;i<Data.size();i++){
    	double d=Data.at(i);
    	unsigned int j=0;
    	while (minX+j*stepsize<d && j<numbins){
    		j++;
    	}

    	if(j<numbins){counts[j]++;;}

    }

    for (unsigned int i=1;i<=numbins-1;i++){
    	datafile<<(bins.at(i)-stepsize/2)<<" "<<counts.at(i)<<std::endl;
    }

    datafile.close();
    bins.clear();
    counts.clear();
    return 0;
}; // end of function histogram


