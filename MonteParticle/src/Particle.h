#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "Interaction.h"

// This is a typedef for a random number generator.
// Try boost::mt19937 or boost::ecuyer1988 instead of boost::minstd_rand
//typedef boost::mt19937 base_generator_type;

typedef boost::numeric::ublas::bounded_vector<double, 3> Vector;
const double pi=3.14159265358979323846;

//typeof boost::numeric::ublas::zero_vector<double> (3) Zero;
//int GlobalCounter=0;

using namespace std;
      namespace ublas=boost::numeric::ublas;
  
class particle{
      protected:
    	  	  Vector tempPosition;  // to keep track of future position, modified with reflection or transmission
    	      Vector Phat; // a unit vector in the direction of the momentum
              std::vector<double> E; //Energy
              std::vector<Vector> Position; // Position of particle in Angstroms, +z is down from metal vacuum boundary
              std::vector<string> interaction;
              //unsigned long long OID,ID; // ID of the originating particle and this particle
              double OE; // Energy of the particle that created this particle
              double theta, phi; // direction of motion, Physics spherical coordinate standard see Boas page 261.  theta from +z phi from +x
              double phase; // wave phase in radians, may be used later for diffraction effects.
              std::string particle_type; // electron for most every thing
      
      public: 
              unsigned long long OID,ID; // ID of the originating particle and this particle

              // constructors  
              particle(const double& En,const Vector& R=boost::numeric::ublas::zero_vector<double> (3),const double& T=0, const double& P=0,string Int_type="Pr",string p_type="el");
              particle(const particle oldparticle,const double& En,const Vector& M=boost::numeric::ublas::zero_vector<double> (3),string Int_type="Pr", string p_type="el");

              // Default destructor and copy constructor
              
              void save2file(ofstream *file_ptr,ofstream *file_ptr2, int last);
              double getE(void); // returns the latest Energy of the particle.
              double gettheta(void);
              Vector getPosition(void);
              Vector getDirection(void);
              void move(double r);  // moves the particle a distance r in the direction of motion
              void Delta_direction(const double& D_theta, const double& D_phi); // changes the direction by adding the the amounts D_theta and D_phi to values of theta and phi
              void New_direction(const double& N_theta, const double& N_phi); // replaces old values of theta and phi with passed values
              void Reflected(void);  // this function is not currently fully operational.  It dosn't use input, it simply assumes that the plane of reflection is the xy plane, ie the only surface is the one surface between +z and -z
              void BoundaryCross(double V1,double V2); // Modifies the particle when it crosses a potential boundary from V1 to V2.  Like Reflected it only assumes a xy plane as the boundary
              // the modify functions are used for changing the properties of this particle when it has an interaction
              void modify(const double& Delta_E,const Vector& M=boost::numeric::ublas::zero_vector<double> (3),string Int_type="inlstc",double Delta_phase=0); // used for inelastic collision when the direction is known.
              void modify(const double& Delta_E,const double& Delta_theta,const double& Delta_phi,string Int_type="inlstc",double Delta_phase=0); // used for inelastic collisions
              void modify(const double& Delta_theta,const double& Delta_phi,string Int_type="elstc",double Delta_phase=0); // used for elastic collisions
              void modify(string Int_type,double Delta_phase,const Vector& M=boost::numeric::ublas::zero_vector<double> (3));
              void modify(); // a dummy function for debugging -- prints info about particle
              //friend class interaction;
                                 
}; // class particle
#endif /*PARTICLE_H_*/
