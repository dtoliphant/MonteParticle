#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <string>
#include <boost/random.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "Interaction.h"
#include "RandomNumbers.h"
#include "interpolation.h" // see http://www.alglib.net/
#include "stdafx.h"

#define deg2rad 0.017453292519943297
#define Pi 3.14159265358979323846


// Constructors, overloaded for different amounts of data for an interaction.
      interaction::interaction(char *file_name, unsigned long seed)
      	  :current(0), size(0), LAMBDA(1.), SIGMA(1.),Energy(55.),
      	   RandGenerator(seed),sMax1d(),s1d(),s2d(),particle_type("\n"),material("\n"),interaction_type("\n"), n_density(0.1)
      {  // loads all interaction varables from a file
    	  //RandNum RandGenerator(seed); // call the constructor for my RandomNumber class RandNum
    	  //getline(file,material);
      	ifstream file;
      	file.open (file_name);
      	double tempE=0,tempS=0;
      	std::vector<double> E; //Energies
      	std::vector<double> SIGMA_list; // total cross section for each energy in square Angstroms

      	//n_density=0.0590073; // in per Angstroms^3	This gets changed on the next line, I put this hear just so the debuger knows it is defined.
      	file>>material>>particle_type>>interaction_type>>n_density; // reads info from file
      	while (! file.eof() ){
      		file>>tempE>>tempS;
      		E.push_back(tempE);
      		SIGMA_list.push_back(tempS);
      	}
      	// get from file when there are more specific types
      	interaction_type="inelstc";   // I need to update this so this and material etc. are all read in from the file.
      	material="Au";
      	particle_type="el";
      	file.close();
      	size=E.size();
      	current=(size/2);
      	SIGMA=SIGMA_list[current];
      	if (SIGMA==0){SIGMA=1;} // just so it doesn't divide by zero.
      	if (n_density==0){n_density=1;} // just so it doesn't divide by zero
      	LAMBDA = 1/(n_density * SIGMA);  // not sure if this is correct formula but the units work and makes sense.

        s1d=RandGenerator.spline(E, SIGMA_list);
        sMax1d=RandGenerator.spline(E, SIGMA_list);


      };
      interaction::interaction(char *file_name,int numAngles, unsigned long seed)
      	  :current(0), size(0), LAMBDA(1.), SIGMA(1.),Energy(55.),
           RandGenerator(seed),sMax1d(),s1d(),s2d(),particle_type("\n"),material("\n"),interaction_type("\n"), n_density(1.)
      {  // loads all interaction variables from a file, for Files that have angular info also
    	std::vector<double> numbers; // a temp storage for getting the numbers from the buffer
    	std::vector<double> E; //Energies
    	std::vector<double> theta; //Angles for 2d spine fit of SIGMA_angle  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! degrees or radians ???
       	std::vector<double> SIGMA_list; // total cross section for each energy in square Angstroms
		std::vector<double> MaxSigAngle; // the max value of d(sigma)/d(theta) distribution, used for random number generator
		std::vector<double> SIGMA_angle; // 2D array to account for Dsigma/Dtheta distribution

    	ifstream file;
      	file.open (file_name);

      	// get from file when there are more specific types
      	interaction_type="elstc";   // I need to update this so this and material etc. are all read in from the file.
      	material="Au";
      	particle_type="el";
      	n_density=0.0590073; // in per Angstroms^3   !!!!!!!!!!!!!!!!!!! this should come from the file as well.  Just a temp fix for now.
      	double temp=0.0;
      	char filebuffer[2048]="..";
      	while (strcmp(filebuffer,"________________________________________")){
			file.getline(filebuffer,2048);
      	}
      	file.getline(filebuffer,2048);
      	linepars(filebuffer,numbers,',');
      	numbers.erase(numbers.begin(),numbers.begin()+2); // removes the first 2 values in the array
      	theta=numbers;
      	numbers.clear();
      	//int i=0;
      	while (!file.eof()){  // gets all the data data from the file.  THis needs some error catching so it doesn't crash if something is a little different.
      		file.getline(filebuffer,2048);
      		//cout<<filebuffer<<endl;
      		linepars(filebuffer,numbers);
      		E.push_back(numbers[0]);
      		SIGMA_list.push_back(numbers[1]);
      		numbers.erase(numbers.begin(),numbers.begin()+2); // removes the first 2 values in the array
      		temp=*max_element(numbers.begin(), numbers.end());  // max_element returns a pointer to the first maximum value in the array  the * dereferences this pointer.
      		MaxSigAngle.push_back(temp);
      		//cout<<"size of SIGMA_angle= "<<SIGMA_angle.size()<<endl;
      		SIGMA_angle.insert(SIGMA_angle.end(), numbers.begin(), numbers.end());
      		numbers.clear();
      	}// while not end of file
      	//cout<<"made it out"<<endl;
      	file.close();
      	size=E.size();
      	current=(size/2);
      	SIGMA=SIGMA_list[current];
      	if (SIGMA==0){SIGMA=1;} // just so it doesn't divide by zero.
      	if (n_density==0){n_density=1;} // just so it doesn't divide by zero.
      	LAMBDA = 1/(n_density * SIGMA);  // not sure if this is correct formula but the units work and makes sense.
      	//cout<<"E= "<<E[0]<<" theta="<<theta[1]<<" SIGMA_list= "<<SIGMA_angle[0+181]<<endl;
      	s1d=RandGenerator.spline(E, SIGMA_list);
      	sMax1d=RandGenerator.spline(E, MaxSigAngle);

      	//cout<<"before spline2d"<<endl;
      	s2d=RandGenerator.spline2d(theta,E,SIGMA_angle); //
      	//cout<<"after spline2d"<<endl;

      };
      interaction::interaction()  // Sets values manually, that is right here! used primarily for testing other parts of the code
	  :current(0), size(0), LAMBDA(1.), SIGMA(1.), Energy(55.),
           RandGenerator(33u), sMax1d(), s1d(),s2d(), particle_type("el"), material("Au"), interaction_type("elstc"), n_density(0.0590073) // in per Angstroms^3
      {
    	  std::vector<double> E; //Energies
    	  std::vector<double> SIGMA_list; // total cross section for each energy in square Angstroms

          for (int i=0;i<30;i++){
          	E.push_back(i*400); // this gives me a good range from 0keV to 12,000 keV
          	SIGMA_list.push_back(6.6); // 6.6 is about the classical cross section of a Au atom in Angstroms
          };
          size=E.size();
          current=(size/2);
          SIGMA=SIGMA_list[current];
          if (SIGMA==0){SIGMA=1;} // just so it doesn't divide by zero.
          if (n_density==0){n_density=1;} // just so it doesn't divide by zero.
          LAMBDA = 1/(n_density * SIGMA);  // not sure if this is correct formula but the units work and makes sense.
          s1d=RandGenerator.spline(E, SIGMA_list);
      };  

      double interaction::str2double(char* string){
          //***********************************************
    	  // this function returns the double numeric value of a string
          // string="-345.25" then it returns -345.25
          //***********************************************
          int j=0,i=0,k;
          double var=0.0,frac=0.0;
          if (string[0]=='\0'){return 0;}
          if (string[0]=='-'){j=1;}
          for (i=j;(string[i]!='\0'&&string[i]!='.');i++){
             var=var*10. + int(string[i])-48;
          }
          if (string[i]!='\0'){i++;}
          for (k=0;string[k+i]!='\0';k++){
             frac=frac*10 + int(string[k+i])-48;
          }
          for (i=0;i<k;i++){
              frac=frac/10;
          }

          var=var+frac;
          if (j==1){var=-var;}
          return var;
      }

      void interaction::linepars(string line ,std::vector<double>& numbers,char delim){ // function to turn a line of data from a file to a array of doubles char is deliminator default is ,
    	 int i=0,q=0;
    	 char OneNum[20]="\0";
    	 while (line[i]!='\0'){ // keep going until the end of the line
    	        			q=0;
    	        			if (line[i]==delim){i++;}
    	        			while (line[i]!=delim && line[i]!='\0'){  // gets each number in the line between ','
    	        				OneNum[q]=line[i];
    	        				i++;q++;
    	        			}

    	        			OneNum[q]='\0';
    	           			numbers.push_back(str2double(OneNum));

    	        			if (line[i]!='\0'){i++;}else{line[i+1]='\0';}
    	  } // while not the end of the line

      }
      double interaction::sigma(double energy,double angle){ // Gives the value of sigma for a given energy and angle. assumes angle in degrees
    	  // This will need to find the correct index for the energy and the angle
    	  // it will then need to interperlate the value of the cross section between values of energy and angle.
    	  // The energy index can be found with the search function.
    	  // the angle index can assume that there are 0-180 degrees.
    	  return abs(spline2dcalc(s2d,energy,angle));
      }

	  double interaction::MaxSig(double energy){
			// gives the max sigma for the distribution sigma(theta) for generation random numbers
			// adds 5% just to make sure the tops don't get cut off.
    	   return abs(spline1dcalc(sMax1d,energy)*1.05); // add 5% just to make sure its above the max
		}
	  double interaction::sigma(double energy){  // returns the total cross section for given energy, also interpolates between values of energy. Changes internal values of SIGMA and LAMBDA
			SIGMA=abs(spline1dcalc(s1d,energy));
			if (SIGMA==0){SIGMA=1;} // just so it doesn't divide by zero.
			if (n_density==0){n_density=1;} // just so it doesn't divide by zero.
			LAMBDA=1/(n_density*SIGMA);
			return SIGMA;
		};
	  double interaction::sigma(){
			return SIGMA;
		};
	  double interaction::lambda(double energy){   // returns the mean free path      for given energy, also interpolates between values of energy
		  sigma(energy);
		  if (SIGMA==0){SIGMA=1;} // just so it doesn't divide by zero.
		  if (n_density==0){n_density=1;} // just so it doesn't divide by zero.
		  return 1/(n_density*SIGMA);
		};
	  double interaction::lambda(){
			return LAMBDA;
		};

	  double interaction::H(double En,double C, double x){
	  	// Gives the Energy loss to a secondary electron.  See Jablonski 2005 SIA 37 Pg 863
		if ((En*En*(C+x*x)*(C+x*x))==0){cout<<"something wrong in interaction.H"<<endl; return 0;}else{
			return 2*C*(C+En*En)*x/(En*En*(C+x*x)*(C+x*x));
		}
	  }

	  void   interaction::process(particle& p){      // for elastic interactions
			// random number for Phi
			// uses the random sampaling method to get a random theta
		  //cout<<"elastic"<<endl;
		  double Dtheta=0;
		  double DPhi=0;
		  double Enrg=p.getE();
		  double Msig=MaxSig(Enrg);

		  Dtheta=deg2rad*RandGenerator.RSMy(s2d,0,180,Msig,Enrg); // A random value chosen from the pdf
		  //Dtheta=0;  // This is wrong, delete it and it should all work fine.
		  DPhi=RandGenerator.R(0,2*Pi);
		  Vector e1=p.getDirection();
		  double x=e1(0); double y=e1(1); double z=e1(2);
		  Vector e2=e1;
		  e2(0)=-y; e2(1)=x; e2(2)=0;
		  Vector e3=e2; e3(0)=-x*z; e3(1)=-y*z; e3(2)=x*x+y*y;
		  if (y==0&&x==0){
			  e2(0)=1;e3(1)=1;
		  }else{	// checks to make sure e2 and e3 are not just the zero vector if z=1 and x and y are 0
			  // now normalize each e vector.
			  e1=e1/(sqrt(x*x+y*y+z*z));
			  e2=e2/(sqrt(x*x+y*y));
			  e3=e3/(sqrt( e3(0)*e3(0)+e3(1)*e3(1)+e3(2)*e3(2) ));
		  }
		  double p1=sqrt(Enrg);
		  //cout<<"theta="<<Dtheta/deg2rad<<" Phi="<<DPhi/deg2rad<<"  p1="<<p1<<"  Energy="<<Enrg<<endl;
		  //cout<<"e1="<<e1<<" e2="<<e2<<" e3="<<e3<<endl;

		  Vector M_old=p1*e1;  //This is really momentum.
		  Vector M2=cos(Dtheta)*e1 + sin(Dtheta)*cos(DPhi)*e2 +sin(Dtheta)*sin(DPhi)*e3;
		  M2=p1*M2;//sqrt(M2(0)*M2(0)+M2(1)*M2(1)+M2(2)*M2(2));



		  //cout<<"M_old= "<<M_old<<"   P^2="<<M_old(0)*M_old(0)+M_old(1)*M_old(1)+M_old(2)*M_old(2)<< endl;
		  //cout<<"M2= "<<M2<<" M2_mag= "<<M2(0)*M2(0)+M2(1)*M2(1)+M2(2)*M2(2)<<endl;
		  //cout<<"theta between ="<< acos((M2(0)*M_old(0)+M2(1)*M_old(1)+M2(2)*M_old(2))/(sqrt((M_old(0)*M_old(0)+M_old(1)*M_old(1)+M_old(2)*M_old(2))*(M2(0)*M2(0)+M2(1)*M2(1)+M2(2)*M2(2)))))/deg2rad<<endl<<endl;



		  //Dtheta=RandGenerator.R(0,180)*deg2rad;  //Just a random value for Dtheta
		  //cout<<"Dtheta= "<<Dtheta/deg2rad<<endl;
		  p.modify(interaction_type,0,M2); // (interaction type, phase change, momentum)  Angles are in radians.
	  };
	  void   interaction::process(std::vector<particle>& PList){ // for inelastic interactions
			// The first one modifies the old electron and the push_back creates a new electron with
			// the energy lost from the old.  $$$ What about interactions that take energy away, but
			// don't produce electrons!, plasmons for example.  This is really very incomplete in many ways.   Lots more
			// thought needs to go into all of this!!
		  	//
		  	// Things this does:  Fined the Energy lost, fine the angle the two particles go off at.

		  //cout<<"inelastic"<<endl;
		  double tempE=PList.back().getE();
		  double C=1643;  //eV as recommended in See Jablonski 2005 SIA 37 Pg 863
		  double maxY=3*sqrt(3)*(C+tempE*tempE)/(8*tempE*tempE*sqrt(C)); //0.1453*(20+tempE*tempE)/(tempE*tempE); maximum of H() pdf if C=20 real max equation is ymax=3*sqrt(3)*(C+En^2)/(8*En^2*sqrt(C))  see function H
		  double deltaE=RandGenerator.RSM(tempE,C,H,0.0,tempE,maxY); // randomly assigns a energy loss from a distribution called H() See Jablonski 2005 SIA 37 Pg 863
		  if (tempE==0){tempE=0.001;}  // just a check to make sure it doesn't divide by zero
		  double Er=abs((tempE-deltaE)/tempE);
		  double deltaT=acos(sqrt(Er));     // uses conservation of momentum with the Energy loss to calculate theta.
		  //deltaT=0; // also wrong.  Please delete !!!$$$
		  double deltaP=RandGenerator.R(0,2*Pi);


		  //$$$ This is working but I need to go through it with a fine comb to see that it is conserving momentum properly
		  Vector e1=PList.back().getDirection();
		  double x=e1(0); double y=e1(1); double z=e1(2);
		  Vector e2=e1;
		  e2(0)=-y; e2(1)=x; e2(2)=0;
		  Vector e3=e2; e3(0)=-x*z; e3(1)=-y*z; e3(2)=x*x+y*y;
		  if (y==0&&x==0){
			  e2(0)=1;e3(1)=1;
		  }else{// checks to make sure e2 and e3 are not just the zero vector if z=1 and x and y are 0
			  // now normalize each e vector.
			  e1=e1/(sqrt(x*x+y*y+z*z));
			  e2=e2/(sqrt(x*x+y*y));
			  e3=e3/(sqrt( e3(0)*e3(0)+e3(1)*e3(1)+e3(2)*e3(2) ));
		  }


		  double p1=sqrt(tempE),p2=sqrt(abs(tempE-deltaE)), p3=sqrt(deltaE);

		  //cout<<"theta="<<deltaT/deg2rad<<" Phi="<<deltaP/deg2rad<<"  p1="<<p1<<"  Energy="<<tempE<<endl;
		  //cout<<"e1="<<e1<<" e2="<<e2<<" e3="<<e3<<endl;

		  	Vector M_old=p1*e1;
		  	Vector M2=M_old; M2=p2*(cos(deltaT)*e1 + sin(deltaT)*cos(deltaP)*e2 +sin(deltaT)*sin(deltaP)*e3);
		  	Vector M3=M_old-M2;
		  	//double junk;
			  //cout<<"M_old= "<<M_old<<"  P^2="<<M_old(0)*M_old(0)+M_old(1)*M_old(1)+M_old(2)*M_old(2)<<" p1^2="<<p1*p1<<endl;
			 // cout<<"M2= "<<M2<<" M2_mag^2= "<<M2(0)*M2(0)+M2(1)*M2(1)+M2(2)*M2(2)<<" p2^2="<<p2*p2<<endl;
			 // cout<<"M3= "<<M3<<" M3_mag^2= "<<M3(0)*M3(0)+M3(1)*M3(1)+M3(2)*M3(2)<<" p3^2="<<p3*p3<<endl;
			 // cout<<"theta between M_old and M2 ="<< acos((M2(0)*M_old(0)+M2(1)*M_old(1)+M2(2)*M_old(2))/(sqrt((M_old(0)*M_old(0)+M_old(1)*M_old(1)+M_old(2)*M_old(2))*(M2(0)*M2(0)+M2(1)*M2(1)+M2(2)*M2(2)))))/deg2rad<<endl;
			 // cout<<"theta between M_3 and M2 ="<< acos((M2(0)*M3(0)+M2(1)*M3(1)+M2(2)*M3(2))/(sqrt((M3(0)*M3(0)+M3(1)*M3(1)+M3(2)*M3(2))*(M2(0)*M2(0)+M2(1)*M2(1)+M2(2)*M2(2)))))/deg2rad<<endl<<endl;

			 // cin>>junk;

			PList.back().modify(deltaE,M2,interaction_type,0); // Modify the old particle


			//Vector M=(sqrt(tempE)*M_old-sqrt(abs(tempE-deltaE))*M_new)/sqrt(deltaE); // This is the momentum
			//M3=M3/sqrt((M3(0)*M3(0)+M3(1)*M3(1)+M3(2)*M3(2)));
			PList.push_back(particle(PList.back(),deltaE,M3,interaction_type,particle_type)); // Creates a new particle
		};
