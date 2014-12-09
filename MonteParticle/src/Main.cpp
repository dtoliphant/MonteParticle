#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <math.h>
#include <boost/random.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "Interaction.h"
#include "Particle.h"
#include "RandomNumbers.h"

/* Electron track:  This program tracks details of incident electrons and the secondary electrons they create as they interact with a material.
 * Every individual electron-material interaction is logged and can be saved to a file for further analysis.  Details about how an electron will interact
 * is first loaded from files.
 *
 *
 * More work needs to be done:
 * Comments started with a $$$ are comments about things that need to be worked on
 */

using namespace std;
      namespace ublas=boost::numeric::ublas;

inline double Reflection(double theta, double Energy, int side ){
       double V1;
       double V2;
       double n,c;
       if (side==0){
          V1=5.+4.;V2=0;   // $$$
       }else{
          V1=0;V2=5.+4.;  // $$$
       };
       //cout<<"reflected!"<<endl;
       //$$$ need to include a if statement to handle total internal reflection.  past that angle this equation will give imaginary values
       n=sqrt((Energy-V1)/(Energy-V2));
       c=sqrt(pow((cos(theta)),2)-1+1/pow(n,2));
       return abs(pow((cos(theta)-c),2)/pow((cos(theta)+c),2));
};
double str2double(char* string){
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
inline bool crossed_boundary(Vector P1, Vector P2){  // just checks to see if the z component of position has changed from + to - or - to +
	return ( (abs(P1(2))==P1(2)) xor (abs(P2(2))==P2(2)));

};
/******************************************************************************************/
/******************************************************************************************/
int main(int argc, char *argv[]){
	std::vector<interaction> AuEl;  // An array of different types of interactions.  Elastic should always be first, then the inelastic interactions.

	char filebuffer[2048]="..";
	char tempbuffer[200]="..",numbuffer[100]="...";
	ifstream inputfile;

	inputfile.open("InPutFile.txt");
	inputfile.getline(filebuffer,2048);
	while (!inputfile.eof()){
		if (filebuffer[0]=='-'){
			char s=filebuffer[1];
			switch (s){
			case 'e':
				//an elastic interaction file needs to be loaded
				 AuEl.push_back(interaction("AllData_aoUnits.txt",181,5u));  // elastic  $$$ the second parameter (numAngles) isn't used or needed anymore.  It is still here to keep it using the correct constructor that will include angle info for the cross section.  Files must have a zero energy or Program may crash.  $$$ Should put a check in interaction class constructor that checks this and fixes it if it doesn't
				break;
			case 'i':
				//an inelastic interaction file needs to be loaded
				int i=3,j=0;
				unsigned int seed=10;
				while (filebuffer[i]!=' '){tempbuffer[j]=filebuffer[i];i++;j++;}
				tempbuffer[j]='\0';
				j=0;
				while (filebuffer[i]!=' '){numbuffer[j]=filebuffer[i];i++;j++;}
				numbuffer[j]='\0';
				seed=int(str2double(numbuffer));

				AuEl.push_back(interaction(tempbuffer,seed));  // inelastic, first argument is file with sigma(energy) or sigma(energy,angle) info.  Last argument is seed for random number generator.
				break;

			}
		}

	}






  // set up interaction class list
  AuEl.push_back(interaction("AllData_aoUnits.txt",181,5u));  // elastic  $$$ the second parameter (numAngles) isn't used or needed anymore.  It is still here to keep it using the correct constructor that will include angle info for the cross section.  Files must have a zero energy or Program may crash.  $$$ Should put a check in interaction class constructor that checks this and fixes it if it doesn't
  AuEl.push_back(interaction("sigma_i_Au.txt",12u));  // inelastic, first argument is file with sigma(energy) or sigma(energy,angle) info.  Last argument is seed for random number generator.

  RandNum RandClass(33u);  //Random numbers class to use for random numbers from pdf or just rand between 0..1, or other.  argument is seed.

  // set up data files that will store the results.
  ofstream datafile,datafile2;
  datafile.open ("DataFile.txt");
  datafile2.open ("DataFile2.txt");

  Vector temp, temp2;
  std::vector<particle> el;  // defines an array of electron classes called el.  Each class represents an electron and stores all of its history.
  el.reserve(12);        // reserve memory for at least 12 electron classes in el

  double sigma_t=0,maxE;
  int I_type=0, material=0; // 0 for I_type must be elastic type.

  unsigned long TotalNumIncElectrons=100; // determines the number of incident electrons 69 electrons/s at 900eV  27 electrons/s for 2keV,
  cout<<"input the total number of electrons you want to run"<<endl;
  cin>>TotalNumIncElectrons;
  cout<<"input the initial incident electron energy"<<endl;
  cin>>maxE;

  time_t rawtime;
  time (&rawtime);  // used for timing to see how long a run took.
  double timeleft=0.,timestart=clock()/CLOCKS_PER_SEC;

// Starts the work of making the calculations!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  /******************************************************************************************/
  for (unsigned long i=0;i<TotalNumIncElectrons;i++){
	  // a progress indicator
	  if (fmod(i,1000.0)==0){
		  timeleft=(clock()/CLOCKS_PER_SEC-timestart)*(TotalNumIncElectrons/(i+.01) -1);
		  cout<<100.0*i/TotalNumIncElectrons<<"%    estimated time left= "<<timeleft<< "s "<<timeleft/60<<"min  "<<timeleft/3600<<"hrs \n"<<endl;
	  }
	  material=0; // 0 is vacuum, 1 is in the bulk, $$$ will need to create a new class for material to handle more complicated geometry
      el.push_back(particle(RandClass.Normal(maxE,10))); // initializes a particle with a random energy uses defaults in partial constructor for other parameters.
      while (!el.empty()){ //continues until the electron has come to a stop and all of its secondaries have also.
// First Diamond in Algorithm F(E,z)<L
    	  while (el.back().getE()>4.0){  // checks energy for minimum value  $$$ this value should be read in from external start-up file.
    		  I_type=0;  //0 is first interaction class an means it is an elastic interaction
    		  double tempE=el.back().getE();
    		  sigma_t=0;
// Rectangle Calculate sigmas
    		  for (unsigned int j=0; j<AuEl.size(); j++){  // calculate sigma total .. in the process the interaction class will save the individual sigmas for the last energy given
    			  sigma_t=sigma_t+AuEl[j].sigma(tempE); //el.back().getE());
    		  }
    		  while (I_type==0){ // while it is an elastic type collision, always runs through the first time
    			  double rand= RandClass.R(0,1);
    			  double sigma_part=0;
    		      //if (el.back().ID>78816000){cout<<"3 "<<endl;}
    			  while (rand>=sigma_part){  // Step by step adds up interaction sigmas determines what type of interaction will be used
    				  sigma_part=sigma_part+(AuEl[I_type].sigma()/sigma_t);  // this is making sigma_part greater than 1 1.0158164 so that there is an infinite loop!!!!!!!!!!!!!!!!
    				  I_type++;
    			  }
    			  I_type--;
    			  double s; // path length
    		      s=-AuEl[I_type].lambda()*log(RandClass.R(0,1));  // s=-l*ln(rand)  is a random distance with appropriate mean free path lambda
    			  temp=el.back().getPosition();  // saves position to use for cross boundary check
    			  el.back().move(s);  // calculate new position
    			  temp2=el.back().getPosition();
    			  bool cross;
    			  cross=crossed_boundary(temp,el.back().getPosition());
    			  if (cross==1){  //crossed boundary
    				  if (RandClass.R(0,1)<= Reflection(el.back().gettheta(),el.back().getE(), material) ){ // assumes boundary is only a plane normal to z
    					  //  reflected
    					  // changes theta and position to account for reflection
        				el.back().Reflected();
    				  }else{  // didn't reflect
        				material=!material;  // changes the material number
        				// change theta, s and energy,  will not change s because its just in vacuum for now will need to change later for more complicated geometry
        				el.back().BoundaryCross(4.0*material,4.0*!material);  // changes E and theta based on boundary and direction, $$$ NEED TO CHANGE THE 4 TO THE CORECT VALUE FOR GOLD
    				  }
    			  } // if crossed boundary
    			  if (material == 0){I_type =-1;} // so it will get out of the while loop for elastic collisions
    			  if (material != 0 && I_type ==0){  // then not in vacuum and an elastic collision
    				  // process elastic collision
    				  AuEl.at(I_type).process(el.back());
    			  }
        		}  //while (I_type==0) // while it is an elastic type collision
    	    	if (material == 0){  // then in vacuum so stop.
    	    		I_type=-1;  // there can be no interaction when its in the vacuum
    	    		break;  // exit and save
    	    	}
    	    // Process inelastic collision
    	    AuEl[I_type].process(el);
      	} //while (el.back().getE()>4.0)  check for exit condition
      	el.back().save2file(&datafile,&datafile2,0); // Save electron then delete it from list
      	el.pop_back();
      	material=1; // set it back in the bulk so that next electron has correct value.
      }  // while to finish up the secondaries;
    } // for  done with all secondaries etc. so saved and create a new primary


  	el.push_back(900); // just junk so I can use save2file function one last time to save a histogram.
  	el.back().move(1);  // have to do this so it trigers the right if statement in save2file ... $$$ fix!
  	el.back().save2file(&datafile,&datafile2, 1);  // the 1 indicates this is the last time this will be called.  So save any data to the file.


    el.clear();  // clears all the memory tied up for el
    AuEl.clear(); // clears all memory tied up for the interactions
    std::cout<<"time running: "<<clock()/CLOCKS_PER_SEC-timestart<<"s"<<std::endl;
    datafile.close();
    datafile2.close();

    return EXIT_SUCCESS;
} // main
