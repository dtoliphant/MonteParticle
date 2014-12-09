#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <boost/random.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "Interaction.h"
#include "Particle.h"

namespace ublas=boost::numeric::ublas;

// a few global variables :
unsigned long long GlobalCounter=0; // used for new particle ID.
unsigned int tempCounter=0; // $$$ probably a better way to do this, but used for a switch to tell if its the last time save2file is called so it acualy saves
double SpetraHistogram[301][301] ={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//long SpetraHistogram2[301][301] ={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
      
//constructor used mostly for incident electrons
particle::particle(const double& En,const Vector& R,const double& T, const double& P,string Int_type,string p_type)
		: tempPosition(R), Phat(boost::numeric::ublas::zero_vector<double> (3)) {

	particle_type=p_type;
	interaction.push_back(Int_type);
    E.push_back(En);
    //tempPosition=R;
    Position.push_back(tempPosition);
    tempPosition(2)=tempPosition(2)-0.001;
    theta=T; phi=P;
    Phat(0)=sin(theta)*cos(phi); Phat(1)=sin(theta)*sin(phi); Phat(2)=cos(theta);
    phase=0.;
    GlobalCounter++;
    //cout<<GlobalCounter<<" primary"<<endl;
    OID=0; OE=0.; ID=GlobalCounter;
              
};
//constructor used mostly for secondary electrons
particle::particle(const particle oldparticle,const double& En,const Vector& M,string Int_type, string p_type)
		: tempPosition(oldparticle.tempPosition), Phat(boost::numeric::ublas::zero_vector<double> (3)) {

	particle_type=p_type; // defalt "el" for electron
    interaction.push_back(Int_type); // defalt "Pr" for primary, other posibilities would desribe details of the particals interaction
    E.push_back(En); // the energy
    //tempPosition=oldparticle.tempPosition;
    Position.push_back(tempPosition);
    
    //M=M/sqrt((M(0)*M(0)+M(1)*M(1)+M(2)*M(2)));


    Phat=M/sqrt((M(0)*M(0)+M(1)*M(1)+M(2)*M(2)));

    phase=0.0;
    theta=acos(Phat[2]);
    //cout<<"Constructor theta="<<theta<<" from Momenum[2]="<<Momentum[2]<<endl;
    phi=acos(Phat[0]/sin(theta));
    //Momentum[0]=sin(theta)*cos(phi); Momentum[1]=sin(theta)*sin(phi); Momentum[2]=cos(theta);

    GlobalCounter++;
    //cout<<GlobalCounter <<" secondary"<<endl;
    OID=oldparticle.ID; OE=oldparticle.E.back()+En; ID=GlobalCounter;
           
};

double particle::getE(void){
	return E.back();
};
double particle::gettheta(void){
	return theta;
};
Vector particle::getPosition(void){
	return tempPosition;
};
Vector particle::getDirection(void){
	return Phat;
};

void particle::move(double r){  // moves the particle a distance r in the direction of motion
	tempPosition=tempPosition+(r*Phat);
};
void particle::Delta_direction(const double& D_theta, const double& D_phi){
	theta=theta+D_theta;
	//cout<<"theta="<<theta;
	if (theta>pi){theta=theta-pi;}
    phi=phi+D_phi;
    //cout<<"  phi="<<phi<<endl;
    if (phi>2*pi){phi=phi-2*pi;}
    Phat[0]=sin(theta)*cos(phi); Phat[1]=sin(theta)*sin(phi); Phat[2]=cos(theta);
};

void particle::New_direction(const double& N_theta, const double& N_phi){
	theta=N_theta;
    phi=N_phi;
    Phat[0]=sin(theta)*cos(phi); Phat[1]=sin(theta)*sin(phi); Phat[2]=cos(theta);
};

void particle::modify(const double& Delta_E,const Vector& M,string Int_type,double Delta_phase){
	E.push_back(E.back()-Delta_E);
    interaction.push_back(Int_type);
    Position.push_back(tempPosition);
    //cout<<Int_type<<endl;
    //M=M/sqrt((M(0)*M(0)+M(1)*M(1)+M(2)*M(2)));
    Phat=M/sqrt((M(0)*M(0)+M(1)*M(1)+M(2)*M(2)));
    theta=acos(Phat[2]);
    if (theta!=0){phi=acos(Phat[0]/sin(theta));}else{phi=0;}

    phase=phase+Delta_phase;
};
void particle::modify(const double& Delta_E,const double& Delta_theta,const double& Delta_phi,string Int_type,double Delta_phase){
	E.push_back(E.back()-Delta_E);
    interaction.push_back(Int_type);
    Position.push_back(tempPosition);
    //cout<<Int_type<<endl;
    Delta_direction(Delta_theta,Delta_phi);
    phase=phase+Delta_phase;
};
void particle::modify(const double& Delta_theta,const double& Delta_phi,string Int_type,double Delta_phase){
	E.push_back(E.back());
    interaction.push_back(Int_type);
    Position.push_back(tempPosition);
    //cout<<Int_type<<endl;
    Delta_direction(Delta_theta,Delta_phi);
    phase=phase+Delta_phase;
    
};
void particle::modify(string Int_type,double Delta_phase,const Vector& M){
	E.push_back(E.back());
    interaction.push_back(Int_type);
    Position.push_back(tempPosition);

    Phat=M/sqrt((M(0)*M(0)+M(1)*M(1)+M(2)*M(2)));
    theta=acos(Phat[2]);
    if (theta!=0){phi=acos(Phat[0]/sin(theta));}else{phi=0;}

    phase=phase+Delta_phase;

};
void particle::modify(){
	cout<<E.size()<<" Energy values"<<endl;
	cout<<Position.size()<<" Position values"<<endl;
	cout<<interaction.size()<<" interaction values"<<endl;
	
}

void particle::save2file(ofstream *file_ptr,ofstream *file_ptr2, int last){ // defalt last=0 ie not the last one, so don't save yet.
	// This function is only called when the electron is done, ie got to an energy less than minimum or exited the material.
	// so tempPosition is really last position that is measured.

// outputs the last position or first position of primaries or secondaries.
/*	double r=sqrt(tempPosition(0)*tempPosition(0)+tempPosition(1)*tempPosition(1));
	if (tempPosition(2)>0 && OID==0){

		*file_ptr<<r<<" "<<tempPosition(2)<<"\n";
		//*file_ptr<<tempPosition(0)<<" "<<tempPosition(1)<<" "<<tempPosition(2)<<"\n";
		//Vector start=*Position.begin();
		//cout<<Position.size()<<endl;
		//*file_ptr2<<start(0)<<" "<<start(1)<<" "<<start(2)<<"\n";

	}else{
		if(tempPosition(2)>0){
			*file_ptr2<<r<<" "<<tempPosition(2)<<"\n";
			 //*file_ptr2<<tempPosition(0)<<" "<<tempPosition(1)<<" "<<tempPosition(2)<<"\n";
		}
	}*/

	Vector start=*Position.begin();
    if (tempPosition(2)>0||OID!=0){  // just the ones that make it out or just the ones that stay in.
    	 //*file_ptr<<"Ejected: ";
    	 //cout<<"E:";
    	//if (ID>78826178){cout<<"C "<<endl;}
    	 if (tempCounter==0){
    	    for (int i=0;i<301;i++){
    	    	for (int j=0;j<301;j++){
    	    	SpetraHistogram[i][j]=0;
    	    	}
    	    }

    	    tempCounter=1;
    	 }
    	 //if (ID>78826178){cout<<"D "<<endl;}
    	 if (last==1){
    		 cout<<"last ID="<<ID<<endl;
    		 //*file_ptr<<"## ##\n"<<"@incident electron number = \n"<<"E Counts\n";
    		 for (int i=0;i<300;i++){
    			for (int j=0;j<300;j++) {
    				//if (SpetraHistogram[i][j]>0){SpetraHistogram[i][j]=log(SpetraHistogram[i][j]);}  // this is so that I can graph it and still see the much smaller negative values
    			 *file_ptr<<i<<" "<<j<<" "<<SpetraHistogram[i][j]<<"\n";
    			 //*file_ptr2<<i<<" "<<j<<" "<<SpetraHistogram2[i][j]<<"\n";
    			}
    		 }
    		 *file_ptr<<endl;
    		 *file_ptr2<<endl;
    		
    		 cout<<"Saved"<<endl;
    	 }else{
    		// if (ID>78826178){cout<<"E "<<endl;}
    		 double r=sqrt(tempPosition(0)*tempPosition(0)+tempPosition(1)*tempPosition(1));
    		 double rst=sqrt(start(0)*start(0)+start(1)*start(1));
    		 double z=tempPosition(2);
    		 double zst=start(2);
    		 //if (ID>78826178){cout<<"F "<<" r="<<r<<"  rst="<<rst<<"  z="<<z<<"  zst="<<zst<<endl;}
    		 if (r<300.&&z<300.){
    			 //SpetraHistogram[(int)E.back()]++;
    			 if (z>0){  // to make sure it is only the electrons that are in the material that are counted.
    				 SpetraHistogram[(int)r][(int)z]--;
    				 //SpetraHistogram2[(int)z]--;
    				 //if (ID>78826178){cout<<"G "<<endl;}
    			 }

    			 if (OID!=0&&rst<300&&zst<300){ // so that it doesn't count the incident electrons.  Just the secondaries and those are all created in the material so no need to check z
    				 SpetraHistogram[(int)rst][(int)zst]++;

    				 //SpetraHistogram2[(int)zst]++;
    				 //if (ID>78826178){cout<<"H "<<endl;}
    			 }
    		 }
    		 //*file_ptr<<Position.back()<<"\n";
    	 }
    	 
    	 
    	 
    	 //*file_ptr<<E.back()<<",";
     }else{
    	 //*file_ptr<<"Absorbed: ";
    	 //cout<<"A:";
     }


     /*
     //*file_ptr<<OID <<","<<OE <<","<<particle_type<<","<<ID <<";";
     cout<<particle_type<<ID<<","<<OID <<",["<<OE<<"]; "<<(char)9;  // char 9 is a tab
     vector<Vector>::iterator P;
     vector<double>::iterator En;
     vector<string>::iterator S;
     Vector tempP;
     En=E.begin();
     S=interaction.begin();
     int i=0;
     for (P=Position.begin() ; P < Position.end(); P++, En++, S++ ){
    	 tempP=*P;
         //*file_ptr <<*S <<","<<*En <<"," << *P <<"; ";
         cout <<*S <<",("<< tempP(0)<<","<<tempP(1)<<","<<tempP(2)<<"),["<<*En<<"]; "<<"\n"<<(char)9;
         i++;
     };
     //*file_ptr <<"Theta: " <<theta*180/pi<<",phi:" <<phi*180/pi<<"; ";
     //cout<<"Theta: " <<theta*180/pi<<",phi:" <<phi*180/pi<<"; ";
     //*file_ptr<<"\n\n"; // change to "\n" when done with debuging
     //cout<<"tempP="<<tempPosition;//<<" Phat="<<Phat<<"" ;
     //cout<<"\n"<<endl;
     //cout<<"\n";
     */
}; // save2file function
void particle::Reflected(void){ // will need to modify this so that it can handle reflecting off of any plane, now just xy plane
	theta=pi-theta;
	Phat[0]=sin(theta)*cos(phi); Phat[1]=sin(theta)*sin(phi); Phat[2]=cos(theta);
	tempPosition(2)=-tempPosition(2);
	
	
	
};
void particle::BoundaryCross(double V1,double V2){ // will need to modify this so that it can handle a boundary of any plane, now just xy plane
	double n;

	if ( (E.back()-V2)!=0) {n=sqrt(abs((E.back()-V1)/(E.back()-V2)));}else{n=1;}
	theta=asin(n*sin(theta));
	Phat[0]=sin(theta)*cos(phi); Phat[1]=sin(theta)*sin(phi); Phat[2]=cos(theta);
	E.back()=E.back()+(V1-V2);
	
};
