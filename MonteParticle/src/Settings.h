#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>
using namespace std;

class settings
{
	public:
		//Constructor
		settings(string settingFileName);

		//getters setters
		double getSeed();
		void setSeed(double inputSeed);

		string getElastic();
		void setElastic(string inputelastic);

		double getSeedElastic();
		void setSeedElastic(double inputSeedElastic);

		string getInElastic();
		void setInElastic(string inputInElastic);

		double getSeedInElastic();
		void setSeedInElastic(double SeedInElastic);

		double getNumIncident();
		void setNumIncident(double inputNumIncident);

		double getIncidentEnergy();
		void setIncidentEnergy(double inputIncidentEnergy);

		double getIncidentEnergyDev();
		void setIncidentEnergyDev(double inputIncidentEnergyDev);

		double getEnergyThreshold();
		void setEnergyThreshold(double inputEnergyThreshold);

		double getBoundaryEnergy();
		void setBoundaryEnergy(double inputBoundaryEnergy);

		string getOutputFilename();
		void setOutputFileName(string inputOutputFilename);
		
		bool setSpectra();
		void getSpectra(bool inputSpectra);

		bool getSpectraDouble();
		void setSpectraDouble(bool inputSpectraDouble);

		bool getAngle();
		void setAngle(bool inputSetAngle);

		bool getChargeDist();
		void setChargeDist(bool inputChargeDist);


	private:
		//Data types need to be checked
		double seed; //seed used for main
		string elastic; //First doubleeraction filename
		double seedElastic;
		string inElastic; //First inelastic doubleeraction filename
		double seedInElastic;
		
		double numIncident; //total number of incident electrons
		double incidentEnergy; //Incident electron energy in eV
		double incidentEnergyDev; // Standard deviation of incident electron energy

		double energyThreshold; //Energy where the program no longer tracks electrons again in eV
		double boundaryEnergy; //Bondary potential energy V

		//Output setting variables
		string outputFileName;

		bool spectra;
		bool spectraDouble;
		bool angle;
		bool chargeDist;




};

#endif