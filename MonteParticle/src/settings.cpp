#include "Settings.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

/*--------------------------------------------------------------------------------
Class constructor:
pulls in data from the text setting file passed to it one line at a time. If the
line begins with the '%' skip it, otherwise start adding each character of the line
to a buffer until '%' or the end of the line is reached. The buffer is pushed into
a vector and then the vector is parsed using a switch statement to populate the
member data
--------------------------------------------------------------------------------*/
/*settings::settings(string settingFileName){
	char fileBuffer[2048] = " ";
	char tempBuffer[200] = " ";
	ifstream inputFile;
	inputFile.open(settingFileName.c_str());
	
	while (inputFile.eof()){
		
 
	}
	inputFile.close();
}*/

settings::settings(string settingFileName){
	//Data types need to be checked
	seed = 0.0;
	elastic = " ";
	seedElastic = 0.0;
	inElastic = " ";
	seedInElastic = 0.0;

	numIncident = 0.0;
	incidentEnergy = 0.0;
	incidentEnergyDev = 0.0;

	energyThreshold = 0.0;
	boundaryEnergy = 0.0;

	//Output setting variables
	outputFileName = " ";

	spectra = false;
	spectraDouble = false;
	angle = false;
	chargeDist = false;
}

//getters and setters
double settings::getSeed(){
	return seed;
}
void settings::setSeed(double inputSeed){
	seed = inputSeed;
}

string settings::getElastic(){
	return elastic;
}
void settings::setElastic(string inputelastic){
	elastic = inputelastic;
}

double settings::getSeedElastic(){
	return seedElastic;
}
void settings::setSeedElastic(double inputSeedElastic){
	seedElastic = inputSeedElastic;
}

string settings::getInElastic(){
	return inElastic;
}
void settings::setInElastic(string inputInElastic){
	inElastic = inputInElastic;
}

double settings::getSeedInElastic(){
	return seedInElastic;
}
void settings::setSeedInElastic(double inputSeedInElastic){
	seedInElastic = inputSeedInElastic;
}

double settings::getNumIncident(){
	return numIncident;
}
void settings::setNumIncident(double inputNumIncident){
	numIncident = inputNumIncident;
}

double settings::getIncidentEnergy(){
	return incidentEnergy;
}
void settings::setIncidentEnergy(double inputIncidentEnergy){
	incidentEnergy = inputIncidentEnergy;
}

double settings::getIncidentEnergyDev(){
	return incidentEnergyDev;
}
void settings::setIncidentEnergyDev(double inputIncidentEnergyDev){
	incidentEnergyDev = inputIncidentEnergyDev;
}

double settings::getEnergyThreshold(){
	return energyThreshold;
}
void settings::setEnergyThreshold(double inputEnergyThreshold){
	energyThreshold = inputEnergyThreshold;
}

double settings::getBoundaryEnergy(){
	return boundaryEnergy;
}
void settings::setBoundaryEnergy(double inputBoundaryEnergy){
	boundaryEnergy = inputBoundaryEnergy;
}

string settings::getOutputFilename(){
	return outputFileName;
}
void settings::setOutputFileName(string inputOutputFileName){
	outputFileName = inputOutputFileName;
}

bool settings::setSpectra(){
	return spectra;
}
void settings::getSpectra(bool inputSpectra){
	spectra = inputSpectra;
}

bool settings::getSpectraDouble(){
	return spectraDouble;
}
void settings::setSpectraDouble(bool inputSpectraDouble){
	spectraDouble = inputSpectraDouble;
}

bool settings::getAngle(){
	return angle;
}
void settings::setAngle(bool inputAngle){
	angle = inputAngle;
}

bool settings::getChargeDist(){
	return chargeDist;
}
void settings::setChargeDist(bool inputChargeDist){
	chargeDist = inputChargeDist;
}