#include "Settings.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

/*--------------------------------------------------------------------------------
Class constructor:
pulls in data from the text setting file passed to it one line at a time. If the
line begins with the '%' skip it, otherwise start adding each character of the line
to a buffer until '%' or the end of the line is reached. The buffer is pushed into
a vector and then the vector is parsed using a switch statement to populate the
member data
--------------------------------------------------------------------------------*/
settings::settings(string settingFileName){
	ifstream inputFile;
	inputFile.open(settingFileName);
	string line;
	vector<string> container;

	while (getline(inputFile, line)){
		if (line[0] != '%'&& line.size() > 0)
			container.push_back(line);

	}

	for (int i = 0; i < container.size(); i++){
		string temp = container[i];
		
		string buffer;
		istringstream strm;
		strm.str(temp.substr(3));
		
		switch (temp[1]){
		case 's':
			strm >> buffer;
			seed = stod(buffer);
			strm.clear();
			break;
		case 'e':
			strm >> elastic >> buffer;
			seedElastic = stod(buffer);
			strm.clear();
			break;
		case 'i':
			strm >> inElastic >> buffer;
			seedInElastic = stod(buffer);
			strm.clear();
			break;
		case 'n':
			strm >> buffer;
			numIncident = stod(buffer);
			strm.clear();
			break;
		case 'E':
			strm >> buffer;
			incidentEnergy = stod(buffer);
			strm.clear();
			break;
		case 'd':
			strm >> buffer;
			incidentEnergyDev = stod(buffer);
			strm.clear();
			break;
		case 'm':
			strm >> buffer;
			energyThreshold = stod(buffer);
			strm.clear();
			break;
		case 'V':
			strm >> buffer;
			boundaryEnergy = stod(buffer);
			strm.clear();
			break;
		default :
			break;
		}
		buffer = " ";
	}
	inputFile.close();
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