/*
 * CAverage.cpp
 *
 *  Created on: May 4, 2010
 *      Author: radtke
 */
/*
 * Version modified by Johann Hansing. Original Version by Matthias Radtke
 */

#include "headers/CAverage.h"


using namespace std;

CAverage::CAverage() {
}

CAverage::CAverage(string name, string folder, int datapoints)
    : _instantValues (datapoints), _errors (datapoints)      //vectors will automatically be initialized to zero
{
	_points = datapoints;
	_name = name;

	clear();

	// create empty save file..
	ofstream myfile;
	string s;
	ostringstream outStream;
	outStream << folder << "/InstantValues/" << _name << ".txt";
	s = outStream.str();
	_file_instV = s;
	myfile.open(_file_instV.c_str());
	myfile.close();
}


/* adds Value v as instant Value to vector at datapoint n*/
void CAverage::addInstantValue(vector<double> v) {
    for (int i=0; i<_points; i++){
        _instantValues[i] += v[i];
        _errors[i] += v[i]*v[i];
     }
    //Formula from Wikipedia: Standartabweichung: Bearbeitung fuer auflaufende Messwerte + Standardfehler
    //continues saveAverageInstantValues()!
}


void CAverage::clear() {
	_instantValues.clear();
//	_averageValues.clear();
}

/* save instant values to file */

void CAverage::saveAverageInstantValues(int counter) {
    //timeInt is the interval between the times that an average instant value was saved
	if (_points > 0) {
		ofstream myfile;
		myfile.open(_file_instV.c_str());


		for (unsigned int i = 0; i < _points; i++) {
			myfile << (_instantValues[i]/counter) << " " <<
			//the next formula is to calculate the error from the standard deviation
			        //CAREFUL No error for counter = 1, due to divide by zero!
			sqrt((_errors[i] - _instantValues[i] * _instantValues[i] / counter)/(counter*(counter - 1))) <<
			endl;
		}
		myfile.close();

		_instantValues.clear();
	}
}

