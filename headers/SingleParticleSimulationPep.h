
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <boost/filesystem.hpp>
#include "CAverage.h"
#include "CConfiguration.h"


using namespace std;


//Function declarations
template<typename T>
string toString(const T& value){
    ostringstream oss;
    oss << fixed << value;  // << fixed is needed in order to avoid scientific notation
    return oss.str();
}

string numString( float x ){
    // from http://stackoverflow.com/questions/18881854/how-to-output-float-to-cout-without-scientific-notation-or-trailing-zeros
    ostringstream strout ;
    strout << fixed << x ;
    string str = strout.str() ;
    size_t end = str.find_last_not_of( '0' ) + 1 ;
    str=str.erase( end );
    if (str.back()=='.') return str.substr(0, str.length() - 1);
    return str;
}


template <typename T, size_t N>
inline
size_t sizeOfArray( const T(&)[ N ] )
{
  return N;
}


string createDataFolder(string trigger, double timestep, double simtime, double potRange, double potStrength,
                        double particlesize, bool steric, bool ranU, double dvar, double polydiam, string peptide, double uDebye, double uBend){
    //NOTE: Maybe I can leave out dt, as soon as I settled on a timestep
    //NOTE: As soon as I create input-list with variables, I must change this function
    char range[10];
    sprintf(range, "%.3f", potRange);
    //In the definition of folder, the addition has to START WITH A STRING! for the compiler to know what to do (left to right).
    string folder = "sim_data";
    if (trigger!="x") folder += "/" + trigger; 
    folder += "/" + peptide;
    if (ranU) folder += "/ranU";
    if (steric) folder = folder + "/steric"; 
    folder += "/d" + numString(dvar);
    folder = folder
            + "/dt" + numString(timestep)
            + "/t" + numString(simtime)
            + "/a" + numString(polydiam)
            + "/p" + numString(particlesize)
            + "/uD" + numString(uDebye)
            + "/uB" + numString(uBend)
            + "/k" + range
            + "/u" + numString(potStrength);
    boost::filesystem::create_directories(folder);
    boost::filesystem::create_directory(folder + "/InstantValues");
    boost::filesystem::create_directory(folder + "/Coordinates");
    return folder;
}


void settingsFile(string folder, string trigger, double particlesize, double timestep, double runs, double steps, double potStrength, double potRange,
        bool recordMFP, bool steric, bool ranU, double dvar, double polydiam, string peptide, double uDebye, int badsteps, double uBend){
    //Creates a file where the simulation settings are stored
    //MAYBE ALSO INCLUDE TIME AND DATE!!
    ofstream settingsfile;
    settingsfile.open((folder + "/sim_Settings.txt").c_str());
    settingsfile << "Sim dir: " << folder << endl;    
    settingsfile << "trigger: " << trigger << endl;   
    settingsfile << "Peptide " << peptide << endl;
    settingsfile << "TMP " << recordMFP << endl;//" (Bessel)" << endl;  //TODO Bessel!
    settingsfile << "includesteric " << steric << endl;
    settingsfile << "ranU " << ranU  << endl;
    settingsfile << "p " << particlesize << endl;
    settingsfile << "dt " << timestep  << endl << "runs " << runs << endl << "steps " << steps << endl << "time: " << timestep*steps << endl;
    settingsfile << "k " << potRange << endl << "U_0 " << potStrength << endl;
    settingsfile << "dvar " << dvar << endl;
    settingsfile << "a " << polydiam << endl;
    settingsfile << "uDebye " << uDebye << endl; 
    settingsfile << "uBend " << uBend << endl; 
    settingsfile << "\nBADSTEPS " << badsteps << endl;


    settingsfile.close();
}
