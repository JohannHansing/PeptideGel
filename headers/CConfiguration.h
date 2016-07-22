#ifndef CCONFIGURATION_H_
#define CCONFIGURATION_H_


#include <array>
#include <string>
#include <vector>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <boost/random.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <Eigen/Dense>
#include "CRod.h"
#include "CBead.h"

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// TODO CLEAN THIS UP AFTER INCLUSION OF BEADS! 
///////////////////////////////////////////////////////////////////////////////////////
/* MORE TODO FOR BEADS:
- write trajectories of all beads to file. 
- write energies to file for partition coefficient
*/


#define ifdebug(x) 

using namespace std;

class CConfiguration {
    /*Class where all the configuration variables such as potRange etc. and also most functions for the
     * simulation are stored
     */
private:
    //MISC
    FILE* m_traj_file;
    
    //SCALING
    double _timestep;         //This is the RESCALED timestep! timestep = dt * kT / (frictionCoeffcient * particlesize)
    double _mu_sto;
    double _pradius;     //particle size is most relevant for scaling! (default = 1)
    double _boxsize[3] = {10,10,10};          // ALWAYS define boxsize through particlesize due to scaling!
    double _bdef = 10;     //default boxsize
    double _epsilon;

    //EXPONENTIAL Potential
    double _potRange;         // Avoid values like 10/3 = 3.33333... ALWAYS define the Range of the exp potential through boxsize due to scaling!
    double _potStrength;      // rescaled U_0 for exponential Potential
    double _rodDistance;
    double _cutoffExpSq;
    double _cylLJSq;    //cutoff for Lennard-Jones calculation between cyl and bead (at minimum)
    
    


    //bool Parameters
    bool _potMod;             // true if the modified exponential potential version is used that is not 3*U_0 at the intersections.
    bool _LJPot;              // if true, then LJ interaction is calculated for steric hindrance
    bool _ranU;
    bool _hpi;
    bool _noCyl=false;

    //COUNTERS AND INIT VALUES
    double _boxCoord[3];
    unsigned int _wallcrossings[3]; //counts how many times the particle has crossed the [entry, opposite, side] wall, while
                                            //travelling through the lattice
    int _entryside[3];            //records through which side and in which direction the particle last entered a box. For the face of the cube in the x-z
                                            //plane at y=0, it is entryside[0] = 1, in the x-y plane at z=L it is entryside[1] = -1!
    double _resetpos;
    double _startpos[3];          //Stores where the particle starting position was. This is needed to calculate the mean square displacement
    Eigen::Vector3d _prevpos;           //Stores previous particle position before particle is moved.
    std::array<std::array<double, 16>, 3> _distarr;

    int _min, _max;        // parameters for determining up to which order neighboring rods are considered for the potential

    //Particle parameters
    Eigen::Vector3d _ppos;    //initialize particle position (DO IT LIKE resetpos FOR MOVEPARTICLEFOLLOW/RESET)
    double _f_mob[3];   //store mobility and stochastic force
    double _f_sto[3];
    // rod parameters
    double _polyrad;
    double _polydiamSq;
    
    //BEAD parameters
    vector<CBead> _beads;
    unsigned int _N_beads;
    
    //---------------------------- PEPTIDE ------------------------------
    
    void initPeptide(string peptide){
        //init Beads as array of N beads
        if (peptide=="3block" || peptide=="3rep" || peptide=="3att" || peptide=="3neut"){
            _beads.clear();
            const int len = 3;
            std::array<int,len> signs={1,0,1};
            if (peptide=="3block") signs[0]=-1;
            else if (peptide=="3att"){
                signs[0]=-1;
                signs[2]=-1;
            }
            else if (peptide=="3neut"){
                array<int,len> tmp={0,0,0};
                signs=tmp;
            }
            Eigen::Vector3d beadpos = _ppos;
            for (int i=0;i<len;i++){
                CBead bead = CBead(signs[i],beadpos);
                _beads.push_back(bead);
                beadpos(0) += _r0SP;
                ifdebug(
                    cout << "beadpos\n" << beadpos << endl;
                    cout << "_ppos\n" << _ppos << endl;
                    cout << "beadpos\n" << bead.pos << endl;)
            }
        }
        else if (peptide=="8NNblock" || peptide=="8NNalter" || peptide=="8NNrep" || peptide=="8NNatt"){
            _beads.clear();
            const int len = 8;
            std::array<int,len> signs;
            if (peptide=="8NNblock"){
                array<int,len> tmp={1,1,1,1,-1,-1,-1,-1};
                signs=tmp;
            }
            else if (peptide=="8NNalter"){
                array<int,len> tmp={1,-1,1,-1,1,-1,1,-1};
                signs=tmp;
            }
            else if (peptide=="8NNrep"){
                array<int,len> tmp={1,1,1,1,1,1,1,1};
                signs=tmp;
            }
            else if (peptide=="8NNatt"){
                array<int,len> tmp={-1,-1,-1,-1,-1,-1,-1,-1};
                signs=tmp;
            }
            Eigen::Vector3d beadpos = _ppos;
            for (int i=0;i<len;i++){
                CBead bead = CBead(signs[i],beadpos);
                _beads.push_back(bead);
                beadpos(0) += _r0SP;
                ifdebug(
                    cout << "beadpos\n" << beadpos << endl;
                    cout << "_ppos\n" << _ppos << endl;
                    cout << "beadpos\n" << bead.pos << endl;)
            }
        }
        else if (peptide=="11block" || peptide=="11alter" || peptide=="11rep" || peptide=="11att"){
            _beads.clear();
            const int len = 11;
            std::array<int,len> signs;
            if (peptide=="8NNblock"){
                array<int,len> tmp={1,0,1,0,1,0,-1,0,-1,0,-1};
                signs=tmp;
            }
            else if (peptide=="8NNalter"){
                array<int,len> tmp={1,0,-1,0,1,0,-1,0,1,0,-1};
                signs=tmp;
            }
            else if (peptide=="8NNrep"){
                array<int,len> tmp={1,0,1,0,1,0,1,0,1,0,1};
                signs=tmp;
            }
            else if (peptide=="8NNatt"){
                array<int,len> tmp={-1,0,-1,0,-1,0,-1,0,-1,0,-1};
                signs=tmp;
            }
            Eigen::Vector3d beadpos = _ppos;
            for (int i=0;i<len;i++){
                CBead bead = CBead(signs[i],beadpos);
                _beads.push_back(bead);
                beadpos(0) += _r0SP;
                ifdebug(
                    cout << "beadpos\n" << beadpos << endl;
                    cout << "_ppos\n" << _ppos << endl;
                    cout << "beadpos\n" << bead.pos << endl;)
            }
        }
        else{
            cout << "\nBad peptide name!\nAborting..." << endl;
            abort();
        }
        
         _N_beads = _beads.size();
    }
    
    //---------------------------- POTENTIALS ------------------------------
    double _uLJ;
    
    //SPRING potential
    double _r0SP; //equilirbium distance of beads for peptide is defined by LJ pot minimum
    double _kappaSP = 100;
    double _uSpring;
    double _uCylTot;
    
    //DEBYE potential
    double _uDebye;
    
    //BENDING potential
    double _uBend = 3; // This determines how stiff the peptide is. Judging from looking at VMD simulations, _uBend=3 seems appropriate
    
    void calcLJPot(const double rSq, double& U, double& Fr, double stericSq){
        //Function to calculate the Lennard-Jones Potential
        double  por6 = pow((stericSq / rSq),3); //por6 stands for "p over r to the power of 6" . The 2 comes from the fact, that I need the particle radius, not the particle size
        ifdebug(if (4. * ( por6*por6 - por6 + 0.25 ) > 50) cout << "Very large LJ!!"<< endl;)
        U += 4. * ( por6*por6 - por6 + 0.25 );
        Fr +=  24. / ( rSq ) * ( 2. * por6*por6 - por6 );
    }

    void calcSpringPot(const double& r, double &U, double &Fr) { 
        //double u = _kappaSP * pow(r - _r0SP, 2);
        //if (u > 10000) cout << "U = " << u << "r = " <<  r << endl;
        U = (_kappaSP * pow(r - _r0SP, 2));
        Fr += (_kappaSP * 2. * (_r0SP/r - 1.));
    }    
    
    void calcDebyePot(const double r, double& U, double& Fr, int att_rep){
        // att_rep is positive or negative prefactor that determines the sign of the interaction
        double utmp = att_rep * _uDebye * exp(-r / _potRange) / r;
        U = utmp;
        Fr += utmp * (_potRange + r)/(_potRange*r*r);
        //Fr += utmp * (1./ (_potRange * r) + 1./(r*r) );
    }
    
    void addBendingPot(const double r, double& U, double& Fr){
        // simple repulsive potential to create stiffness for peptide
        U = -_uBend * r / _pradius;
        Fr += _uBend/(r * _pradius);
    }


    boost::mt19937 *m_igen;                      //generate instance of random number generator "twister".
    double zerotoone(){
        boost::uniform_01<boost::mt19937&> dist(*m_igen);
        return dist();
    }

    double atob(double a, double b){
        ifdebug(if (a >= b) cout << "a " << a << "\nb " << b << endl;)
        boost::variate_generator<boost::mt19937&, boost::random::uniform_real_distribution<>> ran_gen(*m_igen, boost::random::uniform_real_distribution<>(a, b));
        return ran_gen();
    }

    int ran_sign(){
    // the variate generator uses _igen (int rand number generator),
    // samples from uniform integer distribution 0, 1
        boost::variate_generator<boost::mt19937&, boost::uniform_int<>> zeroone(*m_igen, boost::uniform_int<>(0, 1));
	return (zeroone() * 2) - 1; //this calculation makes value either -1 or 1
    }

    
    double getfixb(){
        // helper function to fix boxsize
        return _bdef;
    }

    
    //************** Rand ***************    
    double ran_norm(){
        boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > ran_gen(
                *m_igen, boost::normal_distribution<double>(0, _dvar));
        return ran_gen(); 
    }
    
    double _dvar = 1.;
   
    
    std::array<std::array<std::array<CRod, 4>, 4>, 3> _drods; // array to store polymer rods in simulation box, the outermost array stores polymers that are parallel to the same axis
    void initRand(){
        // Initialize system with random rod displacement d
        bool overlaps;
        double xipos, xjpos;
        double tmp[4] = {-_bdef,0.,_bdef,2*_bdef};
        for (int axis=0;axis<3;axis++){//axis 0 is x axis.
            int i,j;
            i=axis+1;
            if (i==3) i=0;
            j=3-(i+axis);
            for (int abcd=0;abcd<4;abcd++){
                for (int efgh=0;efgh<4;efgh++){
                    overlaps=true;
                    while (overlaps){
                        xipos = tmp[abcd] + ran_norm();
                        xjpos = tmp[efgh] + ran_norm();
                        overlaps= testTracerOverlap(i, j, xipos, xjpos);
                        //cout << "Repeat?";
                    }
                    CRod newRod = CRod(axis, xipos, xjpos, _ranU, m_igen );
                    _drods[axis][abcd][efgh] = newRod;
                }
            }
        }
        ifdebug(prinRodPos(0);)
    }
    
    void updateRand(int crossaxis,int exitmarker){
        //exitmarker is -1 for negative direction, or 1 for positive
        //delete all polymers orthogonal to crossaxis, that are outside the box now
        //update other polymer positions
        bool overlaps;
        double cellInterval_ai, cellInterval_aj;
        int i,j;
        i=crossaxis+1;
        if (i==3) i =0;
        j=3-(i+crossaxis);
        // rotate around rods in cells abc and def and reassign
        if (exitmarker == 1){
            for (int abcd=0; abcd<4;abcd++){
                for (int efgh=0; efgh<4;efgh++){
                    _drods[i][abcd][efgh].coord[crossaxis] -= _bdef;
                    _drods[j][abcd][efgh].coord[crossaxis] -= _bdef;
                }
            }
            rotate_left(_drods[j]);
            cellInterval_ai = - _bdef;
            cellInterval_aj = - _bdef;
            for (int abcd=0;abcd<4;abcd++){
                rotate_left(_drods[i][abcd]);
                overlaps=true;
                while (overlaps){
                    _drods[i][abcd][3].coord[crossaxis] = 2*_bdef + ran_norm();
                    _drods[i][abcd][3].coord[j] = cellInterval_aj + ran_norm();
                    overlaps= testTracerOverlap(crossaxis, j, _drods[i][abcd][3].coord[crossaxis], _drods[i][abcd][3].coord[j])
                        //TODO overlap
                        || testRodOverlap(i, crossaxis, j, _drods[i][abcd][3].coord[crossaxis], _drods[i][abcd][3].coord[j]);
                    //cout << "Repeat?";
                }
                overlaps=true;
                while (overlaps){
                    _drods[j][3][abcd].coord[crossaxis] = 2*_bdef + ran_norm();
                    _drods[j][3][abcd].coord[i] = cellInterval_ai + ran_norm();
                    overlaps= testTracerOverlap(crossaxis, i, _drods[j][3][abcd].coord[crossaxis], _drods[j][3][abcd].coord[i])
                    //TODO overlap
                        || testRodOverlap(j, crossaxis, i, _drods[j][3][abcd].coord[crossaxis], _drods[j][3][abcd].coord[i]);
                }
                cellInterval_aj+=_bdef;
                cellInterval_ai+=_bdef;
            }
        }
        else{
            // shift positions of rods
            for (int abcd=0; abcd<4;abcd++){
                for (int efgh=0; efgh<4;efgh++){
                    _drods[i][abcd][efgh].coord[crossaxis] += _bdef;
                    _drods[j][abcd][efgh].coord[crossaxis] += _bdef;
                }
            }
            rotate_right(_drods[j]);
            cellInterval_ai = - _bdef;
            cellInterval_aj = - _bdef;
            for (int abcd=0;abcd<4;abcd++){
                rotate_right(_drods[i][abcd]);
                // new rod positions
                overlaps=true;
                while (overlaps){
                    _drods[i][abcd][0].coord[crossaxis] = -_bdef + ran_norm();
                    _drods[i][abcd][0].coord[j] = cellInterval_aj + ran_norm();
                    overlaps= testTracerOverlap(crossaxis, j, _drods[i][abcd][0].coord[crossaxis], _drods[i][abcd][0].coord[j])
                        //TODO overlap
                        || testRodOverlap(i, crossaxis, j, _drods[i][abcd][0].coord[crossaxis], _drods[i][abcd][0].coord[j]);
                }
                overlaps=true;
                while (overlaps){
                    _drods[j][0][abcd].coord[crossaxis] = -_bdef + ran_norm();
                    _drods[j][0][abcd].coord[i] = cellInterval_aj + ran_norm();
                    overlaps= testTracerOverlap(crossaxis, i, _drods[j][0][abcd].coord[crossaxis], _drods[j][0][abcd].coord[i])
                        //TODO overlap
                        || testRodOverlap(j, crossaxis, i, _drods[j][0][abcd].coord[crossaxis], _drods[j][0][abcd].coord[i]);
                }
                cellInterval_aj+=_bdef;
                cellInterval_ai+=_bdef;
            }
        }
        ifdebug(
            if (testOverlap()){
                cout << "\nERROR still overlap after newrod init!" << endl;
                //abort();
            }
        )
    }


    bool testTracerOverlap(int i, int j, double ri, double rj){
        for (auto & bead :  _beads){
            //make this a little bigger, so that it's out of reach of LJ
            if ((pow( bead.pos(i) - ri , 2 ) + pow( bead.pos(j) - rj , 2 )) < 1.13*_cylLJSq){
                return true;
            }
        }
        return false;
    }

    // bool testTracerOverlap(int i, int j, double ri, double rj){
//         return ((pow( _ppos(i) - ri , 2 ) + pow( _ppos(j) - rj , 2 )) < 1.13*_cylLJSq);//make this a little bigger, so that it's out of reach of LJ
//     }
    
    //TODO overlap
    bool testRodOverlap(int rodaxis, int i, int j, double ri, double rj){
        double distSq;
        if (_polyrad==0) return false;
        for (auto & vrods : _drods[rodaxis]){
            for (auto & rod :  vrods){
                //rods axis need to have distance of a least a (i.e. polydiam)
                distSq = pow(rod.coord[i] - ri,2) + pow(rod.coord[j] - rj,2);
                if (  (distSq < _polydiamSq) && (distSq > 0.00001)  ){//The second clause is to avoid testing overlap with itself
                    return true;
                }
            }
        }
        return false;
    }

    void prinRodPos(int axis){
        for (int irod=0;irod<_drods[axis].size();irod++){
            for (int jrod=0;jrod<_drods[axis].size();jrod++){
                double rx =_drods[axis][irod][jrod].coord[0];
                double ry =_drods[axis][irod][jrod].coord[1];
                double rz =_drods[axis][irod][jrod].coord[2];
                cout << ",[" << rx << "," << ry << "," << rz << "]";
            }
        }
        cout << "]," << endl;
    }

    bool tracerInBox(){
        for (int ax = 0;ax<3;ax++){
            if ((_ppos(ax) < -0.1*_boxsize[ax]) || (_ppos(ax) > 1.1*_boxsize[ax])){
                cout << "\nax " << ax << "\n_boxsize[ax] " << _boxsize[ax] << "\n_ppos(ax) " << _ppos(ax) << endl;
                return false;
            }
        }
        return true;
    }




    template<typename T, size_t N>
    void rotate_left(std::array<T,N> & arr){
        // rotation to the left
        std::rotate(arr.begin(), arr.begin() + 1, arr.end());
    }
    template<typename T, size_t N>
    void rotate_right(std::array<T,N> & arr){
        // rotation to the  right
        std::rotate(arr.rbegin(), arr.rbegin() + 1, arr.rend());
    }



private:
    void setRanNumberGen(double seed);
    void countWallCrossing(int crossaxis, int exitmarker);
    void calculateExpHPI(const double r, double& U, double& Fr);
    void calculateExpPotential(const double r, double& U, double& Fr, int sign);
    void calculateExpPotentialMOD(const double r, double& U, double& Fr, int index);
    void modifyPot(double& U, double& Fr, double dist);
    void calcLJPot(const double r, double &U, double &dU);
    void initPosHisto();





public:
    CConfiguration();
    CConfiguration(string trigger, double timestep,  double potRange,  double potStrength, 
        double psize, const bool posHisto, const bool steric, const bool ranU, double dvar, double polydiam, string peptide, double uDebye, double uBend);
    void updateStartpos();
    int makeStep();
    void checkBoxCrossing();
    void calcStochasticForces();
    void calcMobilityForces();
    void calcBeadInteraction();
    void calc_1RODONLY_MobilityForces();
    void calc_ZRODONLY_MobilityForces();
    void calc_YZRODONLY_MobilityForces();
    void saveXYZTraj(string name, int move, string flag);

    vector<double> getUpot(){
        vector<double> upot_vec(_N_beads);
        for (int i=0; i<_N_beads; i++) {
             upot_vec[i] = _beads[i].upot;
        }
        return upot_vec; 
    }
    int get_Nbeads(){
        return _N_beads;
    }    
    
    double getDisplacement();
    unsigned int getwallcrossings(int i){ return _wallcrossings[i]; }
    bool testOverlap();
    void moveBack();
    //void addHistoValue();
    //void printHistoMatrix(string folder);
    //void positionHistogram(double block, double possq, double pposprev, int l, int posHisto[]);
    std::vector<double> getppos(){ // returns pointer to current particle position array
    	std::vector<double> pos (3);
    	for (int i = 0; i < 3; i++){
            pos[i] = _ppos(i) + _boxCoord[i];
    	}
    	return pos;
    }
    std::vector<double> getppos_rel(){ // returns pointer to current particle position array
    	std::vector<double> pos (3);
    	for (int i = 0; i < 3; i++){
            pos[i] = _ppos(i);
    	}
    	return pos;
    }
    void writeDistances(ostream& distancesfile, unsigned int stepcount) {
    // So far this only writes the tracer particle position
	distancesfile << fixed << stepcount * _timestep << "\t";
        for (auto & arr : _distarr){
            for (auto & dist : arr){
                distancesfile << dist << " ";
            }
        }
        distancesfile << endl;
    }

    void saveBeadsTraj(ofstream &beadstrajfile, double simtime){
        beadstrajfile << fixed << simtime;
        for (auto & bead :  _beads){
            beadstrajfile << fixed << " \t" << bead.pos(0) << " " << bead.pos(1) << " " << bead.pos(2);
        }
        beadstrajfile << endl;
    }



};



#endif /* CCONFIGURATION_H_ */
