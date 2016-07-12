#include "headers/CConfiguration.h"


using namespace std;
using namespace Eigen;


CConfiguration::CConfiguration(){
}

CConfiguration::CConfiguration(double timestep,  double potRange,  double potStrength,
        double psize, const bool posHisto, const bool steric, const bool ranU, double dvar, double polydiam, string peptide, double uDebye){
    setRanNumberGen(0);
    _potRange = potRange;
    _potStrength = potStrength;
    _pradius = psize/2.;
    //TODO overlap
    _polyrad = polydiam/2.;//TODO test
    _polydiamSq = polydiam*polydiam;
    _cylLJSq = pow(_pradius + _polyrad,2);
    _cutoffExpSq = pow(8*_potRange,2);
    _timestep = timestep;
    _LJPot = (steric == false) && (psize != 0);
    _ranU = ranU;
    _dvar = dvar;
    _uDebye = uDebye;
    _mu_sto = sqrt( 2 * _timestep );                 //timestep for stochastic force
    //Spring interaction
    _r0SP = 1.122462 * 2*_pradius;
    for (int i = 0; i < 3; i++){
        _ppos(i) = _bdef/2.;
        _startpos[i] = _ppos(i);
        _entryside[i] = 0;
        _wallcrossings[i] = 0;
        _boxCoord[i] = 0;
        _prevpos(i) = _ppos(i);
    }


    cout << "init Rand..." << endl;
    initRand();

    initPeptide(peptide);

}

void CConfiguration::updateStartpos(){
    //This function is used if the particle should keep moving after each run, and not start at _resetpos again, like in the first run
    //This (hopefully) will give better averages without having to spend a lot of steps in the beginning of each run to get away from _resetpos
    for (int i = 0; i < 3; i++){
    _startpos[i] = _ppos(i) + _boxCoord[i];
    }
}




//TODO pep
void CConfiguration::makeStep(){
    //move the particle according to the forces and record trajectory like watched by outsider
    _prevpos = _beads[0].pos;
    for (auto & bead :  _beads){
        for (int i = 0; i < 3; i++){
            const double disp = _timestep * bead.f_mob(i) + _mu_sto * bead.f_sto(i);
            bead.pos(i) += disp;
            if (abs(disp) > 5){
                cout << "**** Way too big jump!\n"
                        << "\nuDebye = " << _uDebye 
                            << " -- uspring = " << _uSpring << " -- uLJ beads = " << _uLJ <<  endl;
            }
            else if (std::isnan(bead.pos(i))){
                cout << "axis " << i << "\n_prevpos " << _prevpos(i) << "\n_ppos " << bead.pos(i) << endl;
                cout << "f_mob[i] " << bead.f_mob(i) << "\nf_sto[i] " << bead.f_sto(i) << endl;
                abort();
            }
        }
    }
    _ppos = _beads[0].pos;
}

void CConfiguration::checkBoxCrossing(){
    //should the particle cross the confinement of the cube, let it appear on the opposite side of the box
    int exitmarker = 0;
    for (int i = 0; i < 3; i++){
        exitmarker =0;
        if (_ppos(i) < -0.05*_boxsize[i]){//Only create new rod config if the particle has crossed the border by a certain fraction
            _ppos(i) += _boxsize[i];
            _boxCoord[i] -= _boxsize[i];
            exitmarker = -1;
        }
        else if (_ppos(i) > 1.05*_boxsize[i]){//Only create new rod config if the particle has crossed the border by a certain fraction
            _ppos(i) -= _boxsize[i];
            _boxCoord[i] += _boxsize[i];
            exitmarker = 1;
            if (_ppos(i) > 10){
                cout << "Bad _ppos after boxcrossing. Aborting!" << endl;
                cout << "_ppos b4 " << _ppos(i) + _boxsize[i] << endl;
                cout << "_ppos after " << _ppos(i) << endl;
                cout << "\naxis " << i << "\n_boxsize[ax] " << _boxsize[i] << endl;
                abort();
            }
        }
        if (exitmarker!=0){
            updateRand(i,exitmarker);
            ifdebug(
                cout << "[["<<i<<"," << exitmarker <<"] ";
                prinRodPos(0); // cout print rod pos!
            )
            if (_ranU){
                for (auto & vrods : _drods[i]){
                    for (auto & rod :  vrods){
                        rod.shiftSigns(exitmarker);
                    }
                }
            }
            //shift bead positions
            for (auto & bead :  _beads){
                bead.pos(i) -= exitmarker * _boxsize[i];
            }
        }
    }
    //TODO del
    if (!tracerInBox()){
        cout << "\nTRACER NOT IN BOX!!!" << endl;
        abort();
    }
}






void CConfiguration::calcStochasticForces(){

    // the variate generator uses m_igen (int rand number generator),
    // samples from normal distribution with standard deviation 1 (later sqrt(2) is multiplied)
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > ran_gen(
            *m_igen, boost::normal_distribution<double>(0, 1));

    for (auto & bead :  _beads){
        for (int i = 0; i < 3; i++) {
            bead.f_sto(i) = ran_gen();
        }
    }
}


void CConfiguration::calcMobilityForces(){
    //calculate mobility forces from potential Epot - Unified version that includes 2ndOrder if k is larger than or equal 0.2 b , except if ranPot is activated.
    double Epot = 0;
    double rcSq = 1.25992 * _cylLJSq;
    for (auto & bead :  _beads){
        bead.upot = 0.;
        double r_i = 0, r_k = 0;
        array<double,4> r_is, r_ks;
        std::array<double, 16> ri_arr, rk_arr, rSq_arr;
        double r_absSq;
        double utmp = 0, frtmp = 0, uCylTot = 0;     //temporary "hilfsvariables"
        //reset mobility forces to zero
        bead.f_mob = Vector3d::Zero();
    

        for (int i = 0; i < 3; i++){
            unsigned int cnt=0;// counter to loop over array indices
            int k = i + 1;   //k always one direction "further", i.e. if i = 0 = x-direction, then k = 1 = y-direction
            if ( k == 3 ) k = 0;
            int plane = 3 - (i+k); //this is the current plane of the cylindrical coordinates
            int n = 0;     // reset counter for index of next rod in plane  n = 0, 1, 2, 3 -> only needed for ranPot
            double z1, z2, z1inv;
            if (_ranU){
                z1 = 0.25 * _boxsize[plane];
                z2 = _boxsize[plane] - z1;   //z is in cylindrical coordinates. This indicates above/below which value the exp potential is modifed for random signs.
                z1inv = 1./z1;
            }

            for (int abcd=0;abcd<4;abcd++){
                for (int efgh=0;efgh<4;efgh++){
                    r_i = bead.pos(i) - _drods[plane][abcd][efgh].coord(i);
                    r_k = bead.pos(k) - _drods[plane][abcd][efgh].coord(k);
                    ri_arr[cnt]=(r_i);
                    rk_arr[cnt]=(r_k);
                    rSq_arr[cnt]=( r_i * r_i + r_k * r_k);
                    cnt++;
                }
            }
            //TODO distarr[i] = rSq_arr; //store distances for writing them to a file later
            for (int j=0;j<cnt;j++){
                const double rSq = rSq_arr.at(j);
                calculateExpPotential(rSq, utmp, frtmp, bead.sign);

            
                if (_ranU){
                    int abcd = j/4; // This could be made more efficient by replacing the j loop with two loops abcd efgh and a counter j
                    int efgh = j%4;
                    int sign = _drods[plane][abcd][efgh].signs[1];
                    //cout << "abcd " << abcd << "efgh " << efgh << " sign: " << sign << endl;
                    utmp *= sign;
                    frtmp *= sign;
                    if (bead.pos(plane) > z2){
                        if (! _drods[plane][abcd][efgh].samesign[1]){
                            bead.f_mob(plane) += utmp * z1inv;              //this takes care of the derivative of the potential modification and resulting force
                            modifyPot(utmp, frtmp, (_boxsize[plane] - bead.pos(plane)) * z1inv);
                        }
                    }
                    else if (bead.pos(plane) < z1){
                        if (! _drods[plane][abcd][efgh].samesign[0]){
                            bead.f_mob(plane) -= utmp * z1inv;              //this takes care of the derivative of the potential modification and resulting force
                            modifyPot(utmp, frtmp, bead.pos(plane) * z1inv);
                        }
                    }
                }

                if (_LJPot && ( rSq < rcSq )) calcLJPot(rSq, utmp, frtmp, _cylLJSq);

                //TODO del
                if (utmp > 100){
                    cout << "utmp " << utmp  <<"\nr " << sqrt(rSq) << "\nindex j "  << j << "\nplane " << plane << endl;
                    cout << "ri " << ri_arr[j] << "\nrk " << rk_arr[j] << endl;
                }

                uCylTot += utmp;
                bead.upot += utmp;
                bead.f_mob(i) += frtmp * ri_arr[j];
                bead.f_mob(k) += frtmp * rk_arr[j];
            }
        }
        //ifdebug(cout << "bead.f_mob \n" << bead.f_mob << endl;)
        ifdebug(cout << "uCylTot " << uCylTot <<  endl; )
    }
    
    if (_N_beads != 1) calcBeadInteraction();
}


void CConfiguration::calcBeadInteraction(){//TODO pep
    const double beadLJSq = 4*_pradius*_pradius;
    _uDebye = 0; _uLJ = 0; _uSpring = 0;
    int cnt=0;
    Vector3d faddtmp;
    // Eigen::Vector3d _rij_vec_arr;
    // vector<double> _rijSq_arr;
    // for (int i=0; i < _N_beads; i++){// Test if adding this loop to the next one increases speed
//         for (int j=i+1; j < _N_beads; j++){
//             _rij_vec_arr[cnt] = _beads[j].pos - _beads[i].pos;
//             _rijSq_arr[cnt] = _rij_vec_arr[cnt].squaredNorm();
//
//             cnt++;
//         }
//     }
    // loop to calc potentials and forces
    for (int i=0; i < _N_beads; i++){
        for (int j=i+1; j < _N_beads; j++){
            double frtmp = 0;
            double utmp = 0;
            // const double rijSq = _rijSq_arr[cnt];
            // const Vector3d rvec = _rij_vec_arr[cnt];
            const Vector3d rvec = _beads[j].pos - _beads[i].pos;
            const double rijSq = rvec.squaredNorm();
            double rij = 0; 
            int debyeTrig = _beads[i].sign * _beads[j].sign;
            
            // LENNARD-JONES
            if ( rijSq <  1.25992 * beadLJSq ) {
                calcLJPot(rijSq, utmp, frtmp, beadLJSq);
                _uLJ += utmp;
            }
            
            // DEBYE
            if (debyeTrig != 0){// only calculate debye if both particles are charged.
                rij = sqrt(rijSq);
                calcDebyePot(rij, utmp, frtmp, debyeTrig);
                _uDebye += utmp;
            }
            
            // SPRING
            if (j==i+1){
                // only calc square root if it has not already been calculated before for Debye
                if (rij==0) rij = sqrt(rijSq);
                calcSpringPot(rij, utmp, frtmp);
                _uSpring += utmp;
            }
            faddtmp = frtmp * rvec;
            _beads[i].upot += utmp;
            _beads[j].upot += utmp;
            _beads[i].f_mob += - faddtmp ;
            _beads[j].f_mob += faddtmp;
            
            // This needs to come last!!!
            cnt++;
        }
        //ifdebug(cout << "beadInteraction: bead[i].f_mob \n" << _beads[i].f_mob << endl;)
        ifdebug(cout<< "Energies:\n"  <<
             "\nuDebye = " << _uDebye << " -- uspring = " << _uSpring << " -- uLJ beads = " << _uLJ <<  endl; )
    }
}

// So far this guy only writes the positions of the beads along the peptide
void CConfiguration::saveXYZTraj(string name, int move, string flag) {
    Eigen::Map<Eigen::Vector3d> boxCoordinates(_boxCoord);
    Vector3d rtmp;
    if(flag=="w") {    //write to new file
        /*if(m_traj_file!=NULL) {
            fclose(m_traj_file);
        }*/
        m_traj_file = fopen(name.c_str(), flag.c_str());
        if(m_traj_file==NULL) {
            cout << "error creating trajfile" << endl;
        }
    }

    fprintf(m_traj_file, "%d\n%s (%8.3f %8.3f %8.3f) t=%u \n", _N_beads, "sim_name", _boxsize[0], _boxsize[1], _boxsize[1], move);


    // Beads
    for (auto & bead :  _beads) {
        rtmp = bead.pos;//+boxCoordinates;
        fprintf(m_traj_file, "%3s%9.3f%9.3f%9.3f \n","H", rtmp(0), rtmp(1),  rtmp(2));
    }

    //fflush(m_traj_file);

    if(flag=="c") {    //close file
        if(m_traj_file!=NULL) { fclose(m_traj_file); }
    }
}









void CConfiguration::setRanNumberGen(double seed){
    if (seed == 0) {
        m_igen = new boost::mt19937(static_cast<unsigned int>(time(NULL)));
        cout << "random seed is time!" << endl;
    } else {
        m_igen = new boost::mt19937(static_cast<unsigned int>(seed));
        cout << "random seed is " << seed << endl;
    }
}



void CConfiguration::moveBack(){
    //moves particle back to previous position
    _ppos = _prevpos;
}
//****************************POTENTIALS**********************************************************//



void CConfiguration::calculateExpPotential(const double rSq, double& U, double& Fr, int sign){
    //function to calculate an exponential Potential U = U_0 * exp(-1 * r * k)
    // k is the interaction range. U_0 is the strength of the potential
    //which is attractive if direction = -1, and repulsive if direction = 1
    //The potential is weighted with kT!
    if (rSq < _cutoffExpSq && sign != 0){
        const double r = sqrt(rSq);
        U = sign * _potStrength * exp(-1.0 * r / _potRange);
        Fr = U / (_potRange * r);  //This is the force divided by the distance to the rod!
    }
    else{
        U=0;
        Fr=0;
    }
}


void CConfiguration::calculateExpHPI(const double r, double& U, double& Fr){
    cout << "\n\n....calculateExpHPI: NOTHING HERE ..." << endl;
    abort();
//	double u = _hpi_u * exp( - r / _hpi_k);
//	U += u;
//	Fr += u / (_hpi_k * r);
}


void CConfiguration::calculateExpPotentialMOD(const double r, double& U, double& Fr, int plane){
    cout << "NOT DEFINED!!!!!!!" << endl;
}

void CConfiguration::modifyPot(double& U, double& Fr, double weight){
    //function to modify the potential according to the distance along the polymer axis to the next neighbor,
    //in case the next neighboring polymer part is of opposite sign
    //the weight is determined by the distance to the point where the sign changes, divided by the boxsize: weight = 4 * dist / boxsize, such that it's 1 for dist=boxsize/4
//    cout << "NOT DEFINED" << endl;
    U *= weight;
    Fr *= weight;
}

//****************************STERIC HINDRANCE****************************************************//

bool CConfiguration::testOverlap(){//TODO define for 
    //Function to check, whether the diffusing particle of size psize is overlapping with any one of the rods (edges of the box)
    //most if borrowed from moveParticleAndWatch()
    cout << "ERROR! Not defined yet !!!!" << endl;
    abort();

    bool overlaps = false;
    double r_i = 0, r_k = 0;
    double r_abs = 0;

    int i,j;
    for (int axis=0;axis<3;axis++){
        for (int abc=0;abc<_drods[axis].size();abc++){
            for (int def=0;def<_drods[axis].size();def++){
                i = axis +1;
                if (i==3) i=0;
                j=3-(i+axis);
                if (testTracerOverlap(i, j, _drods[axis][abc][def].coord(i), _drods[axis][abc][def].coord(j))){
                    ifdebug(cout << "Overlap for\naxis" << axis << "\nabc " << abc << "\ndef " << def << endl;)
                    return true;
                }
            }
        }
    }
    
    return overlaps;
}





//****************************OLD CODE****************************************************

//****************************POS HISTOGRAM****************************************************//

// void CConfiguration::initPosHisto(){
//     _posHistoM.resize(100);
//     for (int i = 0; i < 100; i++){
//         _posHistoM[i].resize(100);
//         for (int j = 0; j < 100; j++){
//             _posHistoM[i][j].resize(100, 0);  //initialize all the 100*100*100 matrix elements to zero!
//         }
//     }
// }
//
// void CConfiguration::addHistoValue(){
//     //adds a value to the position histogram
//     int x = _ppos(0) / _boxsize * 100;     //CAREFUL: THIS CAN'T BE DONE AT A POINT WHERE X MIGHT BE ZERO!!!
//     int y = _ppos(1) / _boxsize * 100;
//     int z = _ppos(2) / _boxsize * 100;
//     if ((x < 0) || (y < 0) || (z < 0) || (x > 99) || (y > 99) || (z > 99)){
//         cout << "The Position Histogram function 'conf.addHisto()' is in a bad place, since there is a negative position _ppos()" << endl;
//         cout << "The current position is: " << _ppos(0) << " " << _ppos(1) << " " << _ppos(2) << endl;
//     }
//     _posHistoM[x][y][z] += 1;
// }
//
// void CConfiguration::printHistoMatrix(string folder){
//     //function to print the positionHistogram to a file called posHistoMatrix.txt
//     //The elements of the matrix are M[x][y][z]. First the z is counted from 1 to 100 in one row, then the y follows, then, after the 'X', the 100 x elements.
//
//     ofstream matrixfile;
//     matrixfile.open((folder + "/InstantValues/posHistoMatrix.txt").c_str());
//     int maxval = 0;
//
//     for (int i = 0; i < 100; i++){
//         for (int j = 0; j < 100; j++){
//             for (int k = 0; k < 100; k++){
//                 matrixfile << _posHistoM[i][j][k] << " ";
//                 if (maxval < _posHistoM[i][j][k] ) maxval = _posHistoM[i][j][k];
//             }
//             matrixfile << endl;
//         }
//         matrixfile << "X" << endl;
//     }
//     matrixfile << "MaxVal " << maxval;   //THis does not affect the grid files, since data is only copied to them when a "X" line comes!
//
//
//     matrixfile.close();
// }



void CConfiguration::calc_1RODONLY_MobilityForces(){   //see OldCode
}



void CConfiguration::calc_ZRODONLY_MobilityForces(){   //see OldCode
}



void CConfiguration::calc_YZRODONLY_MobilityForces(){   //see OldCode
}

double CConfiguration::getDisplacement(){
    double d = 0;
    for (int i = 0; i < 3; i++){
        d += pow((_timestep * _f_mob[i] + _mu_sto * _f_sto[i]), 2);
    }
    return sqrt(d);
}

