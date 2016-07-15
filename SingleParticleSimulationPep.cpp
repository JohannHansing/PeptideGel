#include "headers/SingleParticleSimulationPep.h"

using namespace std;



int main(int argc, const char* argv[]){
    //Main includes the iteration loops for the simulation

    //NOTE: so far wallcrossings is added for all runs!!! makes it kind of incorrect, since after each run, ppos is reset.
    //NOTE: so far saving Instant Values for each tenth step!

    //TRIGGERS:
    string trigger = argv[1];
    bool ranRod = (strcmp(argv[2] , "true") == 0 ) ;
    bool rand = (strcmp(argv[3] , "true") == 0 ) ;
    bool recordMFP = (strcmp(argv[4] , "true") == 0 ) ;
    bool recordPosHisto = (strcmp(argv[5] , "true") == 0 ) ;
    bool includeSteric = (strcmp(argv[6] , "true") == 0 ) ;  // steric 2
    bool ranU = (strcmp(argv[7] , "true") == 0 ) ;
    string peptide = argv[8];
    int boolpar = 8;
    ifdebug(cout << "copied bools. ";)

    // Checking for correct structure of input arguments
    for (int k= 1; k < argc; k++ ) cout << "parameter " << k << " " << argv[k] << endl;
    for (int b_i=2; b_i<boolpar; b_i++){
        if (!((strcmp(argv[b_i] , "true") == 0 )  || (strcmp(argv[b_i] , "false") == 0 ))){
            cerr << "Error; Bool parameter " << b_i << " is not either 'true' or 'false'!" << endl;
            exit(1);
        }
    }

    int runs = atoi( argv[boolpar+1] );                       // Number of Simulation runs to get mean values from
    double timestep = atof( argv[boolpar+2] );
    int simtime = atoi( argv[boolpar+3] );                   // simulation time
    int instantvalues = 200;
    unsigned int steps;

    double particlesize = atof( argv[boolpar+4] );
    double urange = atof( argv[boolpar+5] );
    double ustrength = atof( argv[boolpar+6] );
    double dvar = atof( argv[boolpar+7] );
    double polydiam = atof( argv[boolpar+8] );
    double uDebye = atof( argv[boolpar+9] );
    unsigned int saveInt;
    int instValCount = 0;                             //Counter for addInstantValue
    int badsteps = 0;                                 //Count how many times bead displacement had to be adjusted

    steps = simtime/timestep;
    saveInt = steps/instantvalues;
    const int trajout = (int)(10/timestep);
    

    ifdebug(cout << "copied  params. ";)

    cout << "--- PEPTIDE is " << peptide << endl;
    if (trigger!="noCyl" && trigger!="x"){
        cout << "ERROR: Bad trigger!" << endl;
        abort();
    }

    //initialize instance of configuration
    CConfiguration conf = CConfiguration(trigger, timestep, urange, ustrength, particlesize, recordPosHisto, 
                            includeSteric, ranU, dvar,polydiam, peptide, uDebye);
    ifdebug(cout << "created CConf conf. ";)
    

    //Create data folders and print location as string to string "folder"
    string folder = createDataFolder(trigger, timestep, simtime, urange, ustrength, particlesize, includeSteric, ranU, 
                             dvar,polydiam, peptide, uDebye);
    ifdebug(cout << "created folder. ";)
    cout << "writing to folder " << folder << endl;


    //initialize averages
    CAverage energyU = CAverage("Upot", folder, conf.get_Nbeads());
    ifdebug(cout << "created CAverage files. ";)


    unsigned int stepcount = 0;
    ofstream trajectoryfile;
    trajectoryfile.open((folder + "/Coordinates/trajectory.txt").c_str());
    
    ofstream beadstrajfile;
    beadstrajfile.open((folder + "/Coordinates/beadstraj.txt").c_str());
    beadstrajfile << "# Trajectories for peptide " << peptide << endl;
    
    
    ofstream distancesfile;
    // TODO distancefile
    // NEED TO DEFINE BEHAVIOR OF SAVING DISTANCES IN CCONF
    //distancesfile.open((folder + "/Coordinates/squareDistances.txt").c_str());
    
    settingsFile(folder, trigger, particlesize, timestep, runs, steps, ustrength, urange, recordMFP, includeSteric, ranU,  
                    dvar, polydiam, peptide, uDebye, badsteps);
                    
    //create .xyz file to save the trajectory for VMD
    string traj_file = folder + "/Coordinates/single_traj.xyz";
    conf.saveXYZTraj(traj_file,0,"w");
    


    //cout << "Starting Run Number: " << simcounter << " out of " << totalsims << endl;
    cout << "Starting Simulation!" << endl;


// **************START OF RUNS-LOOP
    for (int l = 0; l<runs; l++){

        conf.updateStartpos();

        for (int i = 0; i < steps; i++){  //calculate stochastic force first, then mobility force!!


            conf.calcStochasticForces();

            conf.calcMobilityForces();
            conf.calcBeadInteraction();


            if (((i+1)%saveInt) == 0){       //saving Instant Values for each saveInt'th step!
                energyU.addInstantValue(conf.getUpot());
                instValCount += 1;
            }

            //move particles
            badsteps += conf.makeStep();

            /*
            if (includeSteric && conf.testOverlap()) conf.moveBack();   //TODO steric2
            else conf.checkBoxCrossing();
            */



                //TODO steric
            while (includeSteric && conf.testOverlap()){
                conf.moveBack();
                conf.calcStochasticForces();
                conf.makeStep();
            }
            conf.checkBoxCrossing(); //check if particle has crossed the confinement of the box

            stepcount++;
            if (stepcount%trajout == 0) {
                std::vector<double> ppos = conf.getppos();
                trajectoryfile << fixed << stepcount * timestep << "\t" << ppos[0] << " " << ppos[1] << " " << ppos[2] << endl;
                ifdebug(cout << stepcount * timestep << "\t" << ppos[0] << " " << ppos[1] << " " << ppos[2] << endl;)
                //TODO pass distancefile to function in conf.
                //if (stepcount%(100*trajout) == 0) conf.writeDistances( distancesfile, stepcount);
                //if (stepcount%(10*trajout) == 0){
                    conf.saveBeadsTraj(beadstrajfile, stepcount * timestep);
                //}
            }
            
            if (((i+1)%100 == 0) && (l == 0)){       //Save the first trajectory to file
                //TODO XYZ traj
                //conf.saveXYZTraj(traj_file, i, "a");                    // TODO change back ((i+1)%XXX == 0) to 100
            }
        }
        if (l==0) conf.saveXYZTraj(traj_file, steps, "c"); // Close XYZ traj_file
        
        if (l%50 == 0){
            cout << "run " << l << endl;
            energyU.saveAverageInstantValues(instValCount);
            settingsFile(folder, trigger, particlesize, timestep, runs, steps, ustrength, urange, recordMFP, includeSteric, ranU,  
                            dvar, polydiam, peptide, uDebye, badsteps);
        }
        
        if (badsteps>10000){//if about one out of 100,000 steps is registered as badstep, assuming that the code is running for at least 500 runs
            ofstream warnfile;
            warnfile.open((folder + "/WARNING.txt").c_str());
            warnfile << "High number of badsteps: " << badsteps << endl;
            warnfile.close();
        }


    }//----------END OF RUNS-LOOP


    cout << "Simulation Finished" << endl;

    trajectoryfile.close();
    beadstrajfile.close();
    // TODO distancefile
    //distancesfile.close();

    return 0;
}
