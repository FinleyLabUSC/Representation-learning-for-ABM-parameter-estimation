#include "Environment.h"

Environment::Environment(std::string saveFld): mt((std::random_device())()) {
    /*
     * initialize a simulation environment
     * -----------------------------------
     * loads the three parameter files
     *  cell parameters
     *  environment parameters
     *  recruitment parameters
     *
     * set environment variables to their respective values
     */

    saveDir = saveFld;
    loadParams();

    cd8RecRate = recParams[0];
    cd8Ratio = recParams[1];
    recDist = recParams[2];

    if(envParams[1] == 1){
        threeD = 1.0;
    } else{
        threeD = 0.0;
    }
    simulationDuration = envParams[0];

    tumorCenter = {0,0,0};
    tumorRadius = 0;

    steps = 0;

    dt = 0.005;
    cd82rec = 0;
}

void Environment::simulate(double tstep) {
    /*
     * initializes and runs a simulation
     * ---------------------------------
     * place initial tumor
     * run simulation loop
     *  recruit immune cells
     *  run cell functions
     * ends once time limit is reached or there are no more cancer cells
     */

    //cell_list.push_back(Cell({0,0,0}, 0, cellParams, "cancer", threeD));
    int radiiCells = 5;
    cell_list.push_back(Cell({0,0,0}, 0, cellParams, "cancer", threeD, 0.0));
    int q = 1;
    for(int i=1; i<radiiCells; ++i){
        double circumfrence = 2*i*cellParams[8][0]*3.1415;
        double nCells = circumfrence/cellParams[8][0];
        for(int j=0; j<nCells; ++j){
            double x = i*cellParams[8][0]*cos(2*3.1415*j/nCells);
            double y = i*cellParams[8][0]*sin(2*3.1415*j/nCells);
            cell_list.push_back(Cell({x,y,0}, q, cellParams, "cancer", threeD, 0.0));
            q++;
        }
    }

    tumorSize();

    std::cout << "starting simulation...\n";
    while(tstep*steps/24 < simulationDuration) {
        recruitImmuneCells(tstep);
        runCells(tstep);

        if (fmod(steps * tstep, 24) == 0) {
            // update every simulation day
            tumorSize();
        }

        steps += 1;
        printStep(steps * tstep);
        if (fmod(steps * tstep, 24) == 0) {
            // save every simulation day
            save(tstep);
        }

        int numC = 0;
        for (auto &c: cell_list) {
            if (c.type == 0) {
                numC++;
            }
        }
        if (numC == 0) {
            break;
        }
    }
    tumorSize();
    save(tstep);
}
