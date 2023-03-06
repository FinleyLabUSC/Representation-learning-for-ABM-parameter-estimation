#include "Cell.h"

void Cell::initializeCD8Cell(std::vector<std::vector<double> > &cellParams) {
    state = 1;

    double diameter = cellParams[11][1];

    mu = cellParams[0][1];
    kc = cellParams[1][1];
    damping = cellParams[2][1];
    maxOverlap = cellParams[3][1]*diameter;
    deathProb = cellParams[4][1];
    migrationSpeed = cellParams[5][1];
    killProb = cellParams[6][1];
    influenceRadius = cellParams[7][1];
    std::normal_distribution<double> infilDist(0.0,cellParams[8][1]/3);
    //std::uniform_real_distribution<double> infilDist(0.0, cellParams[8][1]);
    infiltrationDistance = fabs(infilDist(mt));
    migrationBias = cellParams[9][1];
    influenceDec = cellParams[10][1];
    radius = diameter/2.0;

    rmax = 1.5*diameter;
}

void Cell::pdl1Inhibition(std::array<double, 3> otherX, double otherRadius, double otherpdl1, double dt) {
    // inhibition via direct contact

    if(type != 1){return;}
    if(state == 2){return;}

    double distance = calcDistance(otherX);
    if(distance <= radius+otherRadius){
        std::uniform_real_distribution<double> dis(0.0,1.0);
        if(dis(mt) < probTime(otherpdl1, dt)){
            state = 2;
            killProb = 0;
            influenceRadius *= influenceDec;
            migrationSpeed = 0.0;
        }
    }
}
