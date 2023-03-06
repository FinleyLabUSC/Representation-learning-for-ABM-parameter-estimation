#include "Cell.h"

void Cell::initializeCancerCell(std::vector<std::vector<double>> &cellParams) {
    state = 0;
    canProlif = true;

    double diameter = cellParams[6][0];

    mu = cellParams[0][0];
    kc = cellParams[1][0];
    damping = cellParams[2][0];
    maxOverlap = cellParams[3][0]*diameter;
    divProb = cellParams[4][0];
    deathProb = cellParams[5][0];
    radius = diameter/2.0;
    hypoxicL = cellParams[7][0];
    hypoxicP = cellParams[8][0];

    rmax = 1.5*diameter;

    targetLocation = {1e6,1e6};
}

void Cell::dieFromCD8(std::array<double, 3> otherX, double otherRadius, double kp, double dt) {
    /*
     * die from CTL based on a probability
     * contact required
     */
    if(type != 0){return;}

    if(calcDistance(otherX) <= radius+otherRadius){
        std::uniform_real_distribution<double> dis(0.0,1.0);
        if(dis(mt) < probTime(kp, dt)){
            state = -1;
        }
    }
}

void Cell::dieFromHypoxia(std::array<double, 3> tumorCenter, double tstep) {
    /*
      hypoxicL is lambda in an exp OR distance for a uniform prob
    */
    if(state != 0){return;}

    /*double prob = hypoxicP*exp(-hypoxicL*calcDistance(tumorCenter));
    std::uniform_real_distribution<double> dis(0.0,1.0);
    if(dis(mt) < probTime(prob, tstep)){
        state = -1;
    }*/
    
    std::uniform_real_distribution<double> dis(0.0,1.0);
    if(calcDistance(tumorCenter) < hypoxicL && dis(mt) < probTime(hypoxicP, tstep)){
	state = -1;
    }
}
