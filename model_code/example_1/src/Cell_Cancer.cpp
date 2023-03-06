#include "Cell.h"

void Cell::initializeCancerCell(std::vector<std::vector<double>> &cellParams) {
    state = 0;
    canProlif = true;

    double diameter = cellParams[8][0];

    mu = cellParams[0][0];
    kc = cellParams[1][0];
    damping = cellParams[2][0];
    maxOverlap = cellParams[3][0]*diameter;
    divProb = cellParams[4][0];
    deathProb = cellParams[5][0];
    pdl1WhenExpressed = cellParams[6][0];
    pdl1Shift = cellParams[7][0];
    radius = diameter/2.0;

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

void Cell::inherit(double pd) {
    /*
     * daughter cells have the same EMT and PD-L1 as the mother cell
     */
    if(type != 0){return;}

    pdl1 = pd;
}

void Cell::gainPDL1(double dt) {
    if(type != 0){return;}

    // induced by ifn-y secreting cells
    // posInfluence is Th + active CD8
    double posInfluence = influences[1];
    std::uniform_real_distribution<double> dis(0.0,1.0);
    if(dis(mt) < probTime(posInfluence, dt)){
        pdl1 += pdl1Shift;
    }
    pdl1 = std::min(pdl1,pdl1WhenExpressed);
}