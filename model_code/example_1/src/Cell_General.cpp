#include "Cell.h"

/*
 * BIOLOGICAL AND MODELING INFO
 * ----------------------------
 * - cell forces taken from Osborne 2017
 */

/*
 * CELL TYPES
 * 0 - cancer
 * 1 - CD8
 *
 * CELL STATES
 * -1 - dead
 * 0 - alive (cancer)
 * 1 - active (CD8)
 * 2 - suppressed (CD8)
 *
 */

// ********************
// INITIALIZE CELL TYPE
Cell::Cell(std::array<double, 3> loc, int idx, std::vector<std::vector<double>> &cellParams, std::string cellType, double threeDimensional, double time): mt((std::random_device())()) {
    if(cellType == "cancer"){
        type = 0;
    } else if(cellType == "CD8"){
        type = 1;
    }

    x = loc;
    originalX = loc;
    id = idx;
    threeD = threeDimensional;
    timeBorn = time;

    /*
     * initialize parameters as 0
     * then initialize only the relevant parameters
     */
    radius = 0;
    compressed = false;
    currentOverlap = 0;
    divProb = 0;
    deathProb = 0;
    canProlif = false;
    mu = 0;
    kc = 0;
    damping = 0;
    maxOverlap = 0;
    rmax = 0;
    currentForces = {0,0,0};
    migrationSpeed = 0;
    targetLocation = {0,0,0};
    infiltrationDistance = 0;
    migrationBias = 0;
    influenceRadius = 0;
    pdl1 = 0;
    pdl1Shift = 0;
    pdl1WhenExpressed = 0;
    killProb = 0;
    for(int i=0; i<influences.size(); ++i){
        influences[i] = 0;
    }
    influenceDec = 0;
    state = 0;

    // for influence distance, assume a soft-cutoff where p(distance) = probTh
    probTh = 0.001;

    if(cellType == "cancer"){
        initializeCancerCell(cellParams);
    } else if(cellType == "CD8"){
        initializeCD8Cell(cellParams);
    } else{
        std::cout << "Cell type requested: " << cellType << std::endl;
        throw std::runtime_error("Cell::Cell -> unavailable cell type");
    }
}
// ********************

// ************************
// BIO-MECHANICAL FUNCTIONS
// ---------------
// FORCE FUNCTIONS
std::array<double, 3> Cell::attractiveForce(std::array<double, 3> dx, double otherRadius) {
    double dxNorm = calcNorm(dx);
    std::array<double, 3> dxUnit = {dx[0]/dxNorm, dx[1]/dxNorm,  dx[2]/dxNorm};
    double sij = radius + otherRadius;

    double scaleFactor = mu*(dxNorm - sij)*exp(-kc*(dxNorm - sij)/sij);
    double F0 = dxUnit[0]*scaleFactor;
    double F1 = dxUnit[1]*scaleFactor;
    double F2 = dxUnit[2]*scaleFactor;

    return {F0, F1, F2};
}

std::array<double, 3> Cell::repulsiveForce(std::array<double, 3> dx, double otherRadius) {
    double dxNorm = calcNorm(dx);
    std::array<double, 3> dxUnit = {dx[0]/dxNorm, dx[1]/dxNorm, dx[2]/dxNorm};
    double sij = radius + otherRadius;

    double scaleFactor = mu*sij*log10(1 + (dxNorm - sij)/sij);
    double F0 = dxUnit[0]*scaleFactor;
    double F1 = dxUnit[1]*scaleFactor;
    double F2 = dxUnit[2]*scaleFactor;

    return {F0, F1, F2};
}

void Cell::calculateForces(std::array<double, 3> otherX, double otherRadius, int &otherType) {
    double distance = calcDistance(otherX);
    if(distance < rmax){
        std::array<double, 3> dx = {(otherX[0]-x[0]),
                                    (otherX[1]-x[1]),
                                    (otherX[2]-x[2])};
        if(distance < (radius + otherRadius)){
            std::array<double, 3> force = repulsiveForce(dx, otherRadius);
            currentForces[0] += force[0];
            currentForces[1] += force[1];
            currentForces[2] += force[2];
        } else if(type == 0 && otherType == 0){ // attraction if both are cancer cells
            std::array<double, 3> force = attractiveForce(dx, otherRadius);
            currentForces[0] += force[0];
            currentForces[1] += force[1];
            currentForces[2] += force[2];
        }
    }
}

void Cell::resolveForces(double dt) {
    x[0] += (dt/damping)*currentForces[0];
    x[1] += (dt/damping)*currentForces[1];
    x[2] += (dt/damping)*currentForces[2];

    resetForces();
}

void Cell::resetForces() {
    /*
     * resets forces with a slight randomizing factor
     */
    double D = 1;
    std::uniform_real_distribution<double> dis(-D, D);
    currentForces = {dis(mt),dis(mt),dis(mt)*threeD};
}

void Cell::neighboringCells(std::array<double, 3> otherX, int otherID){
    /*
     * determine which cells are within 2*maximum interaction distance
     */
    double dis = calcDistance(otherX);
    if(dis <= 10*rmax){
        neighbors.push_back(otherID);
    }
}

// OVERLAP FUNCTIONS
void Cell::calculateOverlap(std::array<double, 3> otherX, double otherRadius) {
    double distance = calcDistance(otherX);
    if(distance < radius + otherRadius){
        currentOverlap += radius + otherRadius - distance;
    }
}

void Cell::resetOverlap() {
    currentOverlap = 0;
}

void Cell::isCompressed() {
    compressed = currentOverlap > maxOverlap;
    resetOverlap();
}
// ************************

// *****************************
// SPECIFIC BIOLOGICAL FUNCTIONS
// -----------------------
// CELL BEHAVIOR FUNCTIONS
std::array<double, 4> Cell::proliferate(double dt) {
    // positions 0, 1, and 2 are cell location
    // position 3 is boolean didProliferate?
    if(!canProlif){return {0,0, 0,0};}

    std::uniform_real_distribution<double> dis(0.0, 1.0);
    if(dis(mt) < probTime(divProb, dt)){
        // place daughter cell a random angle away from the mother cell
        std::uniform_real_distribution<double> rd(-1.0,1.0);
        std::array<double, 3> dx = {rd(mt),
                                      rd(mt),
                                      rd(mt)*threeD};
        double norm = calcNorm(dx);
        return{radius*(dx[0]/norm)+x[0],
               radius*(dx[1]/norm)+x[1],
               radius*(dx[2]/norm)+x[2],
               1};
    } else{
        return {0,0, 0,0};
    }
}

void Cell::age(double dt) {
    /*
     * cells die based on a probability equal to 1/lifespan
     */
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    if(dis(mt) < probTime(deathProb, dt)){
        state = -1;
    }
}

void Cell::migrate(double dt, std::vector<std::array<double, 3>> edgeCells, std::array<double, 3> tumorCenter) {
    /*
     * biased random-walk towards their target
     *
     * bias is done by generating a random vector, then adding it to the correct direction with scaling
     */
    std::uniform_real_distribution<double> dis(-1,1);

    std::array<double, 3> dx = {targetLocation[0] - x[0],
                                targetLocation[1] - x[1],
                                targetLocation[2] - x[2]};

    std::array<double, 3> randomVector = {dis(mt),
                                          dis(mt),
                                          dis(mt)};

    double norm = calcNorm(dx);
    double normRV = calcNorm(randomVector);
    for(int i=0; i<dx.size(); ++i){
        dx[i] /= norm;
        randomVector[i] /= normRV;
        dx[i] = migrationBias*dx[i] + (1 - migrationBias)*randomVector[i];
    }
    dx[2] *= threeD;
    norm = calcNorm(dx);

    // determine distance from tumor edge
    /*
     * first, find nearest edge cell
     */
    int idx = 0;
    double distFromEdge = 1e6;
    for(int i=0; i<edgeCells.size(); ++i){
        double dist = calcDistance(edgeCells[i]);
        if(dist < distFromEdge){
            distFromEdge = dist;
            idx = i;
        }
    }
    double edgeDistFromCenter = calcNorm({edgeCells[idx][0] - tumorCenter[0],
                                          edgeCells[idx][1] - tumorCenter[1],
                                          edgeCells[idx][2] - tumorCenter[2]});

    // if migrating into the tumor, stop if the infiltration distance is reached
    double distanceFromCenter = calcDistance(tumorCenter);
    for(int i=0; i<x.size(); ++i){
        x[i] += dt*migrationSpeed*(dx[i]/norm)*(edgeDistFromCenter - distanceFromCenter < infiltrationDistance*edgeDistFromCenter);
    }
    if(fabs(x[0]) > 1e10 || fabs(x[1]) > 1e10){
        std::cout << "Error\n";
        std::cout << "Type: " << type << std::endl;
        std::cout << "X: " << x[0] << std::endl;
        std::cout << "original X: " << originalX[0] << " " << originalX[1] << std::endl;
        std::cout << "Mig Speed: " << migrationSpeed << std::endl;
        std::cout << "dx/norm: " << dx[0]/norm << std::endl;
        std::cout << "dist from center: " << distanceFromCenter << std::endl;
        throw std::runtime_error("migration");
    }
}

void Cell::migrationTarget(std::array<double, 3> tumorCenter) {
    if(type != 0) {
        // immune cells migrate towards the tumor
        targetLocation = tumorCenter;
    }
}

void Cell::prolifState() {
    /*
     * cancer cells and CD8 can proliferate
     * right now, CD8 proliferation prob is set to 0, however leaving it in for future changes
     */
    if(type == 0){
        canProlif = !(state == -1 || compressed);
    }
}

// CELL INFLUENCE
void Cell::addInfluence(std::array<double, 3> otherX, double otherInfluence, int otherState) {
    /*
     * determine influence based on distance for each cell state
     *
     * I believe I'm handling the probabilities correctly
     * totalProb = 1 - (1-p1)*(1-p2)*...*(1-pn)
     * the commented out way of just summing the probabilities is probably incorrect
     */
    if(otherState == -1){return;}

    //influences[otherState] *= calcInfDistance(calcDistance(otherX), otherInfluence);
    influences[otherState] = 1 - (1 - influences[otherState])*(1 - calcInfDistance(calcDistance(otherX), otherInfluence));
}

void Cell::clearInfluence() {
    for(int i=0; i<influences.size(); ++i){
        influences[i] = 0;
    }
}
// *****************************

// ***************
// OTHER FUNCTIONS
// ----------------------
// MATHEMATICAL FUNCTIONS
double Cell::calcDistance(std::array<double, 3> otherX) {
    double d0 = (otherX[0] - x[0]);
    double d1 = (otherX[1] - x[1]);
    double d2 = (otherX[2] - x[2]);

    return sqrt(d0*d0 + d1*d1 + d2*d2);
}

double Cell::calcInfDistance(double dist, double xth) {
    /*
     * calculate influence using an exponential decay based on distance from cell center
     */
    double alpha = -log2(probTh);
    double lambda = alpha*0.693/xth;

    return exp(-lambda*dist);
}

double Cell::calcNorm(std::array<double, 3> dx){
    return sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
}

double Cell::probTime(double pInit, double tstep) {
    /*
     * assumes an initial probability at a 1 hr timestep
     *
     * p0 = initial probability
     * dt0 = initial timestep
     * dtn = new timestep
     * steps = dtn/dt0
     * p = 1 - (1 - p0)^steps
     */

    return 1 - pow((1 - pInit), tstep);
}

// BOOK-KEEPING
void Cell::updateID(int idx) {
    id = idx;
}
// ***************