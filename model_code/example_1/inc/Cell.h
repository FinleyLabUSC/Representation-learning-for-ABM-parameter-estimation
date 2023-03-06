#ifndef IMMUNE_MODEL_CELL_H
#define IMMUNE_MODEL_CELL_H

#include <array>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <iostream>

class Cell{
public:
    /*
     * FUNCTIONS
     */

    // initialization
    Cell(std::array<double, 3> loc, int idx, std::vector<std::vector<double>> &cellParams, std::string cellType, double threeDimensional, double time);
    void initializeCancerCell(std::vector<std::vector<double>> &cellParams);
    void initializeCD8Cell(std::vector<std::vector<double>> &cellParams);

    // force functions
    std::array<double, 3> attractiveForce(std::array<double, 3> dx, double otherRadius);
    std::array<double, 3> repulsiveForce(std::array<double, 3> dx, double otherRadius);
    void calculateForces(std::array<double, 3> otherX, double otherRadius, int &otherType);
    void resolveForces(double dt);
    void resetForces();
    void neighboringCells(std::array<double, 3> otherX, int otherID);

    // overlap functions
    void calculateOverlap(std::array<double, 3> otherX, double otherRadius);
    void resetOverlap();
    void isCompressed();

    // cell behavior functions
    std::array<double, 4> proliferate(double dt);
    void age(double dt);
    void migrate(double dt, std::vector<std::array<double, 3>> edgeCells, std::array<double, 3> tumorCenter);
    void migrationTarget(std::array<double, 3> tumorCenter);

    // cell influences
    void addInfluence(std::array<double, 3> otherX, double otherInfluence, int otherType);
    void clearInfluence();

    // CD8 specific
    void pdl1Inhibition(std::array<double, 3> otherX, double otherRadius, double otherpdl1, double dt);

    // cancer specific
    void prolifState();
    void dieFromCD8(std::array<double, 3> otherX, double otherRadius, double kp, double dt);
    void inherit(double pd);
    void gainPDL1(double dt);

    // other functions
    double calcDistance(std::array<double, 3> otherX);
    void updateID(int idx);
    double calcInfDistance(double dist, double xth);
    static double calcNorm(std::array<double, 3> dx);
    static double probTime(double pInit, double dt);

    /*
     * PARAMETERS
     */

    // location
    std::array<double, 3> x;
    double threeD;
    std::array<double, 3> originalX;

    // physical properties
    double radius;
    bool compressed;
    double currentOverlap;
    std::vector<int> neighbors;

    // age, division, and lifespan
    double divProb;
    double deathProb;
    bool canProlif;

    // force properties
    double mu;
    double kc;
    double damping;
    double maxOverlap;
    double rmax;
    std::array<double, 3> currentForces;

    // migration
    double migrationSpeed;
    std::array<double, 3> targetLocation;
    double infiltrationDistance;
    double migrationBias;

    // cancer properties
    double pdl1Shift;

    // interactions with other cells
    double influenceRadius;
    double pdl1;
    double pdl1WhenExpressed;
    std::array<double, 3> influences;
    double probTh;

    // T cell killing
    double killProb;
    double influenceDec;

    // identification
    int id;
    int type;
    int state;
    double timeBorn;

private:
    std::mt19937 mt;
};

#endif //IMMUNE_MODEL_CELL_H
