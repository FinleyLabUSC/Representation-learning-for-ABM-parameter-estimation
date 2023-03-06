#ifndef IMMUNE_MODEL_ENVIRONMENT_H
#define IMMUNE_MODEL_ENVIRONMENT_H

#include <vector>
#include <algorithm>
#include <random>
#include "Cell.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <omp.h>

class Environment{
public:
    Environment(std::string saveFld);
    void simulate(double tstep);

private:
    void runCells(double tstep);
    void neighborInfluenceInteractions(double tstep);
    void internalCellFunctions(double tstep);
    void recruitImmuneCells(double tstep);
    std::array<double, 3> recruitImmuneWhole();

    void save(double tstep);
    void loadParams();

    void calculateForces(double tstep);

    void printStep(double time);
    void tumorSize();

    double probTime(double pInit, double tstep);

    double dt;
    double threeD;

    // cell lists
    std::vector<Cell> cell_list;
    std::vector<std::array<double, 3>> edgeCells;

    // parameter lists
    std::vector<std::vector<double>> cellParams;
    std::vector<double> recParams;
    std::vector<double> envParams;

    std::string saveDir;
    int steps;
    double cd8RecRate;
    double cd8Ratio;
    double recDist;
    double cd82rec;
    double tumorRadius;
    std::array<double, 3> tumorCenter;

    // environment params
    double simulationDuration;

    std::mt19937 mt;
};

#endif //IMMUNE_MODEL_ENVIRONMENT_H
