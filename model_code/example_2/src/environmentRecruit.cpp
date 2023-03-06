#include "Environment.h"

void Environment::recruitImmuneCells(double tstep) {
    // recruitment is scaled by number of cancer cells

    int numT8 = 0;
    int numC = 0;

    for(auto &c : cell_list){
        if(c.type == 1){
            numT8++;
        } else if(c.type == 0){
            numC++;
        }
    }

    double cd82c = static_cast<double>(numT8)/ static_cast<double>(numC);
    //double ratio = std::max(0.0, (1 - cd82c/cd8Ratio));
    cd82rec += tstep*cd8RecRate*static_cast<double>(numC)*static_cast<double>(cd82c < cd8Ratio);//*ratio;
    while (cd82rec >= 1) {
        std::array<double, 3> recLoc = recruitImmuneWhole();
        cell_list.emplace_back(recLoc, static_cast<int>(cell_list.size()), cellParams, "CD8", threeD, static_cast<double>(steps)*tstep/24);
        cd82rec -= 1;
    }
}

std::array<double, 3> Environment::recruitImmuneWhole() {
    /*
     * cells enter a distance, d, away from the tumor radius based on an exponential distribution
     * cells enter d away from a random edgeCell, such that the cell, edgeCell, and tumor center form a straight line
     */

    std::uniform_int_distribution<int> ec(0, edgeCells.size()-1);
    std::array<double, 3> x = edgeCells[ec(mt)];
    std::array<double, 3> dx = {x[0] - tumorCenter[0],
                                x[1] - tumorCenter[1],
                                x[2] - tumorCenter[2]};
    double norm = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
    double alpha = -log2(0.01);
    double lambda = alpha*0.693/recDist;
    std::exponential_distribution<double> loc(lambda);
    double distance = std::max(loc(mt), recDist);

    std::array<double, 3> recLoc = {distance*(dx[0]/norm) + x[0],
                                    distance*(dx[1]/norm) + x[1],
                                    distance*(dx[2]/norm) + x[2]};

    if(fabs(recLoc[0]) > 1e10 || fabs(recLoc[1]) > 1e10 || fabs(recLoc[2]) > 1e10){
        std::cout << "rl0: " << recLoc[0] << std::endl
        << "rl1: " << recLoc[1] << std::endl
        << "rl2: " << recLoc[2] << std::endl
        << "edge cell 0: " << x[0] << std::endl
                << "edge cell 1: " << x[1] << std::endl
                << "edge cell 2: " << x[2] << std::endl
                << "norm: " << norm << std::endl
                << "distance: " << distance << std::endl;

        throw std::runtime_error("recruiting cells");
    }

    return recLoc;
}