#include "Environment.h"

void Environment::neighborInfluenceInteractions(double tstep) {

    /*
     * FIRST LOOP
     * - determine neighbors
     * - determine influences on a cell
     *
     * SECOND LOOP
     * - determine cancer cell neighbor distances
     *
     * THIRD LOOP
     * - differentiate macrophages and CD4
     * - set CD8 killProb
     * - inhibit CD8
     *
     * FOURTH LOOP
     * - cancer cell migratory
     * - gain PD-L1
     * - CD8 kill cancer cell
     */

#pragma omp parallel for
    for(int i=0; i<cell_list.size(); ++i){
        cell_list[i].neighbors.clear();
        cell_list[i].clearInfluence();
        for(auto &c : cell_list){
            // assume that a cell cannot influence itself
            if(cell_list[i].id != c.id){
                cell_list[i].neighboringCells(c.x, c.id);
                cell_list[i].addInfluence(c.x, c.influenceRadius, c.state);
            }
        }
    }

#pragma omp parallel for
    for(int i=0; i<cell_list.size(); ++i){
        if(cell_list[i].type == 1 && cell_list[i].state == 1){
            for(auto &c : cell_list[i].neighbors){
                if(cell_list[c].type == 0){
                    cell_list[i].pdl1Inhibition(cell_list[c].x, cell_list[c].radius, cell_list[c].pdl1, tstep);
                }
            }
        }
    }

#pragma omp parallel for
    for(int i=0; i<cell_list.size(); ++i){
        if(cell_list[i].type == 0){
            cell_list[i].gainPDL1(tstep);
            // die from neighboring CD8
            for(auto &c : cell_list[i].neighbors){
                if(cell_list[c].type == 1 && cell_list[c].state == 1){
                    cell_list[i].dieFromCD8(cell_list[c].x, cell_list[c].radius, cell_list[c].killProb, tstep);
                }
            }
        }
    }
}

void Environment::calculateForces(double tstep) {
    /*
     * 1. Calculate total force vector for each cell
     * 2. Resolve forces on each cell
     * 3. Determine current overlap for each cell
     * 4. Determine if each cell is compressed
     */

    // divide tstep into smaller steps for solving
    // only solve forces between neighboring cells to improve computation time
    int Nsteps = static_cast<int>(tstep/dt);

    // determine migration target
#pragma omp parallel for
    for(int i=0; i<cell_list.size(); ++i){
        cell_list[i].migrationTarget(tumorCenter);
    }

    // iterate thru Nsteps, calculating and resolving forces between neighbors
    // also includes migration
    for(int q=0; q<Nsteps; ++q){
        // migrate first
#pragma omp parallel for
        for(int i=0; i<cell_list.size(); ++i){
            cell_list[i].migrate(dt, edgeCells, tumorCenter);
        }

        // calc forces
#pragma omp parallel for
        for(int i=0; i<cell_list.size(); ++i){
            for(auto &c : cell_list[i].neighbors){
                cell_list[i].calculateForces(cell_list[c].x, cell_list[c].radius, cell_list[c].type);
            }
        }

        // resolve forces
#pragma omp parallel for
        for(int i=0; i<cell_list.size(); ++i){
            cell_list[i].resolveForces(dt);
        }
    }

    // calculate overlap for cancer cells and CD8
#pragma omp parallel for
    for(int i=0; i<cell_list.size(); ++i){
        if(cell_list[i].type == 0 || cell_list[i].type == 3){
            for(auto &c : cell_list[i].neighbors){
                if(cell_list[c].type == 0){
                    cell_list[i].calculateOverlap(cell_list[c].x, cell_list[c].radius);
                }
                if(cell_list[c].type == 3 && cell_list[i].type == 3){
                    cell_list[i].calculateOverlap(cell_list[c].x, cell_list[c].radius);
                }
            }
            cell_list[i].isCompressed();
            cell_list[i].prolifState();
        }
    }
}

void Environment::internalCellFunctions(double tstep) {
    /*
     * cell death via aging
     * cell proliferation
     * remove cell if out of bounds
     */
    int numCells = cell_list.size();
    for(int i=0; i<numCells; ++i){
        cell_list[i].age(tstep);
        if(cell_list[i].type == 0){
            std::array<double, 4> newLoc = cell_list[i].proliferate(tstep);
            if(newLoc[3] == 1){
                cell_list.push_back(Cell({newLoc[0], newLoc[1], newLoc[2]}, cell_list.size(), cellParams, "cancer", threeD, static_cast<double>(steps)*tstep/24));
                cell_list[cell_list.size() - 1].inherit(cell_list[i].pdl1);
            }
        }
        if(cell_list[i].type == 1){
            std::array<double, 4> newLoc = cell_list[i].proliferate(tstep);
            if(newLoc[3] == 1){
                cell_list.push_back(Cell({newLoc[0], newLoc[1], newLoc[2]}, cell_list.size(), cellParams, "CD8", threeD, static_cast<double>(steps)*tstep/24));
            }
        }
    }

    // remove dead cells
    std::vector<Cell> new_cell_list;
    for(auto & cell : cell_list){
        if(cell.state != -1){
            new_cell_list.push_back(cell);
        }
    }
    cell_list = new_cell_list;

    // shuffle cell list
    std::shuffle(std::begin(cell_list), std::end(cell_list), mt);
    for(int i=0; i<cell_list.size(); ++i){
        cell_list[i].updateID(i);
    }
}

void Environment::runCells(double tstep) {
    neighborInfluenceInteractions(tstep);
    calculateForces(tstep);
    internalCellFunctions(tstep);
}