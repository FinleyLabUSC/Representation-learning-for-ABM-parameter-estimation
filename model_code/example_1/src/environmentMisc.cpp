#include "Environment.h"

void Environment::tumorSize() {
   double avgX = 0;
   double avgY = 0;
   double avgZ = 0;
   int numC = 0;
   for(auto &c : cell_list){
       if(c.type == 0) {
           avgX += c.x[0];
           avgY += c.x[1];
           avgZ += c.x[2];
           numC++;
       }
   }

   avgX /= static_cast<double>(numC);
   avgY /= static_cast<double>(numC);
   avgZ /= static_cast<double>(numC);

   tumorCenter = {avgX, avgY, avgZ};

   double dist = 0;
   for(auto & cell : cell_list){
       if(cell.type == 0){
           dist = std::max(dist, cell.calcDistance(tumorCenter));
       }
   }

   tumorRadius = dist;

   edgeCells.clear();
   for(int i=0; i<cell_list.size(); ++i){
       if(cell_list[i].type != 0){continue;}

       double radius = cell_list[i].radius;
       std::array<double, 3> x = cell_list[i].x;
       std::array<double, 3> dx = {x[0] - tumorCenter[0],
                                   x[1] - tumorCenter[1],
                                   x[2] - tumorCenter[2]};
       double norm = cell_list[i].calcNorm(dx);
       if(norm < 0.75*tumorRadius){continue;}
       dx[0] /= norm;
       dx[1] /= norm;
       dx[2] /= norm;
       std::array<double, 3> nx = {x[0] + 4*radius*dx[0],
                                   x[1] + 4*radius*dx[1],
                                   x[2] + 4*radius*dx[2]};
       bool free = true;
       for(int j=0; j<cell_list.size(); ++j){
           if(i == j){continue;}
           if(cell_list[j].type != 0){continue;}
           if(cell_list[j].calcDistance(nx) < 2*cell_list[j].radius){
               free = false;
               break;
           }
       }

       if(free){
           edgeCells.push_back(x);
       }
   }
}

double Environment::probTime(double pInit, double tstep) {
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