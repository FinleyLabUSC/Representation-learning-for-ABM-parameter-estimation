#include "Environment.h"

void Environment::printStep(double time) {
    int numT8 = 0;
    int numT8s = 0;
    int numC = 0;

    for(auto &cell : cell_list){
        if(cell.type == 1){
            if(cell.state == 1){numT8++;}
            if(cell.state == 2){numT8s++;}
        } else if(cell.type == 0){
            numC++;
        }
    }
    std::cout << "************************************\n"
              << "Time (d): " << time/24 << std::endl
              << "Cancer: " << numC << std::endl
              << "CD8: " << numT8 << " " << numT8s << std::endl;
}