#include <iostream>
#include "Environment.h"

int main(int argc, char **argv) {
    std::string folder = argv[1];
    std::string paramSet = argv[2];
    std::string set = argv[3];

    std::string saveFld = "./"+folder+"/simulation_"+paramSet+"/set_"+set;
    std::string str = "mkdir -p "+saveFld;

    /*if(folder == "trainingData"){
        str = "mkdir -p ./trainingData/set_" + set;
        saveFld = "./trainingData/set_" + set;
    } else{
        str = "mkdir -p ./"+folder+"_"+paramSet+"/set_"+set;
        saveFld = "./"+folder+"_"+paramSet+"/set_"+set;
    }*/

    const char *command = str.c_str();
    std::system(command);

    str = "python3 genParams.py "+paramSet+" "+saveFld;
    command = str.c_str();
    std::system(command);

    double start = omp_get_wtime();
    Environment model(saveFld);
    model.simulate(0.25);
    double stop = omp_get_wtime();
    std::cout << "Duration: " << (stop-start)/(60*60) << std::endl;

    return 0;
}
