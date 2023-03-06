#include "Environment.h"

void Environment::loadParams() {
    std::ifstream dataCP(saveDir+"/params/cellParams.csv");
    std::string line;
    while(std::getline(dataCP, line)){
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<double> parsedRow;
        while(std::getline(lineStream, cell, ',')){
            parsedRow.push_back(std::stod(cell));
        }
        cellParams.push_back(parsedRow);
    }
    dataCP.close();

    std::ifstream dataRP(saveDir+"/params/recParams.csv");
    while(std::getline(dataRP, line)){
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<double> parsedRow;
        while(std::getline(lineStream, cell, ',')){
            parsedRow.push_back(std::stod(cell));
        }
        recParams.push_back(parsedRow[0]);
    }
    dataRP.close();

    std::ifstream dataEP(saveDir+"/params/envParams.csv");
    while(std::getline(dataEP, line)){
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<double> parsedRow;
        while(std::getline(lineStream, cell, ',')){
            parsedRow.push_back(std::stod(cell));
        }
        envParams.push_back(parsedRow[0]);
    }
    dataEP.close();
}

void Environment::save(double tstep) {

    std::ofstream myfile;

    int numCancer = 0;
    for(auto & cell : cell_list){
        if(cell.type == 0 && cell.state != -1){
            numCancer++;
        }
    }

    int c8 = 0;
    for(auto &c : cell_list){
        if(c.type == 1){
            c8++;
        }
    }

    double time = steps*tstep/24;
    myfile.open(saveDir+"/outputs.csv");
    myfile << time << "," << numCancer << "," << c8
           << "," << tumorCenter[0] << "," << tumorCenter[1] << "," << tumorCenter[2] << "," << tumorRadius << std::endl;
    myfile.close();

    myfile.open(saveDir+"/cancerAlive.csv");
    for(auto &cell : cell_list){
        if(cell.type == 0 && cell.state == 0) {
            myfile << cell.x[0] << "," << cell.x[1] << "," << cell.x[2] << "," << cell.radius << std::endl;
        }
    }
    myfile.close();

    myfile.open(saveDir+"/cancerDead.csv");
    for(auto &cell : cell_list){
        if(cell.type == 0 && cell.state == -1) {
            myfile << cell.x[0] << "," << cell.x[1] << "," << cell.x[2] << "," << cell.radius << std::endl;
        }
    }
    myfile.close();

    myfile.open(saveDir+"/cd8Cells.csv");
    for(auto &cell : cell_list){
        if(cell.type == 1) {
            myfile << cell.x[0] << "," << cell.x[1] << "," << cell.x[2] << "," << cell.radius << "," << cell.infiltrationDistance << std::endl;
        }
    }
    myfile.close();
}
