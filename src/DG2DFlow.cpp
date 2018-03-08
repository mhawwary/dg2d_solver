
#include"SimCase.hpp"
#include"general_tools.h"
#include"solver_tools.h"
#include"../include/error.h"

using namespace std;

int main(int argc, char** argv){


    if (argc < 2) {
        FatalError_exit("ERROR: No input file specified ... ");
        return(0);
    }

    std::string input_fname;
    input_fname = argv[argc-1];  // input file name with directory

    SimCase Case;
    Case.setup(input_fname);  // Setup and parse input parameters
    Case.InitSim();           // Preprocessing steps
    Case.RunSim();            // Main Solution Loop
    //Case.PostProcess(iter);       // Dumping Simulation Post Processing data

    return 0;
}


