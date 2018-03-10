#ifndef SIMCASE_H
#define SIMCASE_H

#include"Mesh.hpp"
#include"SimplyConnectGrid.hpp"
#include"general_tools.h"
#include"solver_tools.h"
#include"SimData.hpp"
#include"SpaceSolver.hpp"
#include"dg_2d_diffus_solver.hpp"
#include"ExplicitTimeSolver.hpp"

using namespace std;

class SimCase {

public:
//  Construction Functions :
  SimCase(void);
  ~SimCase(void);
  void setup(const std::string &input_fname_);
  void setup_output_directories();
  void InitSim();
  void RunSim();
  void PostProcess(const int& iter_);

protected:
  void logo();
  void copy_problem_inputdata();
  void copyFile(const std::string& fileNameFrom
                , const std::string& fileNameTo);


protected:
  //std::string input_fname;  // input file name
  SimData simdata ;
  Mesh     *grid=nullptr;
  MeshData *grid_data=nullptr;
  SpaceSolver *dg_solver_=nullptr;
  ExplicitTimeSolver *time_solver_=nullptr;

};

#endif
