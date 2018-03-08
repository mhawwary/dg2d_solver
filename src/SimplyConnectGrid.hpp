#ifndef SIMPLYCONNECTGRID_H
#define SIMPLYCONNECTGRID_H

#include"global_var.h"
#include"general_tools.h"
#include"../include/error.h"

class SimplyConnectGrid{

    public:

    SimplyConnectGrid(void);
    ~SimplyConnectGrid(void);

    void Setup(const int& oJm_, const int& oIm_
                   ,const double& xo_,const double& xf_
                   ,const double& yo_,const double& yf_
                   ,const char* odump_dir);
    void ConstructGrid();
    const char* getMeshfname(void){
        return meshfnamedir_.c_str();
    }

protected:
    void Initialize();
    void Reset();
    void CartesianGrid();
    void GridData();
    void WriteMesh(std::string& owrite_fname_);
    //void WriteMeshTecplot(std::string& owrite_fname_);

protected:
    char *dump_dir=nullptr;
    std::string meshfnamedir_;  // FullName with directory
    int Im,Jm,Nnodes,Ncells,Nfaces;
    double xo_=-1.0,xf_=1.0;
    double yo_=-1.0,yf_=1.0;
    /*
    -(Im,Jm): max counter/index in x and y dir respectively
    -Nn: no. of nodes
    -Nc: no. of cells
    -Nf: no. of faces
    -(x,y): grid points matrices,
    -nodes: global node ID matrix
    -(normx, normy):
    -faces: global IDs of faces
    -Rcell, Lcell: two arrays containing right and left cells for each face
    */
    double **x=nullptr,**y=nullptr;
    double **nodelist=nullptr;
    double *normx=nullptr,*normy=nullptr;
    int **facelist=nullptr;
    int *Rcell=nullptr,*Lcell=nullptr;

};

#endif
