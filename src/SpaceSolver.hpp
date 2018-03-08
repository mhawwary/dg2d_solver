#ifndef SPACESOLVER_H
#define SPACESOLVER_H

#include"MeshData.h"
#include"SimData.hpp"
#include"quadrature.h"
#include"quadrature_2d.h"
#include"general_tools.h"
#include"global_var.h"
#include"mysparse_matrix.h"

class SpaceSolver{

public:

    double T_period=1;

    //  Construction Functions :
    SpaceSolver(){}
    virtual ~SpaceSolver(){}

    virtual void setup_solver(MeshData*& meshdata_, SimData& simdata_)=0;
    virtual void InitSol()=0;
    virtual void UpdateResid(double *o_resid_, double *o_qn_)=0;
    virtual void Compute_vertex_sol()=0;

    // Input/Output members
    void UpdatePhyTime(const double& dt_){

        phy_time += dt_;

        return;
    }

    void SetPhyTime(const double &time_){

        phy_time=time_;

        return;
    }

    double GetTimeStep(){

        return time_step;
    }

    double GetLastTimeStep(){

        return last_time_step;
    }

    double GetCFL(){

        return CFL;
    }

    int GetNdof(){
        return Ndof;
    }

    int GetNelem(){
        return grid_->Nelem;
    }

    double* GetNumSol(){

        return Qn;
    }

    double* GetVertexNumSol(){

        return Qv;
    }

    double* GetVertexExactSol(){

        return Qv_exact;
    }

    double GetPhyTime(){

        return phy_time;
    }

protected:
    //Element mapping functions:
    virtual void eval_local_xy_coord(const int& oelemID
                                     ,const double& xi_,const double& eta_
                                     ,double& o_x_,double& o_y_)=0;
    virtual double eval_bilinear_mapshapefunc(const double& xi_
                                              ,const double& eta_
                                              ,const int& index_)=0;

    virtual double eval_init_sol(const double& xx_,const double& yy_)=0;
    virtual double eval_1dbasis_poly(const double& xi_, const int& basis_k_)=0;
    virtual double eval_2dbasis_poly(const double& xi_,const double& eta_
                                     , const int& basis_k_)=0;
    virtual double initSol_legendre_proj(const int& eID, const int &basis_id,
                                         const GaussQuad2D & quad_)=0;
    virtual double exactSol_legendre_proj(const int& eID, const int &basis_id,
                                         const GaussQuad2D & quad_)=0;
    virtual double evalSolution(const double* q_
                                , const double& xi_pt_
                                , const double& eta_pt_)=0;

    virtual void CalcTimeStep()=0;
    virtual void ComputeExactSolShift()=0;
    virtual void Compute_projected_exact_sol()=0;
    virtual void Compute_exact_vertex_sol()=0;
    virtual void UpdateResidOneCell(const int& eID_, double* resid_
                                    ,double *o_qn_)=0;

    /* These pointers are passed to the space solver
   *  and are not supposed to be freed in this scope
   */
    MeshData *grid_=nullptr;
    SimData  *simdata_=nullptr;

    //Other variables
    double phy_time=0.0;
    double time_step=1e-5;
    double last_time_step=1e-5;
    double CFL=1.0;
    int Porder=1;
    int Ndof = 1;
    int Ndof_1D=1;

    // Locally defined arrays and can be freed in this scope
    double *Qn=nullptr;      // Nelem * Ndof long
    double *Qex_proj=nullptr; // projected exact solution , Nelem * Ndof long
    double *Qv_exact=nullptr;  // exact, Nnodes long
    double *Qv=nullptr;       // numerical, Nnodes long

    double max_eigen_advec=0.0; // maximum eigenvalue for adevction
    double x_exact_sol_shift=0.;
    double y_exact_sol_shift=0.;
    double x_wave_length_=0.;
    double y_wave_length_=0.;

    int Nquad_1d=1; // 1D Gauss Quadrature rules
    int Nquad_2d=1; // 2D Gauss Quadrature rules
    GaussQuad quad_1d_;
    GaussQuad2D quad_2d_;

    double *Lk_norm2sq=nullptr; // squared norm of Lk's := ||Lk||^{2}_{2}

};

#endif
