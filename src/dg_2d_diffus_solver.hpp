#ifndef DG_2D_DIFFUS_SOLVER_H
#define DG_2D_DIFFUS_SOLVER_H


#include"SpaceSolver.hpp"


class DG_2D_DIFFUS_SOLVER:public SpaceSolver{

public:

    //  Construction Functions:
    DG_2D_DIFFUS_SOLVER(void){}
    virtual ~DG_2D_DIFFUS_SOLVER(void);

    virtual void setup_solver(MeshData*& meshdata_, SimData& simdata_);
    virtual void InitSol();

    //Element mapping functions:
    virtual void eval_local_xy_coord(const int& oelemID
                                     ,const double& xi_,const double& eta_
                                     ,double& o_x_,double& o_y_);
    virtual double eval_bilinear_mapshapefunc(const double& xi_
                                              ,const double& eta_
                                              ,const int& index_);
    virtual void UpdateResid(double *o_resid_, double *o_qn_);
    virtual void UpdateResidOneCell(const int& eID_, double* resid_
                                    ,double *o_qn_);

    virtual void Compute_vertex_sol();
    //   virtual void Compute_cont_sol(){}
    //   virtual double ComputePolyError();
    //   virtual double L1_error_projected_sol();
    //   virtual double L2_error_projected_sol();
    //   virtual double L1_error_average_sol();
    //   virtual double L2_error_average_sol();
    //   virtual double L1_error_nodal_gausspts();
    //   virtual double L2_error_nodal_gausspts();
    //   virtual double L1_error_nodal_gausspts_proj();
    //   virtual double L2_error_nodal_gausspts_proj();
    //   virtual void print_cont_vertex_sol();
    //   virtual void print_average_sol();
    //   virtual void dump_errors(double& L1_proj_sol_,double& L2_proj_sol_
    //                    ,double& L1_aver_sol_,double& L2_aver_sol_
    //                    ,double& L1_nodal_gausspts, double& L2_nodal_gausspts);
    //   virtual void dump_discont_sol();
    //   virtual void dump_timeaccurate_sol();

    virtual double eval_init_sol(const double& xx_,const double& yy_);
    virtual double eval_1dbasis_poly(const double& xi_, const int& basis_k_);
    virtual double eval_2dbasis_poly(const double& xi_,const double& eta_
                                     , const int& basis_k_);

    virtual double initSol_legendre_proj(const int& eID, const int &basis_id,
                                         const GaussQuad2D & quad2d_);
    virtual double exactSol_legendre_proj(const int& eID, const int &basis_id,
                                          const GaussQuad2D & quad2d_);
    virtual double evalSolution( const double* q_
                                 , const double& xi_pt_
                                 , const double& eta_pt_);

    virtual void CalcTimeStep();
    virtual void ComputeExactSolShift();
    virtual void Compute_exact_vertex_sol();
    virtual void Compute_projected_exact_sol();

    //Debug functions:
    void test_local_xy_coord();
    void test_eval_2dbasis();
    void Compute_exact_vertex_sol_usingproj();
    void TestGhostElements();

private:
    void Reset_solver();
    void read_sparsematrix(MYSPARSE_MATRIX& smatrix_, std::string& readfname);
    void load_scheme_matrices(void);
    double eval_exact_trigsol_diffus(const double& xx_,const double& yy_
                                     , const double& tt_);

private:
    double eta_face=1.0;  // penalty parameter;
    // Loaded Matrices for solution update, useful only for cartesian mesh
    MYSPARSE_MATRIX K0,K1,K2,K3;
    MYSPARSE_MATRIX Q0,Q1,Q2,Q3;
    MYSPARSE_MATRIX L,Q;
};

#endif
