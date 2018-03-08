#ifndef SIMDATA_H
#define SIMDATA_H

#include"../include/getpot.h"
#include"global_var.h"
#include"general_tools.h"


struct SimData {

    // Case parameters:
    //-----------------------------
    std::string case_title;    // parent folder for the case results
    std::string mesh_fname;   // mesh file name
    std::string case_postproc_dir;  // postprocessing directory
    std::string grid_type_;  // Cartesian
    int generate_grid_flag_=0;
    int Nnodes_x_ = 1;  // no. of nodes in x-direction
    int Nnodes_y_ = 1;  // no. of nodes in y-direction
    double x0_=0.0;  // x_begin of domain
    double xf_=1.0;  // x_end of domain
    double y0_=0.0;  // y_begin of domain
    double yf_=1.0;  // y_end of domain
    int uniform_=1;  // 0: for nonuniform mesh elements
    int refine_level_=0; // 0: no refinement
    //int Npplot = 1;  // to be depricated to use N_uniform_pts_per_elem_ instead
    //int N_exact_plot_pts=100;   // no. of points for exact solution plotting
    //int N_uniform_pts_per_elem_ = 100; // no. of global equally spaced points for plotting

    // Simulation parameters:
    //-----------------------------
    int unsteady_data_print_flag_ = 1;   // 0: use iter , 1: use time
    int unsteady_data_print_iter_ = 1000;  // iter no. for print of unsteady data
    double unsteady_data_print_time_ =1.0; // time point for print of unsteady data
    int restart_flag=0;      //  0: start a new simulation; 1: restart a previous simulation
    int restart_iter_=0;     // iteration to start from
    std::string Sim_mode;    // normal/test/error_analysis_dt/error_analysis_CFL/CFL_const/dt_const
    std::string case_no_;    // case no for turbulent flows

    // Wave parameters:
    //----------------------------
    double a_wave_x_=1.0;    // wave speed in x-dir
    double a_wave_y_=1.0;    // wave speed in y-dir
    int wave_form_ = 0;  // 0:KannanWang2009, 1:Trigonometric, 2:Gaussian wave
    // ./KannanWang2009 u(x,0) = A * sin ( f x + dphy ) * exp(C * y):
    double wave_freq_= 1.0;  // f
    double wave_amp_ = 1.0;  // A
    double wave_shift = 0.0; // phy
    double exp_exponent_const = 1.0; // C
    /* ../Trigonometric , sine or cosine
     * u(x,y,t) = A .* sin ( fx  x + dphy_x ) .* sin ( fy  y + dphy_y )
     * .*exp(-(fx^2+fy^2) * t) + C
     * */
    std::string x_wave_type_;  // sin/cos
    std::string y_wave_type_;  // sin/cos
    double x_wave_freq_= 2.0;  // fx/pi
    double y_wave_freq_= 2.0;  // fy/pi
    //double wave_amp_ = 1.0;  // A
    double wave_const = 0.0; // C
    double x_wave_phase = 0.0; // dphy_x
    double y_wave_phase = 0.0; // dphy_y
    // ../Gaussian , u(x) = A * exp( M *x^2) :
    double Gaussian_amp_ = 1.0;   // A
    double Gaussian_exponent_ = -40;  // M

    // Space Solver parameters:
    //-----------------------------
    std::string eqn_set;    //  Advection / Diffusion / Advection_Diffusion
    std::string eqn_type_;  // inv_burger / linear_advec / linear_diffus
    int poly_order_=0;      // polynomial order
    std::string schemes_data_dir_;    // parent folder for schemes data
    // ./advec_eqn:
    double upwind_param_=1.0;     // Beta for upwind-central
    // ./heat_eqn:
    double penalty_param_=1.0;    // Epsilon for relaxation
    double thermal_diffus=1.0;    // kappa
    std::string diffus_scheme_type_;   // BR2 / SIP / .....

    // Time Solver parameters:
    //--------------------------------
    int calc_dt_flag=1; // 1: specify CFL and calc dt, 0: specify dt and calc CFL
    int calc_dt_adv_diffus_flag=0;
    /* 0: based on advection,
     * 1: based on diffusion,
     * 2: based on combined advection-diffusion */
    double dt_ = 1e-3;  // dt time step
    double t_init_ = 0.0;  // initial time
    double t_end_ =1e20;  // end time
    double maxIter_ = 1e10; // maximum number of iterations
    double CFL_    = 1.0;   // CFL no.
    double Nperiods = 1.0; // no. of periods for simulation
    int end_of_sim_flag_=0;  // 1: use max_iteration as a stopping criteria if not converged or diverged
    // ./explicit:
    int RK_order_=0;        // Runge-Kutta type (0: euler FT, 2: SSPRK22, 3: SSPRK33)
    //------------------------------------------------------------------

    // Function members:
    //--------------------------
    void Parse(const std::string &fname);

};

#endif
