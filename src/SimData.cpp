#include"SimData.hpp"

void SimData::Parse(const std::string &fname){

    GetPot gp_input(fname.c_str());

    // Case parameters:
    //-----------------------------
    case_title = gp_input("Case/title","test");
    case_postproc_dir = gp_input("Case/postprocess_dir","./results/");
    mesh_fname = gp_input("Case/mesh_file_name","./mesh/mesh.dat");
    grid_type_ = gp_input("Case/grid_type","Cartesian");
    generate_grid_flag_ = gp_input("Case/generate_grid_flag",0);
    Nnodes_x_ = gp_input("Case/Nnodes_x",1);
    Nnodes_y_ = gp_input("Case/Nnodes_y",1);
    x0_ = gp_input("Case/x_begin",-1.0);
    xf_ = gp_input("Case/x_end",1.0);
    y0_ = gp_input("Case/y_begin",-1.0);
    yf_ = gp_input("Case/y_end",1.0);
    uniform_ = gp_input("Case/uniform_grid",1);
    refine_level_ = gp_input("Case/refinement_level",0);
    //N_exact_plot_pts = gp_input("Case/N_exact_plot_pts",100);
    //N_uniform_pts_per_elem_ = gp_input("Case/N_uniform_pts_per_elem",2);
    //Npplot = N_uniform_pts_per_elem_;

    // Simulation parameters:
    //-----------------------------
    unsteady_data_print_flag_=gp_input("Simulation/unsteady_data_print_flag",0);
    unsteady_data_print_iter_=gp_input("Simulation/unsteady_data_print_iter",0);
    unsteady_data_print_time_=gp_input("Simulation/unsteady_data_print_time",1.0);
    restart_iter_ = gp_input("Simulation/restart_iter",0);
    restart_flag = gp_input("Simulation/restart_flag",0);
    Sim_mode = gp_input("Simulation/mode","normal");
    case_no_ = gp_input("Simulation/case_no","00");

    // Wave parameters:
    //----------------------------
    a_wave_x_ = gp_input("wave/x_wave_speed",1.0);
    a_wave_y_ = gp_input("wave/y_wave_speed",1.0);
    wave_form_ = gp_input("wave/wave_form",0);
    // ./KannanWang2009:
    wave_freq_ = gp_input("wave/KannanWang2009/wave_freq",1.0);
    wave_amp_  = gp_input("wave/KannanWang2009/wave_amplitude",1.0);
    wave_shift = gp_input("wave/KannanWang2009/wave_freq_shift",0.0);
    exp_exponent_const = gp_input("wave/KannanWang2009/exp_exponent_const",1.0);
    //../Trigonometric:
    x_wave_type_=gp_input("wave/Trigonometric/x_wave_type","Sin");
    y_wave_type_=gp_input("wave/Trigonometric/y_wave_type","Sin");
    x_wave_freq_=gp_input("wave/Trigonometric/x_wave_freq",2.0);
    y_wave_freq_=gp_input("wave/Trigonometric/y_wave_freq",2.0);
    x_wave_phase=gp_input("wave/Trigonometric/x_wave_phase",0.0);
    y_wave_phase=gp_input("wave/Trigonometric/y_wave_phase",0.0);
    wave_amp_=gp_input("wave/Trigonometric/wave_amplitude",1.0);
    wave_const=gp_input("wave/Trigonometric/wave_const",0.0);
    //../Trigonometric2:
    wave_type_=gp_input("wave/Trigonometric2/wave_type","Sin");
    // ../Gaussian:
    Gaussian_amp_ = gp_input("wave/Gaussian/Gaussian_amplitude",1.0);
    Gaussian_exponent_ = gp_input("wave/Gaussian/Gaussian_exponent",-50.0);
    //../linear
    slope_x = gp_input("wave/linear/slope_x",1.0);
    slope_y = gp_input("wave/linear/slope_y",1.0);
    lin_const_ = gp_input("wave/linear/lin_const_",0.0);

    // Space Solver parameters:
    //-----------------------------
    eqn_set = gp_input("space_solver/eqn_set","Diffusion");
    eqn_type_ = gp_input("space_solver/eqn_type","linear_diffus");
    if(eqn_set=="Advection"){
        if(eqn_type_!="linear_advec")
            FatalError_exit("eqn type is not compatible");
    }else if(eqn_set=="Diffusion"){
        if(eqn_type_!="linear_diffus")
            FatalError_exit("eqn type is not compatible");
    }else if(eqn_set=="Advection_Diffusion"){
        if(eqn_type_!="linear_advec_diffus")
            FatalError_exit("eqn type is not compatible");
    }else{
        FatalError_exit("eqn set is not implemented");
    }
    poly_order_=gp_input("space_solver/polynomial_order",1);
    schemes_data_dir_=gp_input("space_solver/schemes_data_dir","./schemes_data/");
    // ./advec_eqn:
    upwind_param_=gp_input("space_solver/advec_eqn/upwind_param",1.0);
    // ./heat_eqn:
    thermal_diffus
            = gp_input("space_solver/heat_eqn/thermal_diffusivity",1.0);
    penalty_param_
            = gp_input("space_solver/heat_eqn/penalty_param",1.0);
    diffus_scheme_type_ =
            gp_input("space_solver/heat_eqn/diffus_scheme_type","BR2");

    // Time Solver parameters:
    //--------------------------------
    calc_dt_flag = gp_input("time_solver/calculate_dt_flag",1);
    calc_dt_adv_diffus_flag = gp_input("time_solver/calc_dt_adv_diffus_flag",0);
    CFL_ = gp_input("time_solver/CFL_no",1.0);
    dt_ = gp_input("time_solver/dt",1e-9);
    t_init_ = gp_input("time_solver/initial_time",0.0);
    t_end_ = gp_input("time_solver/final_time",1.0);
    maxIter_ = gp_input("time_solver/maximum_iteration",1e9);
    end_of_sim_flag_ = gp_input("time_solver/end_of_simulation_flag",0);
    Nperiods = gp_input("time_solver/no_of_periods",1.0);
    // ./explicit:
    RK_order_=gp_input("time_solver/explicit/RK_order",1);

    return;
}
