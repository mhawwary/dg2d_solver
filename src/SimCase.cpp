#include"SimCase.hpp"
//#include"general_tools.h"

SimCase::SimCase(void){
    emptypointer(grid);
    emptypointer(grid_data);
}

SimCase::~SimCase(void){
    emptypointer(dg_solver_);
    emptypointer(grid);
    emptypointer(grid_data);
    emptypointer(time_solver_);
}

void SimCase::setup(const std::string &input_fname_){

    logo();
    // Parsing input data
    simdata.Parse(input_fname_);
    setup_output_directories();
    copy_problem_inputdata();

    if(simdata.generate_grid_flag_==1){  // Generating Cartesian Grid
        SimplyConnectGrid tempGrid_;
        tempGrid_.Setup(simdata.Nnodes_x_,simdata.Nnodes_y_
                       ,simdata.x0_,simdata.xf_
                       ,simdata.y0_,simdata.yf_
                       ,simdata.case_postproc_dir.c_str());
        tempGrid_.ConstructGrid();
        simdata.mesh_fname = tempGrid_.getMeshfname();
    }

    // Allocating Solvers:
    if(simdata.eqn_set=="Diffusion")
        dg_solver_ = new DG_2D_DIFFUS_SOLVER;
    else
        _notImplemented("Equation set");

    time_solver_ = new ExplicitTimeSolver;

    return;
}

void SimCase::setup_output_directories(){

    // Setting up some directories:
    //---------------------------------
    allocator<char> allchar; // default allocator for char
    char case_dir[70];

    sprintf(case_dir,"p%dRK%d_%s/",simdata.poly_order_
            ,simdata.RK_order_,simdata.diffus_scheme_type_.c_str());

    char *current_working_dir=allchar.allocate(1500);
    getcwd(current_working_dir,1500);
    mkdir(simdata.case_postproc_dir.c_str(),0777);  // parent directory if not already exists
    chdir(simdata.case_postproc_dir.c_str());
    mkdir(simdata.case_title.c_str(),0777);
    chdir(simdata.case_title.c_str());
    mkdir(case_dir,0777);
    chdir(case_dir);

    mkdir("./field_output",0777);
    mkdir("./errors",0777);
    mkdir("./mesh_metrics",0777);
    chdir(current_working_dir);
    simdata.case_postproc_dir+=simdata.case_title+"/"+case_dir;

    cout<<"\n--> Currnet working directory: "
        <<current_working_dir<<endl;
    cout<<"--> Post processing directory: "
        <<simdata.case_postproc_dir.c_str()<<endl;

    emptyarray(current_working_dir);

    return;
}

void SimCase::InitSim(){

    // Preparing MeshData:
    //------------------------
    cout << "Mesh Processing..........";
    grid = new Mesh;
    grid->Read(simdata.mesh_fname);
    grid->generate_meshData();
    std::string write_fname_;
    std::string write_fname_tec;
    write_fname_=simdata.case_postproc_dir+"mesh_metrics/grid_info.dat";
    write_fname_tec=simdata.case_postproc_dir+"mesh_metrics/gridtest.plt";
    grid->WriteMesh(write_fname_);
    grid->WriteMeshTecplot(write_fname_tec);
    grid_data = grid->Release_meshData();
    emptypointer(grid);
    cout<<"done!" <<endl;
    printf("Ntriangles: %d,  Nquads: %d,  Npolygons: %d\n"
           ,grid_data->NtriElem,grid_data->NquadElem,grid_data->Npolygon);
    //Setupp space & time solvers and initialize the solution:
    //-----------------------------------------------------------------
    dg_solver_->setup_solver(grid_data,simdata);
    dg_solver_->InitSol();
    PostProcess(0);
    time_solver_->setupTimeSolver(dg_solver_,&simdata);
//    simdata.dump_python_inputfile();

    return;
}

void SimCase::RunSim(){

    int n_iter_print;
    int local_iter=0;
    double gtime = dg_solver_->GetPhyTime();
    double dt_= dg_solver_->GetTimeStep();
    double dt_last_print=0.0;
    double temp_tol=1e-8;

    if(simdata.unsteady_data_print_flag_==0){      // use iter_print to print
        n_iter_print = simdata.unsteady_data_print_iter_;
        if(n_iter_print<=1)
            FatalError_exit("Warning: iter to print is very small <=1 ");
        dt_last_print = dt_;
        n_iter_print--;

    }else if(simdata.unsteady_data_print_flag_==1){  // use time point to print
        if(simdata.unsteady_data_print_time_ < (dt_ - temp_tol))
            FatalError_exit("Warning:  time to print is less than dt");

        n_iter_print= (int) round( simdata.unsteady_data_print_time_/ dt_) ;

        if((n_iter_print*dt_) > (simdata.unsteady_data_print_time_-temp_tol) ){
            n_iter_print--;
            dt_last_print = simdata.unsteady_data_print_time_ - (n_iter_print * dt_);
        }else if((n_iter_print*dt_) < (simdata.unsteady_data_print_time_+temp_tol) ){
            dt_last_print = simdata.unsteady_data_print_time_ - (n_iter_print*dt_);
        }

        if(n_iter_print<=1)
            FatalError_exit("Warning: iter to print is very small <=1 ");

    }else if(simdata.unsteady_data_print_flag_==2){
        // print using the specified iter_print without dt changing except the last one
        n_iter_print=simdata.unsteady_data_print_iter_;
        dt_last_print=dt_;
    }else{
        FatalError_exit("unsteady data print flag error");
    }

    printf("\nnIter to print unsteady data: %d, dt_last: %1.5e\n"
           ,n_iter_print, dt_last_print);

    // First Solve:
    time_solver_->ComputeInitialResid(dg_solver_->GetNumSol());
    time_solver_->SolveOneStep(dg_solver_->GetNumSol());
    time_solver_->space_solver->UpdatePhyTime(dt_);
    gtime=dg_solver_->GetPhyTime();
    local_iter++;

    if(n_iter_print==1){
        printf("\nIter No:%d, time: %f",time_solver_->GetIter(),gtime);
        PostProcess(time_solver_->GetIter());
        local_iter=0;
    }

    // main solution loop:
    if(simdata.unsteady_data_print_flag_==0
            || simdata.unsteady_data_print_flag_==1){
        while ( gtime < (simdata.t_end_-(1+1e-5)*(dt_+1e-10))){
            time_solver_->SolveOneStep(dg_solver_->GetNumSol());
            time_solver_->space_solver->UpdatePhyTime(dt_);
            gtime=dg_solver_->GetPhyTime();
            local_iter++;

            if(local_iter%n_iter_print==0){
                time_solver_->Set_time_step(dt_last_print);
                time_solver_->SolveOneStep(dg_solver_->GetNumSol());
                time_solver_->space_solver->UpdatePhyTime(dt_last_print);
                gtime=dg_solver_->GetPhyTime();
                printf("\nIter No:%d, time: %1.5f",time_solver_->GetIter(),gtime);
                PostProcess(time_solver_->GetIter());
                time_solver_->Set_time_step(dt_);
                //time_solver_->Reset_iter(time_solver_->GetIter()-1);
                local_iter=0;
            }
        }
        //PostProcess(time_solver_->GetIter());

    }else if(simdata.unsteady_data_print_flag_==2){
        while ( fabs(gtime - simdata.t_end_) > (dt_+temp_tol) ){

            time_solver_->SolveOneStep(dg_solver_->GetNumSol());
            time_solver_->space_solver->UpdatePhyTime(dt_);
            gtime=dg_solver_->GetPhyTime();
            local_iter++;

            if(local_iter%n_iter_print==0){
                printf("\nIter No:%d, time: %1.5f",time_solver_->GetIter(),gtime);
                PostProcess(time_solver_->GetIter());
                local_iter=0;
            }
        }

        // Last iteration:
        dt_last_print = dg_solver_->GetLastTimeStep();
        time_solver_->Set_time_step(dt_last_print);
        time_solver_->SolveOneStep(dg_solver_->GetNumSol());
        if(dt_last_print>=temp_tol)
            time_solver_->space_solver->UpdatePhyTime(dt_last_print);
        else
            time_solver_->space_solver->UpdatePhyTime(dt_);

        PostProcess(time_solver_->GetIter());
    }
    dump_exactsol_field_data(dg_solver_->GetVertexExactSol()
                          ,time_solver_->GetIter()
                          ,simdata.case_postproc_dir,grid_data);
    double L1_proj=dg_solver_->GetL1projSolerror();
    double L2_proj=dg_solver_->GetL2projSolerror();
    std::cout<<"\nL1_proj: "<<L1_proj<<"   "<<"L2_proj: "<<L2_proj<<"\n";
    dg_solver_->dump_errors(L1_proj,L2_proj);
    return;
}

void SimCase::PostProcess(const int& iter_){

    if(iter_==simdata.maxIter_){
//        dump_exactsol_field_data(dg_solver_->GetVertexExactSol()
//                              ,time_solver_->GetIter()
//                              ,simdata.case_postproc_dir,grid_data);
//        double L1_proj=dg_solver_->GetL1projSolerror();
//        double L2_proj=dg_solver_->GetL2projSolerror();
//        std::cout<<"\nL1_proj: "<<L1_proj<<"   "<<"L2_proj: "<<L2_proj<<"\n";
//        dg_solver_->dump_errors(L1_proj,L2_proj);
    }

    dg_solver_->Compute_vertex_sol();
    dump_field_data(dg_solver_->GetVertexNumSol()
                    ,iter_
                    ,simdata.case_postproc_dir,grid_data);

    return;
}

void SimCase::logo() {

    cout<<"                                                                                         "<<endl;
    cout<<"_________________________________________________________________________________________"<<endl;
    cout<<"                                                                                         "<<endl;
    cout<<"                 "<<"  Welcome to the Discontinuous Galerkin solver   "<<"               "<<endl;
    cout<<"                       "<<"  for 2D scalar conservation laws   "<<"                      "<<endl;
    cout<<"                                                                                         "<<endl;
    cout<<"       Author:               Mohammad Alhawwary, PhD. Student                            "<<endl;
    cout<<"  Affiliation:   Aerospace Engineering Department, University of Kansas, USA             "<< endl;
    cout<<"                                                                                         "<<endl;
    cout<<"_______________________________________03/01/2018________________________________________"<<endl;
    return ;
}

void SimCase::copy_problem_inputdata(){

    char dirr[100];

    sprintf(dirr,"%s/case_input.in",simdata.case_postproc_dir.c_str());

    copyFile("./input/case_input.in",dirr);

    return;
}

void SimCase::copyFile(const std::string& fileNameFrom
                       , const std::string& fileNameTo)
{
     std::ifstream in (fileNameFrom.c_str());
     std::ofstream out (fileNameTo.c_str());
     out << in.rdbuf();
     out.close();
     in.close();
}
