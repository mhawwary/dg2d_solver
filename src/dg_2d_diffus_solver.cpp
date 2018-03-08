#include"dg_2d_diffus_solver.hpp"

// Constructor/Destructor/ Setup functions:
//------------------------------------------------

DG_2D_DIFFUS_SOLVER::~DG_2D_DIFFUS_SOLVER(void){

    Reset_solver();
}

void DG_2D_DIFFUS_SOLVER::Reset_solver(){

    emptyarray(Qn);
    emptyarray(Qv_exact);
    emptyarray(Qv);
    emptyarray(Qex_proj);
    emptyarray(Lk_norm2sq);

    return;
}

void DG_2D_DIFFUS_SOLVER::setup_solver(MeshData*& meshdata_, SimData& osimdata_){

    grid_ = meshdata_;
    simdata_ = &osimdata_;

    eta_face = simdata_->penalty_param_;
    Porder = simdata_->poly_order_;
    Ndof_1D = (Porder+1);
    Ndof = Ndof_1D*Ndof_1D;

    // Setup quadrature orders
    switch (Ndof_1D) {
    case 1:  // p0
        Nquad_1d = 1;
        break;
    case 2: // p1
        Nquad_1d = 2;
        break;
    case 3:  // p2
        Nquad_1d = 3;
        break;
    case 4:  // p3
        Nquad_1d = 4;
        break;
    case 5:  // p4
        Nquad_1d = 5;
        break;
    case 6:  // p5
        Nquad_1d = 6;
        break;
    default:
        break;
    }
    quad_1d_.setup_quadrature(Nquad_1d);
    quad_2d_.setup_quadrature(Nquad_1d);
    Nquad_2d=quad_2d_.Nq;

    Lk_norm2sq = new double[Ndof];

    Qn       = new double[grid_->Nelem*Ndof];
    Qex_proj = new double[grid_->Nelem*Ndof];
    Qv_exact = new double[grid_->Nnodes];
    Qv = new double[grid_->Nnodes];

    SetPhyTime(simdata_->t_init_);

    CalcTimeStep();
    load_scheme_matrices();
    ComputeExactSolShift();
    Compute_projected_exact_sol();
    Compute_exact_vertex_sol();

    return;
}

//Basic members:
void DG_2D_DIFFUS_SOLVER::load_scheme_matrices(){

    int i;
    std::string fname_;
    std::string scheme_dir_;

    if(simdata_->diffus_scheme_type_=="BR2")
        scheme_dir_=simdata_->schemes_data_dir_+"BR2/Quad/p"
                +std::to_string(Porder)+"/";
    else if(simdata_->diffus_scheme_type_=="SIPG")
        scheme_dir_=simdata_->schemes_data_dir_+"SIPG/Quad/p"
                +std::to_string(Porder)+"/";

    // Reading mass matrix M
    fname_ = scheme_dir_+"M.dat";
    std::ifstream input;
    open_inputfile_forreading(fname_,input);

    for(i=0; i<Ndof; i++){
        input >> Lk_norm2sq[i];
        Lk_norm2sq[i]=4.0*Lk_norm2sq[i];
        /* this is because the dumped Mass matrix
         * was integrated with multiplication by the
         * Jacobian=h^2/4 and without the h^2
         */
    }
    input.close();

   // Reading scheme update matrices
   fname_ = scheme_dir_+"K0.dat";
   read_sparsematrix(K0,fname_);
   fname_ = scheme_dir_+"K1.dat";
   read_sparsematrix(K1,fname_);
   fname_ = scheme_dir_+"K2.dat";
   read_sparsematrix(K2,fname_);
   fname_ = scheme_dir_+"K3.dat";
   read_sparsematrix(K3,fname_);
   fname_ = scheme_dir_+"L.dat";
   read_sparsematrix(L,fname_);
   fname_ = scheme_dir_+"Q.dat";
   read_sparsematrix(Q,fname_);
   fname_ = scheme_dir_+"Qp0.dat";
   read_sparsematrix(Q0,fname_);
   fname_ = scheme_dir_+"Qp1.dat";
   read_sparsematrix(Q1,fname_);
   fname_ = scheme_dir_+"Qp2.dat";
   read_sparsematrix(Q2,fname_);
   fname_ = scheme_dir_+"Qp3.dat";
   read_sparsematrix(Q3,fname_);

   return;
}

void DG_2D_DIFFUS_SOLVER::read_sparsematrix(MYSPARSE_MATRIX& smatrix_
                                             ,std::string& readfname_){

    double h2=0;
    h2 = 1.0/grid_->elemlist[0].Vc; // 1/h^2
    int i,j;
    double temp_val=0;
    int nnz_row_index_;
    std::ifstream input;
    open_inputfile_forreading(readfname_,input);

    for(i=0; i<Ndof; i++){
        nnz_row_index_=0;
        temp_val=0;
        for(j=0; j<Ndof; j++){
            input >> temp_val;
            if(temp_val!=0){
                temp_val = temp_val * h2; // 1/h^2
                smatrix_.set_value(temp_val,i,j);
                nnz_row_index_++;
                temp_val=0;
            }
        }
        smatrix_.set_nnz_row(nnz_row_index_);
    }
    input.close();
    smatrix_.set_size();

    return;
}

double DG_2D_DIFFUS_SOLVER::eval_1dbasis_poly(const double &o_xi_, const int &basis_k_){

    double Lk_=0.0;

    switch (basis_k_) {
    case 0:
        Lk_ = 1.0;  break;
    case 1:
        Lk_ = o_xi_;  break;
    case 2:
        Lk_ = 0.5 * (3.0 * o_xi_*o_xi_ -1.0);  break;
    case 3:
        Lk_ = 0.5 * (5.0 * pow(o_xi_,3) - 3.0 * o_xi_);  break;
    case 4:
        Lk_ = (35. * pow(o_xi_,4) - 30. * o_xi_*o_xi_ + 3. ) /8. ;  break;
    case 5:
        Lk_ = (63.0 * pow(o_xi_,5) - 70.0 * pow(o_xi_,3) + 15.0 * o_xi_ ) /8. ;  break;
    default:
        FatalError_exit("k basis exceeds poly order");
        break;
    }

    return Lk_;
}

double DG_2D_DIFFUS_SOLVER::eval_2dbasis_poly(const double &o_xi_
                                            , const double &o_eta_
                                            , const int &basis_k_){
    if(basis_k_>=Ndof){
        FatalError_exit("basis_k >= Ndof, eval_2d_basis");
        return 0.0;
    }else{
        int i,j;
        double Lk_=0.0;

        //1D basis functions:
        double o_phi_xi_=0.0;
        double o_phi_eta_=0.0;
        //evaluate 2D basis
        i = int(basis_k_/(Porder+1));  // row 1D position for eta dir
        j = basis_k_%(Porder+1);       // column 1D position for xi dir

        o_phi_xi_ =eval_1dbasis_poly(o_xi_,j);
        o_phi_eta_=eval_1dbasis_poly(o_eta_,i);

        Lk_ = o_phi_xi_*o_phi_eta_;
        return Lk_;
    }
}

double DG_2D_DIFFUS_SOLVER::eval_bilinear_mapshapefunc(const double &o_xi_
                                                     ,const double &o_eta_
                                                     ,const int &o_index_){
    double o_epsi_=0.0;

    switch (o_index_) {
    case 0:
        o_epsi_ = 0.25*(1.0 - o_xi_)*(1.0-o_eta_);
        break;
    case 1:
        o_epsi_ = 0.25*(1.0 + o_xi_)*(1.0-o_eta_);
        break;
    case 2:
        o_epsi_ = 0.25*(1.0 + o_xi_)*(1.0+o_eta_);
        break;
    case 3:
        o_epsi_ = 0.25*(1.0 - o_xi_)*(1.0+o_eta_);
        break;
    default:
        FatalError_exit("Wrong index for mapping functions");
        break;
    }

    return o_epsi_;
}

void DG_2D_DIFFUS_SOLVER::eval_local_xy_coord(const int& o_eID_
                                         ,const double &xi_
                                         , const double &eta_
                                          ,double& o_x_
                                          ,double& o_y_){
    int nodeID_;
    double xn=0.0;
    double yn=0.0;
    double o_epsi_=0.0;
    o_x_=0.0;
    o_y_=0.0;

    register int i;
    for(i=0; i<grid_->elemlist[o_eID_].n_local_nodes; i++){
        nodeID_ = grid_->elemlist[o_eID_].to_node[i] ;
        xn = grid_->Xn[nodeID_];
        yn = grid_->Yn[nodeID_];
        o_epsi_=eval_bilinear_mapshapefunc(xi_,eta_,i);
        o_x_ += o_epsi_ * xn;
        o_y_ += o_epsi_ * yn;
    }

    return;
}

//Solver Computation functions:
void DG_2D_DIFFUS_SOLVER::ComputeExactSolShift(){
    // Preparing shift information
    double a=0.;
    x_wave_length_ = simdata_->xf_ - simdata_->x0_ ;
    a = simdata_->a_wave_x_;  // should be changed to a_wave_x
    x_exact_sol_shift = (a * simdata_->t_end_ );
    y_wave_length_ = simdata_->yf_ - simdata_->y0_ ;
    a = simdata_->a_wave_y_;
    y_exact_sol_shift = (a * simdata_->t_end_ );

    x_exact_sol_shift=0.0;
    y_exact_sol_shift=0.0;

    return;
}

void DG_2D_DIFFUS_SOLVER::InitSol(){

    register int j;
    int k=0;

    if(simdata_->eqn_type_=="linear_diffus"){ // linear diffusion equation
        for(j=0; j<grid_->Nelem; j++)
            for(k=0; k<Ndof; k++)
                Qn[j*Ndof+k] = initSol_legendre_proj(j,k,quad_2d_);
    }else{
        FatalError_exit("eqn_type is not implemented");
    } 

    return;
}

double DG_2D_DIFFUS_SOLVER::initSol_legendre_proj(const int &eID,
                                       const int &basis_k,
                                        const GaussQuad2D &o_quad2d_){

    double xx=0.0,yy=0.0;
    double II=0.0;
    double Qinit_=0.0;
    double Lk_=0.0;

    for (int i=0; i<o_quad2d_.Nq; i++){
        eval_local_xy_coord(eID,o_quad2d_.Gaus_pts[0][i]
                             ,o_quad2d_.Gaus_pts[1][i],xx,yy);
        Qinit_= eval_init_sol(xx,yy);
        Lk_ = eval_2dbasis_poly(o_quad2d_.Gaus_pts[0][i]
                             ,o_quad2d_.Gaus_pts[1][i], basis_k);
        II += o_quad2d_.Gaus_wts[i] * Qinit_ * Lk_ ;
    }
    II = II / Lk_norm2sq[basis_k] ;

    return II;
}

double DG_2D_DIFFUS_SOLVER::exactSol_legendre_proj(const int &eID,
                                       const int &basis_k,
                                        const GaussQuad2D &o_quad2d_){

    double xx=0.0,yy=0.0;
    double II=0.0;
    double Qexx_=0.0;
    double Lk_=0.0;

    for(int i=0; i<o_quad2d_.Nq; i++){
        eval_local_xy_coord(eID,o_quad2d_.Gaus_pts[0][i]
                             ,o_quad2d_.Gaus_pts[1][i],xx,yy);
        xx -= x_exact_sol_shift;
        yy -= y_exact_sol_shift;

        if(simdata_->wave_form_==0)
            Qexx_=eval_init_sol(xx,yy);
        else if(simdata_->wave_form_==1)
            Qexx_= eval_exact_trigsol_diffus(xx,yy,simdata_->t_end_);

        Lk_ = eval_2dbasis_poly(o_quad2d_.Gaus_pts[0][i]
                             ,o_quad2d_.Gaus_pts[1][i], basis_k);

        II += o_quad2d_.Gaus_wts[i] * Qexx_ * Lk_ ;
    }
    II = II / Lk_norm2sq[basis_k] ;

    return II;
}

double DG_2D_DIFFUS_SOLVER::eval_init_sol(const double &xx_, const double&yy_){

    double Qinit_=0.0;

    if(simdata_->wave_form_==0){ //KannanWang2009
        Qinit_ = simdata_->wave_amp_
                * sin(simdata_->wave_freq_*xx_+simdata_->wave_shift)
                * exp(simdata_->exp_exponent_const*yy_);

    }else if(simdata_->wave_form_==1){ //Trigonometric, periodic
        double o_x_wave_=0.0,o_y_wave_=0.0;
        if(simdata_->x_wave_type_=="sin")
            o_x_wave_=sin(simdata_->x_wave_freq_*PI*xx_+simdata_->x_wave_phase);
        else if(simdata_->x_wave_type_=="cos")
            o_x_wave_=cos(simdata_->x_wave_freq_*PI*xx_+simdata_->x_wave_phase);
        else
            _notImplemented("x_wave_type");

        if(simdata_->y_wave_type_=="sin")
            o_y_wave_=sin(simdata_->y_wave_freq_*PI*yy_+simdata_->y_wave_phase);
        else if(simdata_->x_wave_type_=="cos")
            o_y_wave_=cos(simdata_->y_wave_freq_*PI*yy_+simdata_->y_wave_phase);
        else
            _notImplemented("y_wave_type");

        Qinit_ = simdata_->wave_amp_*o_x_wave_*o_y_wave_+simdata_->wave_const;
    }else{
        _notImplemented("wave form");
        FatalError_exit("not found waveform");
    }

    return Qinit_;
}

double DG_2D_DIFFUS_SOLVER::eval_exact_trigsol_diffus(const double &xx_
                                                      , const double &yy_
                                                      , const double &tt_){
    double Qexx_=0.0;
    if(simdata_->wave_form_==1){ //Trigonometric, periodic
        double o_x_wave_=0.0,o_y_wave_=0.0;
        if(simdata_->x_wave_type_=="sin")
            o_x_wave_=sin(simdata_->x_wave_freq_*PI*xx_+simdata_->x_wave_phase);
        else if(simdata_->x_wave_type_=="cos")
            o_x_wave_=cos(simdata_->x_wave_freq_*PI*xx_+simdata_->x_wave_phase);
        else
            _notImplemented("x_wave_type");

        if(simdata_->y_wave_type_=="sin")
            o_y_wave_=sin(simdata_->y_wave_freq_*PI*yy_+simdata_->y_wave_phase);
        else if(simdata_->x_wave_type_=="cos")
            o_y_wave_=cos(simdata_->y_wave_freq_*PI*yy_+simdata_->y_wave_phase);
        else
            _notImplemented("y_wave_type");
        double KK_=pow(simdata_->x_wave_freq_,2)
                + pow(simdata_->y_wave_freq_,2);
        Qexx_ = simdata_->wave_amp_*o_x_wave_*o_y_wave_*exp(-KK_*tt_)
                +simdata_->wave_const;
    }else{
        _notImplemented("wave form");
        FatalError_exit("not found waveform");
    }

    return Qexx_;
}

void DG_2D_DIFFUS_SOLVER::Compute_projected_exact_sol(){

    register int j;
    int k=0;
    for(j=0; j<grid_->Nelem; j++)
        for(k=0; k<Ndof; k++)
            Qex_proj[j*Ndof+k] = exactSol_legendre_proj(j,k,quad_2d_);
    return;
}

void DG_2D_DIFFUS_SOLVER::Compute_exact_vertex_sol(){

    register int i;
    double xx=0.0,yy=0.0;

    for(i=0; i<grid_->Nnodes; i++){
        xx=grid_->Xn[i]-x_exact_sol_shift;
        yy=grid_->Yn[i]-y_exact_sol_shift;
        if(simdata_->wave_form_==1){
            Qv_exact[i] = eval_exact_trigsol_diffus(xx,yy,simdata_->t_end_);
        }else{
            Qv_exact[i] = eval_init_sol(xx,yy);
        }
    }
    return;
}

void DG_2D_DIFFUS_SOLVER::CalcTimeStep(){

    double dx = (simdata_->xf_-simdata_->x0_)/(simdata_->Nnodes_x_-1);
    double dy = (simdata_->yf_-simdata_->y0_)/(simdata_->Nnodes_y_-1);
    double Vol_ = dx * dy;
    double radius_advec_x=0.0,radius_advec_y=0.0
            ,radius_advec_=0.0,radius_diffus_=0.0;

    T_period = (simdata_->xf_ - simdata_->x0_) / simdata_->a_wave_x_;
    radius_diffus_ = simdata_->thermal_diffus / Vol_ ;
    radius_advec_x =  simdata_->a_wave_x_ / dx  ;
    radius_advec_y =  simdata_->a_wave_y_ / dy  ;
    radius_advec_ = 0.5*(radius_advec_x+radius_advec_y); // average radius ( not accurate )

    if(simdata_->calc_dt_flag==1){   // use CFL as input

        CFL = simdata_->CFL_;
        if(simdata_->calc_dt_adv_diffus_flag==0)   // based on advection effect only
            time_step = CFL / radius_advec_ ;
        else if(simdata_->calc_dt_adv_diffus_flag==1)  // based on diffusion effect only
            time_step = CFL /  radius_diffus_ ;
        else if(simdata_->calc_dt_adv_diffus_flag==2)  // based on combined advection and diffusion effects
            time_step = CFL / ( radius_advec_ + radius_diffus_ );
        else
            FatalError_exit("Wrong Calc dt adv diffus flag");

        last_time_step = time_step;
        simdata_->dt_ = time_step;

    }else if(simdata_->calc_dt_flag==0){   // use dt as input

        time_step = simdata_->dt_;
        last_time_step = time_step;

        if(simdata_->calc_dt_adv_diffus_flag==0)
            CFL = time_step * radius_advec_  ;
        else if(simdata_->calc_dt_adv_diffus_flag==1)
            CFL = time_step * radius_diffus_ ;
        else if(simdata_->calc_dt_adv_diffus_flag==2)
            CFL = time_step * ( radius_advec_+ radius_diffus_ );
        else
            FatalError_exit("Wrong Calc dt adv diffus flag");

        simdata_->CFL_ = CFL;

    }else{

        FatalError_exit("Wrong Calc_dt_flag");
    }

    // Determining end of simulation parameters:
    //----------------------------------------------------
    double temp_tol=1e-8;
    if(simdata_->end_of_sim_flag_==0){  // use Nperiods as stopping criteria

        simdata_->t_end_ = simdata_->Nperiods * T_period;

        simdata_->maxIter_ = (int) floor(simdata_->t_end_/time_step);

        if((simdata_->maxIter_ * time_step)
                > (simdata_->Nperiods * T_period) ){

            last_time_step = simdata_->t_end_ - ((simdata_->maxIter_-1) * time_step);

        }else if((simdata_->maxIter_ * time_step)
                 < (simdata_->Nperiods * T_period) ){

            last_time_step = simdata_->t_end_ - (simdata_->maxIter_ * time_step);
        }

    }else if(simdata_->end_of_sim_flag_==1){ // use final time as stopping criteria

        simdata_->Nperiods = simdata_->t_end_/T_period;
        simdata_->maxIter_ = (int) floor(simdata_->t_end_/time_step);

        if((simdata_->maxIter_ * time_step)
                > (simdata_->t_end_-temp_tol) ){

            last_time_step = simdata_->t_end_ - ((simdata_->maxIter_-1) * time_step);

        }else if((simdata_->maxIter_ * time_step)
                 < (simdata_->t_end_+temp_tol) ){

            last_time_step = simdata_->t_end_ - (simdata_->maxIter_ * time_step);
        }

        if(last_time_step<=1e-10){
            last_time_step=time_step;
            simdata_->maxIter_--;
        }

    }else if(simdata_->end_of_sim_flag_==2){  // use Max Iter as stopping criteria

        simdata_->t_end_ = simdata_->maxIter_ * time_step;
        simdata_->Nperiods = simdata_->t_end_/T_period;

    }else{
        FatalError_exit("Wrong end_of_simulation_flag");
    }

    // Screen Output of input and simulation parameters:
    std::cout <<"\n===============================================\n";
    std::cout << "ThermalDiffusiv: "<< simdata_->thermal_diffus << std::endl;
    std::cout << "CFL no.        : "<<CFL<<std::endl;
    std::cout << "time step, dt  : "<<time_step<<std::endl;
    std::cout << "last_time_step : "<<last_time_step<<std::endl;
    std::cout << "input Nperiods : "<<simdata_->Nperiods<<std::endl;
    std::cout << "new   Nperiods : "<<simdata_->t_end_/T_period<<std::endl;
    std::cout << "x_exact_sol_shift: "<<x_exact_sol_shift<<std::endl;
    std::cout << "y_exact_sol_shift: "<<y_exact_sol_shift<<std::endl;
    std::cout << "T_period       : "<<T_period<<std::endl;
    printf("actual_end_time:%1.5f",simdata_->t_end_);
    std::cout <<"\nMax_iter: "<<simdata_->maxIter_<<std::endl;

    std::cout << "\nNumber of Elements: "<< grid_->Nelem<< std::endl;
    std::cout << "Polynomial  order : "<< Porder  << std::endl;
    std::cout << "Runge-Kutta order : "<< simdata_->RK_order_    << std::endl;
    std::cout << "Upwind parameter  : "<< simdata_->upwind_param_<< std::endl;
    std::cout << "Penalty parameter : "<< eta_face << std::endl;
    std::cout << "Poly GaussQuad order  : "<< Nquad_2d << std::endl;
    std::cout << "Jacobian: " <<grid_->elemlist[0].Jc<<std::endl;
    std::cout <<"===============================================\n";
}

double DG_2D_DIFFUS_SOLVER::evalSolution(const double *q_
                                         , const double &xi_pt_
                                         , const double &eta_pt_){
    double Q_=0.0;
    int k;
    for(k=0; k<Ndof; k++)
        Q_ += q_[k] * eval_2dbasis_poly(xi_pt_,eta_pt_,k);

    return Q_;
}

void DG_2D_DIFFUS_SOLVER::UpdateResid(double *o_resid_, double *o_qn_){

    for(register int i=0; i<grid_->Nelem; i++)
        UpdateResidOneCell(i,&o_resid_[i*Ndof],o_qn_);

    return;
}

void DG_2D_DIFFUS_SOLVER::UpdateResidOneCell(const int &eID_, double *resid_
                                             ,double *qn_){

    // note that qn_ here is the vector of total Nelem*Ndof data
    // while resid_ is the current element residual vector of Ndof
    double nx,ny;
    int fID_[4]={0,0,0,0};
    int iR_[4]={0,0,0,0};
    int i,j;
    double Ik;

    for(i=0; i<grid_->elemlist[eID_].n_local_faces; i++){
        fID_[i]=grid_->elemlist[eID_].to_face[i];
        nx = grid_->facelist[fID_[i]].nx;
        ny = grid_->facelist[fID_[i]].ny;
        iR_[i]=grid_->facelist[fID_[i]].Rcell;
        if(iR_[i]==eID_){
            iR_[i]=grid_->facelist[fID_[i]].Lcell;
            nx=-nx;
            ny=-ny;
        }
    }

    for(i=0; i<Ndof; i++)  resid_[i]=0.0;   // setting residual equal to zero

    //current element contribution:
    for(j=0; j<L.nnz_; j++){
        Ik = L.val[j] * qn_[eID_*Ndof+L.ix[1][j]];
        resid_[L.ix[0][j]] +=Ik;
    }
    for(j=0; j<Q.nnz_; j++){
        Ik = eta_face*Q.val[j] * qn_[eID_*Ndof+Q.ix[1][j]];
        resid_[Q.ix[0][j]] +=Ik;
    }

    //face0:
    i=0;
    for(j=0; j<K0.nnz_; j++){
        Ik = K0.val[j] * qn_[iR_[i]*Ndof+K0.ix[1][j]];
        resid_[K0.ix[0][j]] +=Ik;
    }
    for(j=0; j<Q0.nnz_; j++){
        Ik = eta_face*Q0.val[j] * qn_[iR_[i]*Ndof+Q0.ix[1][j]];
        resid_[Q0.ix[0][j]] +=Ik;
    }

    //face1:
    i=1;
    for(j=0; j<K1.nnz_; j++){
        Ik = K1.val[j] * qn_[iR_[i]*Ndof+K1.ix[1][j]];
        resid_[K1.ix[0][j]] +=Ik;
    }
    for(j=0; j<Q1.nnz_; j++){
        Ik = eta_face*Q1.val[j] * qn_[iR_[i]*Ndof+Q1.ix[1][j]];
        resid_[Q1.ix[0][j]] +=Ik;
    }
    //face2:
    i=2;
    for(j=0; j<K2.nnz_; j++){
        Ik = K2.val[j] * qn_[iR_[i]*Ndof+K2.ix[1][j]];
        resid_[K2.ix[0][j]] +=Ik;
    }
    for(j=0; j<Q2.nnz_; j++){
        Ik = eta_face*Q2.val[j] * qn_[iR_[i]*Ndof+Q2.ix[1][j]];
        resid_[Q2.ix[0][j]] +=Ik;
    }
    //face3:
    i=3;
    for(j=0; j<K3.nnz_; j++){
        Ik = K3.val[j] * qn_[iR_[i]*Ndof+K3.ix[1][j]];
        resid_[K3.ix[0][j]] +=Ik;
    }
    for(j=0; j<Q3.nnz_; j++){
        Ik = eta_face*Q3.val[j] * qn_[iR_[i]*Ndof+Q3.ix[1][j]];
        resid_[Q3.ix[0][j]] +=Ik;
    }

    for(i=0; i<Ndof; i++)  // multiplying by gamma
        resid_[i]=resid_[i]*simdata_->thermal_diffus;

    return;
}

void DG_2D_DIFFUS_SOLVER::Compute_vertex_sol(){

    register int i;
    int j,nID;
    double q_temp=0.0;
    double xi_pt_[4]={-1.0,1.0,1.0,-1.0};
    double eta_pt_[4]={-1.0,-1.0,1.0,1.0};

    for(i=0; i<grid_->Nnodes; i++)   Qv[i]=0.0;

    for(i=0; i<grid_->Nelem; i++){
        for(j=0; j<grid_->elemlist[i].n_local_nodes; j++){
            nID = grid_->elemlist[i].to_node[j];
            q_temp = evalSolution(&Qn[i*Ndof],xi_pt_[j],eta_pt_[j]);
            Qv[nID]+=q_temp;
        }
    }

    for(i=0; i<grid_->Nnodes; i++)   Qv[i]=Qv[i]/grid_->Nnode_neighElem[i];

    return;
}

double DG_2D_DIFFUS_SOLVER::L1_error_projected_sol(){

    register int j; int i;
    double L1_error=0.0,elem_error=0.0,II=0.0,VV_=0.0,q_ex,q_n;

    for(j=0; j<grid_->Nelem; j++){

        elem_error=0.0;
        for(i=0; i<quad_2d_.Nq; i++) {
            q_ex = evalSolution(&Qex_proj[j*Ndof]
                    ,quad_2d_.Gaus_pts[0][i], quad_2d_.Gaus_pts[1][i]);
            q_n = evalSolution(&Qn[j*Ndof]
                    ,quad_2d_.Gaus_pts[0][i], quad_2d_.Gaus_pts[1][i]);
            elem_error += quad_2d_.Gaus_wts[i] * fabs(q_ex - q_n);
        }

        II  += (grid_->elemlist[j].Jc * elem_error) ;
        VV_ += grid_->elemlist[j].Vc;
    }
    L1_error = II/VV_;

    return L1_error;
}

double DG_2D_DIFFUS_SOLVER::L2_error_projected_sol(){

    register int j; int i;
    double L2_error=0.0,elem_error=0.0,II=0.0,VV_=0.0,q_ex,q_n;

    for(j=0; j<grid_->Nelem; j++){

        elem_error=0.0;
        for(i=0; i<quad_2d_.Nq; i++) {
            q_ex = evalSolution(&Qex_proj[j*Ndof]
                    ,quad_2d_.Gaus_pts[0][i], quad_2d_.Gaus_pts[1][i]);
            q_n = evalSolution(&Qn[j*Ndof]
                    ,quad_2d_.Gaus_pts[0][i], quad_2d_.Gaus_pts[1][i]);
            elem_error += quad_2d_.Gaus_wts[i] * pow((q_ex - q_n),2);
        }

        II  += (grid_->elemlist[j].Jc * elem_error) ;
        VV_ += grid_->elemlist[j].Vc;
    }
    L2_error = sqrt(II/VV_);

    return L2_error;
}

void DG_2D_DIFFUS_SOLVER::dump_errors(double& L1_proj_sol_,double& L2_proj_sol_){

    char *fname=nullptr;
    fname = new char[150];

    if(simdata_->Sim_mode=="error_analysis_CFL"){
        sprintf(fname,"%serrors/errors_CFL%1.3e_eta%1.2f_t%1.3f.dat"
                ,simdata_->case_postproc_dir.c_str()
                ,CFL
                ,eta_face
                ,simdata_->t_end_);

        FILE* solerror_out=fopen(fname,"at+");
        fprintf(solerror_out, "%d %2.10e %2.10e\n",grid_->Nelem,L1_proj_sol_,L2_proj_sol_);
        fclose(solerror_out);
        emptyarray(fname);

    }else if(simdata_->Sim_mode=="error_analysis_dt"){
        sprintf(fname,"%serrors/errors_N%d_dt%1.3e_eta%1.2f_t%1.3f.dat"
                ,simdata_->case_postproc_dir.c_str()
                ,grid_->Nelem
                ,time_step
                ,eta_face
                ,simdata_->t_end_);

        FILE* solerror_out=fopen(fname,"w");
        fprintf(solerror_out, "%2.10e %2.10e\n",L1_proj_sol_, L2_proj_sol_);
        fclose(solerror_out);
        emptyarray(fname);

        // Dumping all errors in one file as a function of dt:
        //--------------------------------------------------------
        fname = new char[100];
        sprintf(fname,"%serrors/errors_N%d_alldt_eta%1.2f_t%1.3f.dat"
                ,simdata_->case_postproc_dir.c_str()
                ,grid_->Nelem
                ,eta_face
                ,simdata_->t_end_);

        solerror_out=fopen(fname,"at+");
        fprintf(solerror_out, "%1.7e %2.10e %2.10e\n",time_step,L1_proj_sol_,L2_proj_sol_);
        fclose(solerror_out);
        emptyarray(fname);

         // Dumping all errors in one file as a function of Nelem:
         //--------------------------------------------------------
         fname = new char[100];
         sprintf(fname,"%serrors/errors_dt%1.3e_eta%1.2f_t%1.3f.dat"
                 ,simdata_->case_postproc_dir.c_str()
                 ,time_step
                 ,eta_face
                 ,simdata_->t_end_);

         solerror_out=fopen(fname,"at+");
         fprintf(solerror_out, "%d %2.10e %2.10e\n",grid_->Nelem,L1_proj_sol_,L2_proj_sol_);
         fclose(solerror_out);
         emptyarray(fname);

    }else if( simdata_->Sim_mode=="test" || simdata_->Sim_mode=="normal" ){
        sprintf(fname,"%serrors/errors_N%d_CFL%1.3e_eta%1.2f_t%1.3f.dat"
                ,simdata_->case_postproc_dir.c_str()
                ,grid_->Nelem
                ,CFL
                ,eta_face
                ,simdata_->t_end_);

        FILE* solerror_out=fopen(fname,"w");
        fprintf(solerror_out, "%2.10e %2.10e\n",L1_proj_sol_,L2_proj_sol_);
        fclose(solerror_out);
        emptyarray(fname);

    }else if(simdata_->Sim_mode=="error_analysis_eta"){
        sprintf(fname,"%serrors/errors_N%d_CFL%1.3e_eta%1.2f_t%1.3f.dat"
                ,simdata_->case_postproc_dir.c_str()
                ,grid_->Nelem
                ,CFL
                ,eta_face
                ,simdata_->t_end_);

        FILE* solerror_out=fopen(fname,"w");
        fprintf(solerror_out, "%2.10e %2.10e\n",L1_proj_sol_,L2_proj_sol_);
        fclose(solerror_out);
        emptyarray(fname);

        // Dumping all errors in one file as a function of beta:
        //--------------------------------------------------------
        fname = new char[100];
        sprintf(fname,"%serrors/errors_N%d_CFL%1.3e_alleta_t%1.3f.dat"
                ,simdata_->case_postproc_dir.c_str()
                ,grid_->Nelem
                ,CFL
                ,simdata_->t_end_);

        solerror_out=fopen(fname,"at+");
        fprintf(solerror_out, "%1.2f %2.10e %2.10e\n",eta_face,L1_proj_sol_,L2_proj_sol_);
        fclose(solerror_out);
        emptyarray(fname);
    }

    return;
}


//Debug functions:
void DG_2D_DIFFUS_SOLVER::test_local_xy_coord(){

    double xx=0.0,yy=0.0;
    for(int j=0; j<grid_->Nelem; j++){
        std::cout<<"eID:"<<j<<" ";
        for (int i=0; i<quad_2d_.Nq; i++){
            eval_local_xy_coord(j,quad_2d_.Gaus_pts[0][i]
                    ,quad_2d_.Gaus_pts[1][i],xx,yy);

            std::cout<<"xi:"<<quad_2d_.Gaus_pts[0][i]<<" "
                     <<"eta:"<<quad_2d_.Gaus_pts[1][i]<<" "
                     <<"xx:"<<xx<<" "
                     <<"yy:"<<yy<<" "<<std::endl;
        }
        std::cin.get();
    }
    return;
}

void DG_2D_DIFFUS_SOLVER::test_eval_2dbasis(){

    double LL=0.0;
    double xi_pt_[4]={-1.0,1.0,1.0,-1.0};
    double eta_pt_[4]={-1.0,-1.0,1.0,1.0};
    int ii,jj;

    for(int j=0; j<4; j++){
        std::cout<<"xi:"<<xi_pt_[j]<<",   "
                 <<"eta:"<<eta_pt_[j]<<"\n";
        for(int k=0; k<Ndof; k++){
            LL=eval_2dbasis_poly(xi_pt_[j],eta_pt_[j],k);
            std::cout<<"L"<<k<<": "<<LL<<" "<<std::endl;
            ii = int(k/(Porder+1));  // row 1D position for eta dir
            jj = k%(Porder+1);       // column 1D position for xi dir

            std::cout<<"\nk:"<<k<<",  "
                     <<"eta_i="<<ii<<"  "
                     <<"xi_j="<<jj<<std::endl;

        }
    }
    std::cin.get();

    return;
}

void DG_2D_DIFFUS_SOLVER::Compute_exact_vertex_sol_usingproj(){

    register int i;
    int j,nID;
    double q_temp=0.0;
    double xi_pt_[4]={-1.0,1.0,1.0,-1.0};
    double eta_pt_[4]={-1.0,-1.0,1.0,1.0};

    for(i=0; i<grid_->Nnodes; i++)   Qv_exact[i]=0.0;

    for(i=0; i<grid_->Nelem; i++){
        for(j=0; j<grid_->elemlist[i].n_local_nodes; j++){
            nID = grid_->elemlist[i].to_node[j];
            q_temp = evalSolution(&Qex_proj[i*Ndof],xi_pt_[j],eta_pt_[j]);
            Qv_exact[nID]+=q_temp;
        }
    }

    for(i=0; i<grid_->Nnodes; i++)
        Qv_exact[i]=Qv_exact[i]/grid_->Nnode_neighElem[i];

    return;
}

void DG_2D_DIFFUS_SOLVER::TestGhostElements(){

    register int i;

    int iL,iR;
//    double Xc0,Xc1,Yc0,Yc1;

    for(i=0; i<grid_->Nbndfaces; i++){
        iL = grid_->facelist[i].Lcell;
        iR = grid_->facelist[i].Rcell;
        std::cout<<"fID:"<<grid_->facelist[i].ID<<"  "
                <<"iL:"<<iL<<", "<<grid_->elemlist[iL].bnd_type<<"  "
                <<"iR:"<<iR<<", "<<grid_->elemlist[iR].bnd_type<<std::endl;

        std::cin.get();

//        for(k=0; k<Ndof; k++)
//            Qn[iR*Ndof+k] = Qn[iR*Ndof+k];
    }
    return;
}
