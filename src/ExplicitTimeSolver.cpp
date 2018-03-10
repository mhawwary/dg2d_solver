#include"ExplicitTimeSolver.hpp"

ExplicitTimeSolver::ExplicitTimeSolver(){

 return;
}

ExplicitTimeSolver::~ExplicitTimeSolver(){

    Reset_time_solver();

 return;
}

void ExplicitTimeSolver::setupTimeSolver(SpaceSolver *o_space_solver_
                                         , SimData *o_simdata_){
    register int i;

    space_solver = o_space_solver_;
    simdata = o_simdata_;
    Ndof = space_solver->GetNdof();
    Nelem=space_solver->GetNelem();

    resid = new double[Ndof*Nelem];
    dt_ = space_solver->GetTimeStep();
    IterNo=0;

    if(simdata->RK_order_>1){
        q_temp = new double[Nelem*Ndof];

        if(simdata->RK_order_==4){
            resid_temp0 = new double [Nelem*Ndof];
            resid_temp1 = new double [Nelem*Ndof];
            resid_temp2 = new double [Nelem*Ndof];
        }
    }

    return;
}

void ExplicitTimeSolver::SolveOneStep(double *qn_){

    switch (simdata->RK_order_) {

    case 1:

        FwdEuler(qn_);

        break;

    case 2:

        SSPRK22(qn_);

        break;

    case 3:

        SSPRK33(qn_);

        break;

    case 4:

        classicRK4(qn_);

        break;

    default:
        char *ss=nullptr; ss= new char[100];
        sprintf(ss,"RK order of %d ",simdata->RK_order_);
        _notImplemented(ss);
        emptyarray(ss);
        break;
    }

    IterNo++;

//    if(IterNo==simdata->maxIter_-1)
//        dt_ = space_solver->GetLastTimeStep();

    return;
}

void ExplicitTimeSolver::CopyOldSol(double *q_t_, double *qn_){
    register int j;
    for(j=0; j<Nelem*Ndof; j++)  q_t_[j] = qn_[j];
    return;
}

void ExplicitTimeSolver::CopyOldResid(double *resid_t_, double *old_resid_){
    register int j;
    for(j=0; j<Nelem*Ndof; j++)  resid_t_[j] = old_resid_[j];
    return;
}

void ExplicitTimeSolver::FwdEuler(double *q_){
    register int i;
    for(i=0; i<Nelem*Ndof; i++)
            q_[i] = q_[i] + dt_ * resid[i];

    space_solver->UpdateResid(resid,q_);

    return;
}

void ExplicitTimeSolver::SSPRK22(double *q_){

    register int j;
    CopyOldSol(q_temp,q_);  // Copying level n solution and saving it
    // Step1:
    //-----------
    for(j=0; j<Nelem*Ndof; j++)
        q_[j] = q_temp[j] + dt_ * resid[j];

    space_solver->UpdateResid(resid,q_);

    // Step2:
    //------------
    for(j=0; j<Nelem*Ndof; j++)
        q_[j] = 0.5 * ( q_temp[j] +  q_[j] + dt_ * resid[j] );

    space_solver->UpdateResid(resid,q_);

    return;
}

void ExplicitTimeSolver::SSPRK33(double *q_){

    register int j;
    CopyOldSol(q_temp,q_);  // Copying level n solution and saving it

    // Step1:
    //-----------
    for(j=0; j<Nelem*Ndof; j++)
        q_[j] = q_temp[j] + dt_ * resid[j];

    space_solver->UpdateResid(resid,q_);

    // Step2:
    //------------
    for(j=0; j<Nelem*Ndof; j++)
        q_[j] =  (0.75 * q_temp[j] )
                + 0.25 * ( q_[j] + dt_ * resid[j] );

    space_solver->UpdateResid(resid,q_);

    // Step3:
    //--------------
    for(j=0; j<Nelem*Ndof; j++)
        q_[j] =  ( q_temp[j]/3. )
                + 2. * ( q_[j] + dt_ * resid[j] ) / 3.;

    space_solver->UpdateResid(resid,q_);

    return;
}

void ExplicitTimeSolver::classicRK4(double *q_){

    register int j;
    CopyOldSol(q_temp,q_);  // Copying level 0 solution and saving it
    CopyOldResid(resid_temp0,resid);  // Copying level 0 residual and saving it

    // Step1:
    //-----------
    for(j=0; j<Nelem*Ndof; j++)
        q_[j] = q_temp[j] + 0.5 * dt_ * resid[j];

    space_solver->UpdateResid(resid,q_);

    CopyOldResid(resid_temp1,resid);  // Copying level 1 residual and saving it

    // Step2:
    //------------
    for(j=0; j<Nelem*Ndof; j++)
        q_[j] = q_temp[j] + 0.5 * dt_ * resid[j];

    space_solver->UpdateResid(resid,q_);

    CopyOldResid(resid_temp2,resid);  // Copying level 2 residual and saving it

    // Step3:
    //--------------
    for(j=0; j<Nelem*Ndof; j++)
        q_[j] = q_temp[j] + dt_ * resid[j];

    space_solver->UpdateResid(resid,q_);

    // Step4:
    //--------------
    for(j=0; j<Nelem*Ndof; j++)
        q_[j] = q_temp[j] + (dt_/6.) * ( resid_temp0[j] + 2*resid_temp1[j]
                                         + 2*resid_temp2[j] + resid[j] );

    space_solver->UpdateResid(resid,q_);

    return;
}

void ExplicitTimeSolver::ComputeInitialResid(double *qn_){
    space_solver->UpdateResid(resid,qn_);
    return;
}

void ExplicitTimeSolver::Reset_time_solver(){
    emptyarray(resid);
    emptyarray(q_temp);
    emptyarray(resid_temp0);
    emptyarray(resid_temp1);
    emptyarray(resid_temp2);

    return;
}

