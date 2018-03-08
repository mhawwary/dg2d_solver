#include "SimplyConnectGrid.hpp"

SimplyConnectGrid::SimplyConnectGrid(void){
    Reset();
}

SimplyConnectGrid::~SimplyConnectGrid(void){
    Reset();
}

void SimplyConnectGrid::Reset(){
    emptyarray(dump_dir);
    emptyarray(Im,x);
    emptyarray(Im,y);
    emptyarray(Nnodes,nodelist);
    emptyarray(Nfaces,facelist);
    emptyarray(Rcell);
    emptyarray(Lcell);
    emptyarray(normx);
    emptyarray(normy);
}

void SimplyConnectGrid::Setup(const int& oJm_, const int& oIm_
                                           ,const double& oxo_,const double& oxf_
                                           ,const double& oyo_,const double& oyf_
                                           ,const char* oprint_dir){

    dump_dir=new char[500];
    sprintf(dump_dir,"%s",oprint_dir);
    Im=oIm_;  // Nnodes in y-dir
    Jm=oJm_;  // Nnodes in x-dir
    Nnodes=Im*Jm;
    Ncells=(Im-1)*(Jm-1);
    Nfaces=2*Im*Jm-Im-Jm;

    xo_ = oxo_; xf_ = oxf_;
    yo_ = oyo_; yf_ = oyf_;

    register int i,j;
    // Grid coordinates arrays:
    x=new double*[Im];
    y=new double*[Im];
    for (i=0;i<Im;i++)
    {
        x[i]=new double[Jm];
        y[i]=new double[Jm];
    }

    // Nodes Data arrays:
    nodelist=new double*[Nnodes];
    for(i=0;i<Nnodes;i++)
        nodelist[i]=new double[2];

    // Faces Data arrays:
    facelist=new int*[Nfaces];
    for(i=0;i<Nfaces;i++)
        facelist[i]=new int[2];
    Lcell=new int[Nfaces];
    Rcell=new int[Nfaces];
    normx=new double[Nfaces];
    normy=new double[Nfaces];

//    ConstructGrid();

    return;
}

void SimplyConnectGrid::ConstructGrid(){

    CartesianGrid();   // Generate Cartesian Grid
    GridData(); // Grid basic connectivity data

    std::string tempfname_("cartmesh");
    tempfname_ += std::to_string(Im-1) + "x" + std::to_string(Jm-1);
    std::string write_fname_;

    write_fname_= dump_dir + tempfname_;
    WriteMesh(write_fname_);
    //WriteMeshTecplot(write_fname_);

    return;
}

void SimplyConnectGrid::CartesianGrid(){
    register int i,j;   double dx,dy,Lx,Ly;

    Lx = xf_-xo_;
    Ly = yf_-yo_;

    dx=Lx/(Jm-1);
    dy=Ly/(Im-1);

    for(i=0;i<Im;i++)
        for(j=0;j<Jm;j++)
            x[i][j]=(xo_+dx*j);
            //x[i][j]=(xo_+dx*j)/Lx; // normalized

    for(j=0;j<Jm;j++)
        for(i=0;i<Im;i++)
            y[i][j]=(yo_+dy*i);
            //y[i][j]=(yo_+dy*i)/Ly;   // normalized

    return;
}

void SimplyConnectGrid::GridData(){
    register int i=0,j,k=0,n=0,c=0;
    double dx,dy,L_;

    //********************//
    //NODES DATA          //
    //********************//
    for(i=0;i<Im;i++)
        for(j=0;j<Jm;j++){
            nodelist[k][0]=x[i][j];
            nodelist[k][1]=y[i][j];
            k++;
        }
    //******************//
    //  FACES DATA      //
    //******************//
    n=0;k=0;

    // Horizontal facelist
    //====================

    //Bottom boundary facelist ( -1 tag )
    for(j=0;j<Jm-1;j++)
    {
        facelist[n][0]=k;
        facelist[n][1]=k+1;
        Rcell[n]=-1;
        Lcell[n]=n;
        n++; k++;
    }
    k++;
    //Interior facelist ( 0 tag )
    for(i=1;i<Im-1;i++)
    {
        for(j=0;j<Jm-1;j++)
        {
            facelist[n][0]=k;
            facelist[n][1]=k+1;
            Rcell[n]=n+1-Jm;
            Lcell[n]=n;
            n++; k++;
        }
        k++;
    }
    //Top boundary facelist ( -4 tag )
    for(j=0;j<Jm-1;j++)
    {
        facelist[n][0]=k+1;
        facelist[n][1]=k;
        Lcell[n]=n+1-Jm;
        Rcell[n]=-4;
        n++; k++;
    }

    //Vertical facelist
    //====================
    k=0;
    //Left boundary facelist ( -2 tag )
    for(i=0;i<(Im-1);i++)
    {
        facelist[n][0]=k+Jm;
        facelist[n][1]=k;
        Lcell[n]=c;
        Rcell[n]=-2;
        n++; k=k+Jm; c=c+Jm-1;
    }
    //Interior Faces ( 0 tag )
    for(j=1;j<(Jm-1);j++)
    {
        k=j; c=j;
        for(i=0;i<(Im-1);i++)
        {
            facelist[n][0]=k;
            facelist[n][1]=k+Jm;
            Rcell[n]=c;
            Lcell[n]=c-1;
            n++; k=k+Jm; c=c+Jm-1;
        }
    }
    //Right boundary facelist ( -3 tag )
    k=Jm-1; c=Jm-1;
    for(i=0;i<(Im-1);i++)
    {
        facelist[n][0]=k;
        facelist[n][1]=k+Jm;
        Rcell[n]=-3;
        Lcell[n]=c-1;
        n++; k=k+Jm; c=c+Jm-1;
    }
    //*********************//
    //Normal Calculation   //
    //*********************//
    for(n=0;n<Nfaces;n++)
    {
        dx=nodelist[facelist[n][1]][0]-nodelist[facelist[n][0]][0];
        dy=nodelist[facelist[n][1]][1]-nodelist[facelist[n][0]][1];
        L_ = sqrt(dx*dx+dy*dy);
        normx[n]=dy/L_;
        normy[n]=-dx/L_;
    }

    return;
}

void SimplyConnectGrid::WriteMesh(std::string &owrite_fname_){

    meshfnamedir_ = owrite_fname_ + ".in";
    std::ofstream output(meshfnamedir_);

    output << Nnodes << " " << Nfaces << " " << Ncells ;

    register int i;

    for(i=0; i<Nnodes; i++)
        output <<"\n"<< nodelist[i][0] <<" "<< nodelist[i][1];
    for(i=0; i<Nfaces; i++)
        output <<"\n"<< facelist[i][0] <<" "<< facelist[i][1];
    for(i=0; i<Nfaces; i++)
        output <<"\n"<< Lcell[i]    <<" "<< Rcell[i]   ;

    return;
}

