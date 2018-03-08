
#include"solver_tools.h"

double L1norm(const int& Nelem_, const int& Ndof_per_elem
              , double** quantity_diff, double* Vol){

    double L1norm_=0.0;

    register int i; int j;

    for(i=0; i<Nelem_; i++)
        for(j=0; j<Ndof_per_elem; j++)
            L1norm_ += fabs(quantity_diff[i][j]*Vol[i]);

    L1norm_ = L1norm_/(Ndof_per_elem*Nelem_);

    return L1norm_;
}

double L1norm_perdof(const int& Nelem_, double* quantity_diff){

    double L1norm_=0.0;

    register int i;

    for(i=0; i<Nelem_; i++)
            L1norm_ += fabs(quantity_diff[i]);

    L1norm_ = L1norm_/Nelem_;

    return L1norm_;
}

double L2norm(const int& Nelem_, const int& Ndof_per_elem
              , double** quantity_diff, double* Vol){

    double L2norm_=0.0;

    register int i; int j;

    for(i=0; i<Nelem_; i++)
        for(j=0; j<Ndof_per_elem; j++)
            L2norm_ += pow(quantity_diff[i][j]*Vol[i],2);

    L2norm_ = sqrt(L2norm_/(Ndof_per_elem*Nelem_));

    return L2norm_;
}

double L2norm_perdof(const int& Nelem_, double* quantity_diff){

    double L2norm_=0.0;

    register int i;

    for(i=0; i<Nelem_; i++)
            L2norm_ += pow(quantity_diff[i],2);

    L2norm_ = sqrt(L2norm_/Nelem_);

    return L2norm_;
}

void dump_field_data(double *Qv, const int& oiter
                     , std::string& dump_dir_, MeshData*& grid_data){

    register int i; int j;

    char *fname=nullptr; fname=new char[100];

    sprintf(fname,"%s/field_output/contour_data_%d.plt"
            ,dump_dir_.c_str(),oiter);

    FILE*  outfile=fopen(fname,"wt");

    // Printing the header first:
    //-------------------------------------
    fprintf(outfile, "VARIABLES = \"X\",\"Y\",\"u\"\n");

    fprintf(outfile, "ZONE N=%d, E=%d, F=FEPOINT, ET=QUADRILATERAL\n"
            , grid_data->Nnodes_postproc
            , grid_data->NpostProc);

    // Printing the solution variables at the nodes:
    //------------------------------------------------
    for(i=0; i<grid_data->Nnodes_postproc; i++)
        fprintf(outfile,"%e %e %e\n",grid_data->post_Xn[i]
                ,grid_data->post_Yn[i],Qv[i]);


    fprintf(outfile,"\n");

    // printing connectivity information:
    //-------------------------------------
    int node_id;

    for(i=0; i<grid_data->NpostProc; i++)
    {
        for(j=0; j<grid_data->post_proc_elemlist[i].n_local_nodes; j++){
            node_id = grid_data->post_proc_elemlist[i].to_node[j];
            fprintf(outfile, "%d ",node_id+1);
        }

        if(grid_data->post_proc_elemlist[i].n_local_nodes==3)
            fprintf(outfile, "%d",node_id+1);
        fprintf(outfile,"\n");
    }

    fclose(outfile);
    emptyarray(fname);

    return;
}

void dump_exactsol_field_data(double *Qv, const int& oiter
                     , std::string& dump_dir_, MeshData*& grid_data){

    register int i; int j;

    char *fname=nullptr; fname=new char[100];

    sprintf(fname,"%s/field_output/contour_data_exact_%d.plt"
            ,dump_dir_.c_str(),oiter);

    FILE*  outfile=fopen(fname,"wt");

    // Printing the header first:
    //-------------------------------------
    fprintf(outfile, "VARIABLES = \"X\",\"Y\",\"u\"\n");

    fprintf(outfile, "ZONE N=%d, E=%d, F=FEPOINT, ET=QUADRILATERAL\n"
            , grid_data->Nnodes_postproc
            , grid_data->NpostProc);

    // Printing the solution variables at the nodes:
    //------------------------------------------------
    for(i=0; i<grid_data->Nnodes_postproc; i++)
        fprintf(outfile,"%e %e %e\n",grid_data->post_Xn[i]
                ,grid_data->post_Yn[i],Qv[i]);


    fprintf(outfile,"\n");

    // printing connectivity information:
    //-------------------------------------
    int node_id;

    for(i=0; i<grid_data->NpostProc; i++)
    {
        for(j=0; j<grid_data->post_proc_elemlist[i].n_local_nodes; j++){
            node_id = grid_data->post_proc_elemlist[i].to_node[j];
            fprintf(outfile, "%d ",node_id+1);
        }

        if(grid_data->post_proc_elemlist[i].n_local_nodes==3)
            fprintf(outfile, "%d",node_id+1);
        fprintf(outfile,"\n");
    }

    fclose(outfile);
    emptyarray(fname);

    return;
}

