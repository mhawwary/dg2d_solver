#pragma once

#include"general_tools.h"
#include"MeshData.h"

double L1norm(const int& Nelem_, const int& Ndof_per_elem, double** quantity, double* Vol);

double L1norm_perdof(const int& Nelem_, double* quantity);

double L2norm(const int& Nelem_, const int& Ndof_per_elem, double** quantity, double* Vol);

double L2norm_perdof(const int& Nelem_, double* quantity);

void dump_field_data(double *Qv, const int& oiter
                     , std::string& dump_dir_, MeshData*& grid_data_);
void dump_exactsol_field_data(double *Qv, const int& oiter
                     , std::string& dump_dir_, MeshData*& grid_data_);

//void dump_errors(const double& L1_error_, const double& L2_error_);
//void dump_errors(const int& Nelem_
//                 ,const double& L1_error_, const double& L2_error_);
