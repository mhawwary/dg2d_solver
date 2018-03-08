#ifndef MYSPARSE_MATRIX_H
#define MYSPARSE_MATRIX_H

#include "general_tools.h"

/*A simple sparse matrix object
 !* Mohammad Alhawwary 03052018
 */
struct MYSPARSE_MATRIX{

    std::vector<double> val;
    std::vector<int> *ix=nullptr;
    int nnz_; // total number of non_zero elements
    std::vector<int> nnz_row; // total number of non_zero elements in each row
    int row_maxsize;
    int col_maxsize;
    int dim_;  // dimension of the matrix

    MYSPARSE_MATRIX(const int& odim_){
        nnz_=0;
        row_maxsize=0;
        col_maxsize=0;
        dim_=odim_;
        ix = new std::vector<int>[dim_];
    }

    MYSPARSE_MATRIX(){
        nnz_=0;
        dim_=2;
        row_maxsize=0;
        col_maxsize=0;
        ix = new std::vector<int>[dim_];
    }

    ~MYSPARSE_MATRIX(){
        val.clear();
        for(int i=0; i<dim_; i++)
            ix[i].clear();
        nnz_row.clear();
        emptyarray(ix);
    }

    void init(const int& odim_){
        nnz_=0;
        row_maxsize=0;
        col_maxsize=0;
        dim_=odim_;
        ix = new std::vector<int>[dim_];
        return;
    }

    void set_size(void){
        val.shrink_to_fit();
        for(int i=0; i<dim_; i++)
            ix[i].shrink_to_fit();
        nnz_=val.size();
        nnz_row.shrink_to_fit();
        row_maxsize=nnz_row.size();
        return;
    }

    void set_value(const double &value_, const int &ir, const int &ic){
        val.push_back(value_);
        ix[0].push_back(ir);
        ix[1].push_back(ic);
        return;
    }

    void set_nnz_row(const int& o_nnz_row_index_){
        nnz_row.push_back(o_nnz_row_index_);
        return;
    }

    void set_value(const double &value_, const int &ir, const int &ic
                                    , const int &im){

        val.push_back(value_);
        ix[0].push_back(ir);
        ix[1].push_back(ic);
        ix[2].push_back(im);

        return;
    }

    virtual MYSPARSE_MATRIX& operator =(const MYSPARSE_MATRIX& R_obj){
        val = R_obj.val;
        dim_=R_obj.dim_;
        nnz_=R_obj.nnz_;
        nnz_row=R_obj.nnz_row;
        row_maxsize=R_obj.row_maxsize;
        col_maxsize=R_obj.col_maxsize;

        for(int i=0; i<dim_; i++)
            ix[i] = R_obj.ix[i];
        return *this;

    }

};


#endif

