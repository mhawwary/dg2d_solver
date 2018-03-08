#ifndef QUADRATURE_2D_H
#define QUADRATURE_2D_H

#include"general_tools.h"


/*
 * Two-dimensional gauss quadrature object
 * It works with Quads only for now
 * Mohammad Alhawwary, 03052018
 */

struct  GaussQuad2D{

public:

    unsigned int Nq=1;  // order of quadrature which is accurate for 2(Nq-1)
    unsigned int Nq_1d=1;
    std::string elem_type="Quad";  // element type ( only supports quads for now )

    double **Gaus_pts=nullptr;  // Gauss  xi/eta Points
    double *Gaus_wts=nullptr;  // Gauss weights

protected:
    double *pts_1d=nullptr;  // 1d Gauss pts
    double *wts_1d=nullptr;  // 1d Gauss weights

public:
    ~GaussQuad2D(void){
        emptyarray(2,Gaus_pts);
        emptyarray(Gaus_wts);
        emptyarray(wts_1d);
        emptyarray(pts_1d);
    }

    void setup_quadrature(const int &Nq_){

        Nq_1d=Nq_;
        Nq = Nq_1d*Nq_1d;
        wts_1d = new double[Nq_1d];
        pts_1d = new double[Nq_1d];
        Gaus_pts = new double*[2];
        Gaus_pts[0] = new double[Nq];
        Gaus_pts[1] = new double[Nq];
        Gaus_wts = new double[Nq];
        define_gauss_quadrature();

        return;

    }

protected:
    void define_gauss_quadrature(){

        if(Nq_1d==1){
            pts_1d[0] = 0.0 ;
            wts_1d[0] = 2.0 ;

        }else if(Nq_1d==2){
            pts_1d[0] = -sqrt(1.0/3) ;
            pts_1d[1] =  sqrt(1.0/3) ;

            wts_1d[0] =  1.0 ;
            wts_1d[1] =  1.0 ;

        }else if(Nq_1d==3){
            pts_1d[0] = -sqrt(3.0/5) ;
            pts_1d[1] = 0.0  ;
            pts_1d[2] = sqrt(3.0/5) ;

            wts_1d[0] =  5.0/9 ;
            wts_1d[1] =  8.0/9 ;
            wts_1d[2] =  5.0/9 ;

        }else if(Nq_1d==4){
            pts_1d[0] =  -sqrt(3.0/7+2.0/7*sqrt(6.0/5)) ;
            pts_1d[1] =  -sqrt(3.0/7-2.0/7*sqrt(6.0/5)) ;
            pts_1d[2] =   sqrt(3.0/7-2.0/7*sqrt(6.0/5)) ;
            pts_1d[3] =   sqrt(3.0/7+2.0/7*sqrt(6.0/5)) ;

            wts_1d[0] =  (18-sqrt(30))/36.0  ;
            wts_1d[1] =  (18+sqrt(30))/36.0  ;
            wts_1d[2] =  (18+sqrt(30))/36.0  ;
            wts_1d[3] =  (18-sqrt(30))/36.0  ;

        }else if(Nq_1d==5){
            pts_1d[0] =  -(1.0/3)*sqrt(5+2*sqrt(10.0/7)) ;
            pts_1d[1] =  -(1.0/3)*sqrt(5-2*sqrt(10.0/7)) ;
            pts_1d[2] =             0              ;
            pts_1d[3] =   (1.0/3)*sqrt(5-2*sqrt(10.0/7)) ;
            pts_1d[4] =   (1.0/3)*sqrt(5+2*sqrt(10.0/7)) ;

            wts_1d[0] =  (322-13*sqrt(70))/900.0 ;
            wts_1d[1] =  (322+13*sqrt(70))/900.0 ;
            wts_1d[2] =  128.0/225 ;
            wts_1d[3] =  (322+13*sqrt(70))/900.0 ;
            wts_1d[4] =  (322-13*sqrt(70))/900.0;

        }else if(Nq_1d==6){
            pts_1d[0] =  -9.3246951420315202781230155449399e-01L;
            pts_1d[1] =  -6.6120938646626451366139959501991e-01L;
            pts_1d[2] =  -2.3861918608319690863050172168071e-01L;
            pts_1d[3] =  -pts_1d[2];
            pts_1d[4] =  -pts_1d[1];
            pts_1d[5] =  -pts_1d[0];

            wts_1d[0] =  1.7132449237917034504029614217273e-01L;
            wts_1d[1] =  3.6076157304813860756983351383772e-01L;
            wts_1d[2] =  4.6791393457269104738987034398955e-01L;
            wts_1d[3] =  wts_1d[2];
            wts_1d[4] =  wts_1d[1];
            wts_1d[5] =  wts_1d[0];

        }else if(Nq_1d==7){
            pts_1d[0] =  -9.4910791234275852452618968404785e-01L;
            pts_1d[1] =  -7.4153118559939443986386477328079e-01L;
            pts_1d[2] =  -4.0584515137739716690660641207696e-01L;
            pts_1d[3] =   0.0;
            pts_1d[4] =  -pts_1d[2];
            pts_1d[5] =  -pts_1d[1];
            pts_1d[6] =  -pts_1d[0];

            wts_1d[0] =  1.2948496616886969327061143267908e-01L;
            wts_1d[1] =  2.7970539148927666790146777142378e-01L;
            wts_1d[2] =  3.8183005050511894495036977548898e-01L;
            wts_1d[3] =  4.1795918367346938775510204081633e-01L;
            wts_1d[4] =  wts_1d[2];
            wts_1d[5] =  wts_1d[1];
            wts_1d[6] =  wts_1d[0];

        }else if(Nq_1d==8){
            pts_1d[0] =  -9.6028985649753623168356086856947e-01L;
            pts_1d[1] =  -7.9666647741362673959155393647583e-01L;
            pts_1d[2] =  -5.2553240991632898581773904918925e-01L;
            pts_1d[3] =  -1.8343464249564980493947614236018e-01L;
            pts_1d[4] =  -pts_1d[3];
            pts_1d[5] =  -pts_1d[2];
            pts_1d[6] =  -pts_1d[1];
            pts_1d[7] =  -pts_1d[0];

            wts_1d[0] =  1.0122853629037625915253135430996e-01L;
            wts_1d[1] =  2.2238103445337447054435599442624e-01L;
            wts_1d[2] =  3.1370664587788728733796220198660e-01L;
            wts_1d[3] =  3.6268378337836198296515044927720e-01L;
            wts_1d[4] =  wts_1d[3];
            wts_1d[5] =  wts_1d[2];
            wts_1d[6] =  wts_1d[1];
            wts_1d[7] =  wts_1d[0];

        }else{

            std::cout<< "\n Gauss Quadrature with order of:  "
                <<Nq_1d<<"  is not implemented in the quadrature rule file\n";
        }

        setup_2d_quadratures();

        return;
    }

    void setup_2d_quadratures(){

        int i,j;

        for(i=0; i<Nq_1d; i++)
            for(j=0; j<Nq_1d; j++){
                Gaus_pts[0][i*Nq_1d+j]=pts_1d[j];
                Gaus_pts[1][i*Nq_1d+j]=pts_1d[i];
                Gaus_wts[i*Nq_1d+j] = wts_1d[i]*wts_1d[j];
            }

        emptyarray(wts_1d);
        emptyarray(pts_1d);

        return;
    }
};

#endif

