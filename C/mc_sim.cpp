#include "matrix.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include "epf.h"
#include "opt.h"
#include "regress.h"
#include <sstream>
#include <omp.h>
#include <sstream>
#include <limits>

#define _USE_MATH_DEFINES
typedef std::numeric_limits< double > dbl;

int main(int argc, char **argv){


    int thr_num  = 14;
    int n        = 7;
    int order    = 2;
    int ntrials = 1000;

    matrix v(n);
    string type = "quad_toni_rr3";
    string prefix = "../Results/MC/trial_data_rr_";
    string file = "../Results/MC/trial_data_moments.csv";

    matrix lengths = readcsv("../Data/lengths.csv",2,1);

    int nfirms  = lengths(1);
    int   nobs  = lengths(0);
    int nobsper = nobs / nfirms;
    nfirms = 75000/nobsper;
    // Make the matrix of random draws
    matrix randus(nobsper * nfirms + 50 * nfirms, 10);
    int seedn = 1158;
    srand(seedn);
    randus.randu();

    // Storeage for objective functions
      double f, f_0;
    // Set up data and dynamic problem
    data dat1 = data(50, 50 + nobsper, nfirms, 60, 1);
    dynprob res1 = dynprob(24, 64, 61, 15);

    // Initialize interest rate
    res1.pars[0]=1/(1+(1+0.1) *0.02);

    // Use CRS
    res1.pars[1]=1;

    // Initialize taxes on corporate profits
    res1.pars[14] = 0.1;

    // Set up the order of estimation
    res1.porder[0] = 2;
    res1.porder[1] = 4;
    res1.porder[2] = 9;
    res1.porder[3] = 11;
    res1.porder[4] = 21;
    res1.porder[5] = 22;
    res1.porder[6] = 23;

    // Read in data
    cout<< "FILE: " << "../Results/Estim/par_" + type + ".csv" << endl << endl;
    matrix par = readcsv("../Results/Estim/par_" + type + ".csv",n,1);
    cout << "Current parameter: " << endl;
    par.print();

    ostringstream append;
    ofstream simfile, momfile, momfile_kt, momfile_hp, quadfile, quadfile_trans, quadfile_wvar, quadfile_wvarex, momfile_cv, cubefile;
    ofstream vmomfile, vmomfile_kt, vmomfile_hp, vquadfile, vquadfile_wvar, vquadfile_wvarex, vmomfile_cv, vcubefile, vquadfile_trans;
    ofstream vmomfile2, vmomfile_cv2, vcubefile2, vmomfile_kt2, vmomfile_hp2, vquadfile2, vquadfile_wvar2, vquadfile_wvarex2, vquadfile_trans2;
    ofstream vmomquadfile, vmomcubefile, vmomquadfile_cv, vmomquadfile_kt, vmomquadfile_hp, vmomquadfile_wvar, vmomquadfile_kt_wvar, vmomquadfile_hp_wvar, vmomquadfile_wvarex, vmomquadfile_kt_wvarex, vmomquadfile_hp_wvarex, vmomquadfile_trans;

    simfile.open(file.c_str(), ios::app);
    string resr;
    string momname = prefix + "mom.csv";
    string momname_kt = prefix + "mom_kt.csv";
    string momname_hp = prefix + "mom_hp.csv";
    string quadname = prefix + "quad.csv";
    string quadname_wvar = prefix + "quad_wvar.csv";
    string quadname_trans = prefix + "quad_trans.csv";
    string quadname_wvarex = prefix + "quad_wvarex.csv";
    string momname_cv = prefix + "mom_cv.csv";
    string cubename = prefix + "cube.csv";


    string vmomname = prefix + "vmom.csv";
    string vmomname_kt = prefix + "vmom_kt.csv";
    string vmomname_hp = prefix + "vmom_hp.csv";
    string vquadname = prefix + "vquad.csv";
    string vquadname_wvar = prefix + "vquad_wvar.csv";
    string vquadname_wvarex = prefix + "vquad_wvar.csv";
    string vcubename = prefix + "vcube.csv";
    string vmomname_cv = prefix + "vmom_cv.csv";
    string vquadname_trans = prefix + "vquad_trans.csv";

    string vmomname2 = prefix + "vmom_unclust.csv";
    string vmomname_kt2 = prefix + "vmom_kt_unclust.csv";
    string vmomname_hp2 = prefix + "vmom_hp_unclust.csv";
    string vquadname2 = prefix + "vquad_unclust.csv";
    string vquadname_wvar2 = prefix + "vquad_wvar_unclust.csv";
    string vquadname_wvarex2 = prefix + "vquad_wvarex_unclust.csv";
    string vcubename2 = prefix + "vcube_unclust.csv";
    string vmomname_cv2 = prefix + "vmom_cv_unclust.csv";
    string vquadname_trans2 = prefix + "vquad_trans_unclust.csv";

    string vmomnamequad = prefix + "vmom_quad.csv";
    string vmomname_ktquad = prefix + "vmom_kt_quad.csv";
    string vmomname_hpquad = prefix + "vmom_hp_quad.csv";
    string vmomnamequad_wvar = prefix + "vmom_quad_wvar.csv";
    string vmomname_ktquad_wvar = prefix + "vmom_kt_quad_wvar.csv";
    string vmomname_hpquad_wvar = prefix + "vmom_hp_quad_wvar.csv";
    string vmomnamequad_wvarex = prefix + "vmom_quad_wvarex.csv";
    string vmomname_ktquad_wvarex = prefix + "vmom_kt_quad_wvarex.csv";
    string vmomname_hpquad_wvarex = prefix + "vmom_hp_quad_wvarex.csv";
    string vmomnamecube = prefix + "vmom_cube.csv";
    string vmomnamequad_trans = prefix + "vmom_quad_trans.csv";
    string vmomname_cvquad = prefix + "vmom_cv_quad.csv";


    momfile.open(momname.c_str(), ios::app);
    momfile_kt.open(momname_kt.c_str(), ios::app);
    momfile_hp.open(momname_hp.c_str(), ios::app);
    quadfile.open(quadname.c_str(), ios::app);
    quadfile_trans.open(quadname_trans.c_str(), ios::app);
    quadfile_wvar.open(quadname_wvar.c_str(), ios::app);
    quadfile_wvarex.open(quadname_wvarex.c_str(), ios::app);
    momfile_cv.open(momname_cv.c_str(), ios::app);
    cubefile.open(cubename.c_str(), ios::app);


    vmomfile.open(vmomname.c_str(), ios::app);
    vmomfile_kt.open(vmomname_kt.c_str(), ios::app);
    vmomfile_hp.open(vmomname_hp.c_str(), ios::app);
    vquadfile.open(vquadname.c_str(), ios::app);
    vquadfile_wvar.open(vquadname_wvar.c_str(), ios::app);
    vquadfile_trans.open(vquadname_trans.c_str(), ios::app);
    vquadfile_wvarex.open(vquadname_wvarex.c_str(), ios::app);
    vmomfile_cv.open(vmomname_cv.c_str(), ios::app);
    vcubefile.open(vcubename.c_str(), ios::app);


    vmomfile2.open(vmomname2.c_str(), ios::app);
    vmomfile_kt2.open(vmomname_kt2.c_str(), ios::app);
    vmomfile_hp2.open(vmomname_hp2.c_str(), ios::app);
    vquadfile2.open(vquadname2.c_str(), ios::app);
    vquadfile_wvar2.open(vquadname_wvar2.c_str(), ios::app);
    vquadfile_trans2.open(vquadname_trans2.c_str(), ios::app);
    vquadfile_wvarex2.open(vquadname_wvarex2.c_str(), ios::app);
    vmomfile_cv2.open(vmomname_cv2.c_str(), ios::app);
    vcubefile2.open(vcubename2.c_str(), ios::app);


    vmomquadfile.open(vmomnamequad.c_str(), ios::app);
    vmomquadfile_kt.open(vmomname_ktquad.c_str(), ios::app);
    vmomquadfile_hp.open(vmomname_hpquad.c_str(), ios::app);
    vmomquadfile_wvar.open(vmomnamequad_wvar.c_str(), ios::app);
    vmomquadfile_kt_wvar.open(vmomname_ktquad_wvar.c_str(), ios::app);
    vmomquadfile_hp_wvar.open(vmomname_hpquad_wvar.c_str(), ios::app);
    vmomquadfile_wvarex.open(vmomnamequad_wvarex.c_str(), ios::app);
    vmomquadfile_trans.open(vmomnamequad_trans.c_str(), ios::app);
    vmomquadfile_kt_wvarex.open(vmomname_ktquad_wvarex.c_str(), ios::app);
    vmomquadfile_hp_wvarex.open(vmomname_hpquad_wvarex.c_str(), ios::app);
    vmomquadfile_cv.open(vmomname_cvquad.c_str(), ios::app);
    vmomcubefile.open(vmomnamecube.c_str(), ios::app);

    simfile.precision(dbl::digits10), momfile.precision(dbl::digits10), momfile_cv.precision(dbl::digits10), cubefile.precision(dbl::digits10), momfile_kt.precision(dbl::digits10), momfile_hp.precision(dbl::digits10), quadfile.precision(dbl::digits10), quadfile_trans.precision(dbl::digits10), quadfile_wvar.precision(dbl::digits10), quadfile_wvarex.precision(dbl::digits10);
    vmomfile.precision(dbl::digits10), vmomfile_cv.precision(dbl::digits10), vcubefile.precision(dbl::digits10), vmomfile_kt.precision(dbl::digits10), vmomfile_hp.precision(dbl::digits10), vquadfile.precision(dbl::digits10), vquadfile_wvar.precision(dbl::digits10), vquadfile_trans.precision(dbl::digits10), vquadfile_wvarex.precision(dbl::digits10);
    vmomfile2.precision(dbl::digits10), vcubefile2.precision(dbl::digits10), vmomfile_cv2.precision(dbl::digits10), vmomfile_kt2.precision(dbl::digits10), vmomfile_hp2.precision(dbl::digits10), vquadfile2.precision(dbl::digits10), vquadfile_wvar2.precision(dbl::digits10), vquadfile_trans2.precision(dbl::digits10), vquadfile_wvarex2.precision(dbl::digits10);
    vmomquadfile.precision(dbl::digits10), vmomcubefile.precision(dbl::digits10), vmomquadfile_cv.precision(dbl::digits10), vmomquadfile_kt.precision(dbl::digits10), vmomquadfile_hp.precision(dbl::digits10), vmomquadfile_wvar.precision(dbl::digits10), vmomquadfile_trans.precision(dbl::digits10), vmomquadfile_kt_wvar.precision(dbl::digits10), vmomquadfile_hp_wvar.precision(dbl::digits10), vmomquadfile_wvarex.precision(dbl::digits10), vmomquadfile_kt_wvarex.precision(dbl::digits10), vmomquadfile_hp_wvarex.precision(dbl::digits10);

    parameter_update(par, res1);
    string prefix_test = "../Results/Test/";
    res1.vfi(1e-5,1e+16,3000,0,0,0,thr_num);
    res1.policy();

    int start;

    matrix report;

    for(int t = 0; t < ntrials; t++){

        seedn = 18000*t + 18;
        srand(seedn);
        randus.randu();

        report = report_all_mc(res1, dat1.nsim, dat1.nbr, dat1, randus);
        
        simfile << seedn;
        for(int i = 0; i < report.nrows(); i++){
            simfile << ", " << report(i);
        }
        simfile << endl;

        start = 0;

        momfile << report(0);
        for(int i = 1; i < 8; i++){
           momfile << ", " << report(i);
        }
        momfile << endl;

        start += 8;

        momfile_kt << report(0 + start);
        for(int i = 1; i < 18; i++){
           momfile_kt << ", " << report(i + start);
        }
        momfile_kt << endl;

        start += 18;

        momfile_hp << report(0 + start);
        for(int i = 1; i < 18; i++){
           momfile_hp << ", " << report(i + start);
        }
        momfile_hp << endl;

        start += 18;

        quadfile << report(0 + start);
        for(int i = 1; i < 18; i++){
           quadfile << ", " << report(i + start);
        }
        quadfile << endl;

        start += 18;

        quadfile_wvar << report(0 + start);
        for(int i = 1; i < 21; i++){
           quadfile_wvar << ", " << report(i + start);
        }
        quadfile_wvar << endl;

        start += 21;

        quadfile_wvar << report(0 + start);
        for(int i = 1; i < 21; i++){
           quadfile_wvar << ", " << report(i + start);
        }
        quadfile_wvar << endl;

        start += 21;

        cubefile << report(0 + start);
        for(int i = 1; i < 30; i++){
           cubefile << ", " << report(i + start);
        }
        cubefile << endl;

        start += 30;

        momfile_cv << report(0 + start);
        for(int i = 1; i < 11; i++){
           momfile_cv << ", " << report(i + start);
        }
        momfile_cv << endl;

        start += 11;

        quadfile_trans << report(0 + start);
        for(int i = 1; i < 18; i++){
           quadfile_trans << ", " << report(i + start);
        }
        quadfile_trans << endl;

        start += 18;

        vmomfile << report(0 + start);
        for(int i = 1; i < 8 * 8; i++){
           vmomfile << ", " << report(i + start);
        }
        vmomfile << endl;

        start += 8 * 8;

        vmomfile_kt << report(0 + start);
        for(int i = 1; i < 18 * 18; i++){
           vmomfile_kt << ", " << report(i + start);
        }
        vmomfile_kt << endl;

        start += 18 * 18;

        vmomfile_hp << report(0 + start);
        for(int i = 1; i < 18 * 18; i++){
           vmomfile_hp << ", " << report(i + start);
        }
        vmomfile_hp << endl;

        start += 18 * 18;

        vquadfile << report(0 + start);
        for(int i = 1; i < 18 * 18; i++){
           vquadfile << ", " << report(i + start);
        }
        vquadfile << endl;

        start += 18 * 18;

        vquadfile_wvar << report(0 + start);
        for(int i = 1; i < 21 * 21; i++){
           vquadfile_wvar << ", " << report(i + start);
        }
        vquadfile_wvar << endl;

        start += 21 * 21;

        vquadfile_wvarex << report(0 + start);
        for(int i = 1; i < 21 * 21; i++){
           vquadfile_wvarex << ", " << report(i + start);
        }
        vquadfile_wvarex << endl;

        start += 21 * 21;

        vcubefile << report(0 + start);
        for(int i = 1; i < 30 * 30; i++){
           vcubefile << ", " << report(i + start);
        }
        vcubefile << endl;

        start += 30 * 30;

        vmomfile_cv << report(0 + start);
        for(int i = 1; i < 11 * 11; i++){
           vmomfile_cv << ", " << report(i + start);
        }
        vmomfile_cv << endl;

        start += 11 * 11;

        vquadfile_trans << report(0 + start);
        for(int i = 1; i < 18 * 18; i++){
           vquadfile_trans << ", " << report(i + start);
        }
        vquadfile_trans << endl;

        start += 18 * 18;

        vmomfile2 << report(0 + start);
        for(int i = 1; i < 8 * 8; i++){
           vmomfile2 << ", " << report(i + start);
        }
        vmomfile2 << endl;

        start += 8 * 8;

        vmomfile_kt2 << report(0 + start);
        for(int i = 1; i < 18 * 18; i++){
           vmomfile_kt2 << ", " << report(i + start);
        }
        vmomfile_kt2 << endl;

        start += 18 * 18;

        vmomfile_hp2 << report(0 + start);
        for(int i = 1; i < 18 * 18; i++){
           vmomfile_hp2 << ", " << report(i + start);
        }
        vmomfile_hp2 << endl;

        start += 18 * 18;

        vquadfile2 << report(0 + start);
        for(int i = 1; i < 18 * 18; i++){
           vquadfile2 << ", " << report(i + start);
        }
        vquadfile2 << endl;

        start += 18 * 18;

        vquadfile_wvar2 << report(0 + start);
        for(int i = 1; i < 21 * 21; i++){
           vquadfile_wvar2 << ", " << report(i + start);
        }
        vquadfile_wvar2 << endl;

        start += 21 * 21;


        vquadfile_wvarex2 << report(0 + start);
        for(int i = 1; i < 21 * 21; i++){
           vquadfile_wvarex2 << ", " << report(i + start);
        }
        vquadfile_wvarex2 << endl;

        start += 21 * 21;


        vcubefile2 << report(0 + start);
        for(int i = 1; i < 30 * 30; i++){
           vcubefile2 << ", " << report(i + start);
        }
        vcubefile2 << endl;

        start += 30 * 30;

        vmomfile_cv2 << report(0 + start);
        for(int i = 1; i < 11 * 11; i++){
           vmomfile_cv2 << ", " << report(i + start);
        }
        vmomfile_cv2 << endl;

        start += 11 * 11;

        vquadfile_trans2 << report(0 + start);
        for(int i = 1; i < 18 * 18; i++){
           vquadfile_trans2 << ", " << report(i + start);
        }
        vquadfile_trans2 << endl;

        start += 18 * 18;

        vmomquadfile << report(0 + start);
        for(int i = 1; i < 8 * 18; i++){
           vmomquadfile << ", " << report(i + start);
        }
        vmomquadfile << endl;

        start += 8 * 18;

        vmomquadfile_kt << report(0 + start);
        for(int i = 1; i < 18 * 18; i++){
           vmomquadfile_kt << ", " << report(i + start);
        }
        vmomquadfile_kt << endl;

        start += 18 * 18;

        vmomquadfile_hp << report(0 + start);
        for(int i = 1; i < 18 * 18; i++){
           vmomquadfile_hp << ", " << report(i + start);
        }
        vmomquadfile_hp << endl;

        start += 18 * 18;

        vmomquadfile_wvar << report(0 + start);
        for(int i = 1; i < 8 * 21; i++){
           vmomquadfile_wvar << ", " << report(i + start);
        }
        vmomquadfile_wvar << endl;

        start += 8 * 21;

        vmomquadfile_kt_wvar << report(0 + start);
        for(int i = 1; i < 18 * 21; i++){
           vmomquadfile_kt_wvar << ", " << report(i + start);
        }
        vmomquadfile_kt_wvar << endl;

        start += 18 * 21;

        vmomquadfile_hp_wvar << report(0 + start);
        for(int i = 1; i < 18 * 21; i++){
           vmomquadfile_hp_wvar << ", " << report(i + start);
        }
        vmomquadfile_hp_wvar << endl;

        start += 18 * 21;

        vmomquadfile_wvarex << report(0 + start);
        for(int i = 1; i < 8 * 21; i++){
           vmomquadfile_wvarex << ", " << report(i + start);
        }
        vmomquadfile_wvarex << endl;

        start += 8 * 21;

        vmomquadfile_kt_wvarex << report(0 + start);
        for(int i = 1; i < 18 * 21; i++){
           vmomquadfile_kt_wvarex << ", " << report(i + start);
        }
        vmomquadfile_kt_wvarex << endl;

        start += 18 * 21;

        vmomquadfile_hp_wvarex << report(0 + start);
        for(int i = 1; i < 18 * 21; i++){
           vmomquadfile_hp_wvarex << ", " << report(i + start);
        }
        vmomquadfile_hp_wvarex << endl;

        start += 18 * 21;

        vmomcubefile << report(0 + start);
        for(int i = 1; i < 8 * 30; i++){
           vmomcubefile << ", " << report(i + start);
        }
        vmomcubefile << endl;

        start += 8 * 30;


        vmomquadfile_cv << report(0 + start);
        for(int i = 1; i < 11 * 18; i++){
           vmomquadfile_cv << ", " << report(i + start);
        }
        vmomquadfile_cv << endl;

        start += 11*18;

        vmomquadfile_trans << report(0 + start);
        for(int i = 1; i < 8 * 18; i++){
           vmomquadfile_trans << ", " << report(i + start);
        }
        vmomquadfile_trans << endl;

        start += 8 * 18;

        cout << "Trial " << t << ", thread" << omp_get_thread_num() << endl;
	dat1.reset();
    }

    simfile.close();
    momfile.close();
    momfile_cv.close();
    momfile_kt.close();
    momfile_hp.close();
    quadfile.close();
    cubefile.close();
    vmomfile.close();
    vmomfile_kt.close();
    vmomfile_hp.close();
    vmomfile_cv.close();
    vquadfile.close();
    vcubefile.close();
    vmomfile2.close();
    vcubefile2.close();
    vmomfile_cv2.close();
    vmomfile_kt2.close();
    vmomfile_hp2.close();
    vquadfile2.close();
    vmomquadfile.close();
    vmomcubefile.close();
    vmomquadfile_cv.close();
    vmomquadfile_kt.close();
    vmomquadfile_hp.close();
}
