#include <cmath>
#include <omp.h>
#include "matrix.h"
#include "solve.h"
#include "data.h"
#include <iostream>
#include <iomanip>
using namespace std;

void genu(dynprob &prob, matrix &randus, data &dat){
  double mu = prob.mu,
  sigma     = prob.sigma,
  rho       = prob.rho,
  alpha     = prob.alpha,
  delta     = prob.delta;
  double* A = prob.A;
  double* P = prob.P;
  int na = prob.na;
  int nsim = dat.nsim;
  int nf = dat.nf;
  double* P2;
  double* A2;
  int na2 = dat.na2;
  A2 = new double[na2];
  P2 = new double[na2 * na2];
  tauchen(mu, rho, sigma, sigma * prob.tauch_sds, na2, A2, P2);
  //cout << dat.resi << endl;
  #pragma omp parallel for
  for(int f = 0; f < nf; f++){
    int zi  = (int) floor(na2 * randus(nf*nsim + 3*f));
    int zil = zi;
    for(int t = 0; t < nsim; t++){
      zi  = 0;
      

      double psum = 0;
      for(int i = 0; i < na2; i++){
        psum += P2[zil * na2 + i];
        if(psum < randus(f * nsim + t)){
          zi++;
        }
      }
      dat.zs(f * nsim + t, dat.zfeed) = A2[zi];
      
      zil = zi;
      if(zi > na2-1){
        zi = na2-1;
      }
      else if(zi < 0){
        zi = 0;
      }
      dat.is(f * nsim + t, dat.zfeed) = zi;
    }
  }
  #pragma omp barrier
  
  #pragma omp master
  {
      dat.zfeed = dat.zfeed + 1;
  }
   
}

void sim(dynprob &prob, matrix &randus, data &dat){

 //cout << "here " << endl;
  double* kp        = prob.kp;
  double* bp        = prob.bp;
  double* dp        = prob.dp;
  double* A         = prob.A;
  double* P         = prob.P;
  double* k         = prob.k;
  double* b         = prob.b;
  double* vold      = prob.vold;
  double* mktbp     = prob.mktbp;
  //cout << "here.5" << endl;
  double   alpha     = prob.alpha,
  beta      = prob.beta,
  delta     = prob.delta;
  //cout << "here.6" << endl;
  int nk    = prob.nk,
  nb        = prob.nb,
  na        = prob.na;
  //cout << "here.7" << endl;
    int nsim  = dat.nsim;
  //cout << "WIPING" << endl;
  dat.wipe();
  //cout << "this isn't it." << endl;
  int nbr = dat.nbr;
  //cout << "here.8" << endl;
  int i_l     = 4,
  i_p         = 0,
  k_low       = 0,
  b_low       = 0,
  i_low       = 0;
  //cout << "here.85" << endl;
  double k_l  = k[2];
  //cout << "here.85" << endl;
  double  k_p         = 0,
  k_ll        = 0,
  k_lll       = 0;
  double tau_c = prob.tau_c;
  //cout << "here.9" << endl;
  double b_l  = 0;
  //cout << "here.9" << endl;
  double b_p  = 0,
  b_ll        = 0,
  b_lll       = 0;

  //cout << "here.91" << endl;
  double r_p  = 0,
  r_l         = 0,
  r_lll       = 0,
  r_ll        = 0;

  double w_p  = 0,
         w_l  = 0;

  //cout << "here.92" << endl;
  double v_p  = 0,
  v_l         = 0,
  v_ll        = 0;

  //cout << "here.93" << endl;
  double in_l = 0,
  in_ll       = 0;

  //cout << "here.92" << endl;
  double c_p  = 0,
  dtau        = 0,
  c_l         = 0,
  c_ll        = 0;

  //cout << "here.93" << endl;
  double cf_p = 0,
  cf_l        = 0,
  cf_ll       = 0;

  double el_p    = 0;

  double mktb_p  = 0,
         mktb_l  = 0,
         mktb_ll = 0;

  //cout << "here1 " << endl;
  int nf = dat.nf, track = 0;

  //double i_l, i_p;
  double wb, wk;

  int zi = 2,zil;
  matrix kpt(nf*(nsim-nbr)), kt(nf*(nsim-nbr)), klt(nf*(nsim-nbr)), kllt(nf*(nsim-nbr));
  matrix ist(nf*(nsim-nbr));
  matrix rt(nf*(nsim-nbr)), rlt(nf*(nsim-nbr)), rllt(nf*(nsim-nbr));
  matrix bpt(nf*(nsim-nbr)), bt(nf*(nsim-nbr)), blt(nf*(nsim-nbr)), bllt(nf*(nsim-nbr));
  matrix it(nf*(nsim-nbr)),  ilt(nf*(nsim-nbr));
  matrix year(nf*(nsim-nbr)),  firm(nf*(nsim-nbr));
  matrix vt(nf*(nsim-nbr)), vlt(nf*(nsim-nbr));
  matrix deft(nf*(nsim-nbr));
  matrix zt(nf*(nsim-nbr)), zlt(nf*(nsim-nbr));
  matrix  ct(nf*(nsim-nbr)),  clt(nf*(nsim-nbr));

  matrix  mktbt(nf*(nsim-nbr)),  mktbpt(nf*(nsim-nbr)),  mktblt(nf*(nsim-nbr));
  
  if((dat.mul != prob.mu) || (dat.sigmal != prob.sigma) || (dat.rhol != prob.rho)){
    dat.reset();
    dat.mul    = prob.mu;
    dat.sigmal = prob.sigma;
    dat.rhol   = prob.rho;
  }
  
  if(dat.resi >= dat.nrep){
    dat.resi = 0;
  }
  //cout << "genu" << endl;
  if(dat.zfeed < dat.nrep){
    genu(prob, randus, dat);
  }
    
  //cout << "done" << endl;
  double du = 0, dd = 0, dul = 0, dull = 0, dulll = 0, ddl = 0,
  ddll = 0, ddlll = 0;
  int p_i;
  double psum;
  int defa = 0;
  double defdec = 0;
  double z;
  int z_low;
  double wz;

  #pragma omp parallel for private(k_p, k_ll, k_l, c_l, c_p, r_l, b_ll, b_l, b_p, v_l, v_p, dul, du, ddl, dd, cf_l, z, i_p, b_low, k_low, wk, wb, wz, r_p, cf_p, track, in_ll, in_l, mktb_p, r_ll, z_low, i_low)
  for(int f = 0; f < nf; f++){
    b_l = b[(int) floor(nb * randus(nsim * nf + f*3 + 1))];
    k_l = k[(int) floor(nk * randus(nsim * nf + f*3 + 2))];
    for(int t = 0; t < nsim; t++){
      //cout << t << endl;
      z   = dat.zs(f * nsim + t, dat.resi);
      i_p = dat.is(f * nsim + t, dat.resi);
      z_low = 0;
      b_low = 0;

      for(int i = 1; i < na; i++){
        if(z >= A[i]){
          z_low++;
        }
      }

      if(z_low <0)
      {
        z_low = 0;
      }

      if(z_low>na-2){
        z_low = na-2;
      }

      k_low = 0;

      for(int i = 1; i < nb; i++){
        if(b_l >= b[i]){
          b_low++;
        }
      }

      if(b_low <0)
      {
        b_low = 0;
      }

      if(b_low>nb-2){
        b_low = nb-2;
      }

      k_low = 0;

      for(int i = 1; i < nk; i++){
        if(k_l >= k[i]){
          k_low++;
        }
      }

      if(k_low <0)
      {
        k_low = 0;
      }

      if(k_low>nk-2){
        k_low = nk-2;
      }



      i_low = k_low * nb  * na + b_low * na + z_low;

      wb = (b_l - b[b_low])/(b[b_low+1]-b[b_low]);
      wk = (k_l - k[k_low])/(k[k_low+1]-k[k_low]);
      wz = (z - A[z_low])/(A[z_low+1]-A[z_low]);
      
      r_p = z * pow(k_l, alpha);
  
      

    		k_p    = prob.interpolate(kp, i_low, wb, wk, wz);

    		b_p    = prob.interpolate(bp, i_low, wb, wk, wz);


  		c_p    = prob.interpolate(dp, i_low, wb, wk, wz);

  		v_p    = prob.interpolate(vold, i_low, wb, wk, wz);

      mktb_p = prob.interpolate(mktbp, i_low, wb, wk, wz);

      cf_p = 0;
      if(cf_p < 0){
        cf_p = -b_p;
      }

      el_p = 0;
      if(c_p < 0){
        el_p = -c_p;
      }

      if(t>=nbr){
  			track++;
			zt(f*(nsim-nbr) + t-nbr) = z,      zlt(f*(nsim-nbr) + t-nbr) = wk;
  			kpt(f*(nsim-nbr) + t-nbr) = k_p, kt(f*(nsim-nbr) + t-nbr)  = k_l, klt(f*(nsim-nbr) + t-nbr) = k_ll;
        vt(f*(nsim-nbr) + t-nbr)  = v_p, vlt(f*(nsim-nbr) + t-nbr)  = v_l;
        ct(f*(nsim-nbr) + t-nbr)  = c_p, clt(f*(nsim-nbr) + t-nbr)  = c_l;

        rt(f*(nsim-nbr) + t-nbr)  = r_p, rlt(f*(nsim-nbr) + t-nbr)  = r_l, rllt(f*(nsim-nbr) + t-nbr) = r_ll;

        bpt(f*(nsim-nbr) + t-nbr) = b_p, bt(f*(nsim-nbr) + t-nbr)  = b_l, blt(f*(nsim-nbr) + t-nbr) = b_ll;
        mktbpt(f*(nsim-nbr) + t-nbr) = mktb_p, mktbt(f*(nsim-nbr) + t-nbr)  = mktb_l, mktblt(f*(nsim-nbr) + t-nbr) = mktb_ll;

        year(f*(nsim-nbr) + t-nbr) = (double) t, firm(f*(nsim-nbr) + t-nbr) = (double) f;

        it(f*(nsim-nbr) + t-nbr)  = k_p - (1-delta) * k_l, ilt(f*(nsim-nbr) + t-nbr)  = in_l;
        ist(f*(nsim-nbr) + t-nbr) = el_p;
      }

      in_ll = in_l, in_l = k_p - (1-delta) * k_l;
      k_ll = k_l, k_l=k_p;
      c_l=c_p;
      r_ll = r_l, r_l = r_p;
      b_ll = b_l, b_l=b_p;
      mktb_ll = mktb_l, mktb_l=mktb_p;
      v_l=v_p;
      cf_l = cf_p;

    }
  }
  //cout << "here5 " << endl;

  dat.add_variable("kpt", kpt);
  dat.add_variable("kt", kt);
  dat.add_variable("klt", klt);
  dat.add_variable("kllt", kllt);

  dat.add_variable("bpt", bpt, kpt);
  dat.add_variable("bt", bt, kt);
  dat.add_variable("blt", blt, klt);
  dat.add_variable("bllt", bllt, kllt);
  dat.add_variable("dt", bpt, kpt);

  dat.add_variable("rt", rt, kt);
  dat.add_variable("nwt", rt, kt);
  dat.add_variable("rlt", rlt, klt);
  dat.add_variable("rllt", rllt, kllt);
  
  dat.add_variable("it", it, kt);
  dat.add_variable("ilt", ilt, klt);

  dat.add_variable("year", year);
  dat.add_variable("firm", firm);

  dat.add_variable("vt", vt, kt);
  dat.add_variable("vlt", vlt, klt);
  dat.add_variable("zt", zt);
  dat.add_variable("zlt", zlt);

  dat.add_variable("ist", ist);

  dat.add_variable("ct",  ct, kt);
  dat.add_variable("clt", clt, klt);
  dat.add_variable("mktbpt", mktbpt);
  dat.add_variable("mktbt", mktbt);
  dat.add_variable("mktblt", mktblt);
  dat.resi = dat.resi + 1;
  //return dat;

}