#include <cmath>
#include <omp.h>
#include "matrix.h"
#include "solve.h"
#include "data.h"
#include <iostream>
#include <iomanip>
using namespace std;


// Utility for normal inverse
double r8poly_value ( int n, double a[], double x ){
  int i;
  double value;

  value = 0.0;

  for ( i = n-1; 0 <= i; i-- )
  {
    value = value * x + a[i];
  }

  return value;
}


// inverse cumulative normal of p
double normalINV( double p ){
  static double a[8] = {
    3.3871328727963666080,     1.3314166789178437745e+2,
    1.9715909503065514427e+3,  1.3731693765509461125e+4,
    4.5921953931549871457e+4,  6.7265770927008700853e+4,
    3.3430575583588128105e+4,  2.5090809287301226727e+3 };
  static double b[8] = {
    1.0,                       4.2313330701600911252e+1,
    6.8718700749205790830e+2,  5.3941960214247511077e+3,
    2.1213794301586595867e+4,  3.9307895800092710610e+4,
    2.8729085735721942674e+4,  5.2264952788528545610e+3 };
  static double c[8] = {
    1.42343711074968357734,     4.63033784615654529590,
    5.76949722146069140550,     3.64784832476320460504,
    1.27045825245236838258,     2.41780725177450611770e-1,
    2.27238449892691845833e-2,  7.74545014278341407640e-4 };
  static double const1 = 0.180625;
  static double const2 = 1.6;
  static double d[8] = {
    1.0,                        2.05319162663775882187,
    1.67638483018380384940,     6.89767334985100004550e-1,
    1.48103976427480074590e-1,  1.51986665636164571966e-2,
    5.47593808499534494600e-4,  1.05075007164441684324e-9 };
  static double e[8] = {
    6.65790464350110377720,     5.46378491116411436990,
    1.78482653991729133580,     2.96560571828504891230e-1,
    2.65321895265761230930e-2,  1.24266094738807843860e-3,
    2.71155556874348757815e-5,  2.01033439929228813265e-7 };
  static double f[8] = {
    1.0,                        5.99832206555887937690e-1,
    1.36929880922735805310e-1,  1.48753612908506148525e-2,
    7.86869131145613259100e-4,  1.84631831751005468180e-5,
    1.42151175831644588870e-7,  2.04426310338993978564e-15 };
  double q;
  double r;
  static double split1 = 0.425;
  static double split2 = 5.0;
  double value;

  q = p - 0.5;

  if ( abs(q) <= split1 )
  {
    r = const1 - q * q;
    value = q * r8poly_value ( 8, a, r ) / r8poly_value ( 8, b, r );
  }
  else
  {
    if ( q < 0.0 )
    {
      r = p;
    }
    else
    {
      r = 1.0 - p;
    }

    if ( r <= 0.0 )
    {
      value = -1.0;
    }
	else{
		r = sqrt(-log(r));

		if (r <= split2)
		{
			r = r - const2;
			value = r8poly_value(8, c, r) / r8poly_value(8, d, r);
		}
		else
		{
			r = r - split2;
			value = r8poly_value(8, e, r) / r8poly_value(8, f, r);
		}

		if (q < 0.0)
		{
			value = -value;
		}
	}
  }

  return value;
}

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
  srand(dat.seedi*1000000 + 1 + 10000 * dat.resi);   
  //cout << dat.resi << endl;
  //cout << dat.rs(0, 0) << endl;
  #pragma omp parallel for
  for(int f = 0; f < nf; f++){
//    srand(dat.seedi*1000000 + 3*f + 1 + 10000 * dat.resi);
    int zi  = (int) floor(na * randus(nf*nsim + f*3));
    //int zi = na / 2;
    int zil = zi;
    for(int t = 0; t < nsim; t++){
      zi  = 0;
      

      double psum = 0;
      for(int i = 0; i < na; i++){
        psum += P[zil * na + i];
        if(psum < randus(f * nsim + t)){
          zi++;
        }
      }
      dat.zs(f * nsim + t, dat.zfeed) = A[zi];
      
      zil = zi;
      if(zi > na-2){
        zi = na-2;
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

// Simulate model prob given random distribution randus and store in data class data
void sim(dynprob &prob, matrix &randus, data &dat){


  double* kp = prob.kp;
  double* bp        = prob.bp;
  double* cp        = prob.cp;
  double* A         = prob.A;
  double* k         = prob.k;
  double* b         = prob.b;
  double* vold      = prob.vold;
  double* dp        = prob.dp;

  double mu = prob.mu,
  sigma     = prob.sigma,
  rho       = prob.rho,
  alpha     = prob.alpha,
  delta     = prob.delta,
  tau_c     = prob.tau_c,
  r_f         = (1 / prob.beta) - 1;

  int nk    = prob.nk,
  nb        = prob.nb,
  na        = prob.na;

  int nsim  = dat.nsim;
  dat.wipe();
  int nbr = dat.nbr;
  int i_l     = 4,
  i_p         = 0,
  k_low       = 0,
  b_low       = 0;
    
  double k_l  = prob.k[nk / 2],
  k_p         = 0,
  k_ll        = 0,
  k_lll       = 0,
  k_llll      = 0,
  k_lllll     = 0;

  double b_l  = prob.b[nb / 2],
  b_p         = 0,
  b_ll        = 0,
  b_lll       = 0,
  b_llll      = 0,
  b_lllll     = 0;

  double r_p  = 0,
  r_l         = 0,
  r_ll        = 0,
  r_lll       = 0,
  r_llll      = 0,
  r_lllll     = 0;

  double d_p  = 0,
  d_l         = 0,
  d_ll        = 0,
  d_lll       = 0,
  d_llll      = 0,
  d_lllll     = 0;

  double v_p  = 0,
  v_l         = 0,
  v_ll        = 0,
  v_lll       = 0,
  v_llll      = 0;

  double in_l = 0,
  in_ll       = 0,
  in_lll      = 0,
  in_llll     = 0;

  double c_p  = 0,
  dtau        = 0,
  c_l         = 0,
  c_ll        = 0,
  c_lll       = 0,
  c_llll      = 0,
  c_lllll     = 0;

  double cf_p = 0,
  cf_l        = 0,
  cf_ll       = 0,
  cf_lll      = 0,
  cf_llll     = 0;

  if((dat.mul != mu) || (dat.sigmal != sigma) || (dat.rhol != rho)){
    dat.reset();
    dat.mul    = mu;
    dat.sigmal = sigma;
    dat.rhol   = rho;
  }
  
  if(dat.resi >= dat.nrep){
    dat.resi = 0;
  }
  
  if(dat.zfeed < dat.nrep){
    genu(prob, randus, dat);
  }
  
  int nf = dat.nf, track = 0, f = 0;
  //double i_l, i_p;
  double wb, wk, wz, z, zt, zs, wtau, wtaus;
  matrix kpt(nf*(nsim-nbr)), kt(nf*(nsim-nbr)), klt(nf*(nsim-nbr));
  matrix wbt(nf*(nsim-nbr)), wat(nf*(nsim-nbr));
  matrix rt(nf*(nsim-nbr)), rlt(nf*(nsim-nbr));
  matrix vt(nf*(nsim-nbr)), vlt(nf*(nsim-nbr));
  matrix bpt(nf*(nsim-nbr)), bt(nf*(nsim-nbr)), blt(nf*(nsim-nbr));
  matrix  it(nf*(nsim-nbr)),  ilt(nf*(nsim-nbr));
  matrix  et(nf*(nsim-nbr));
  matrix  year(nf*(nsim-nbr)),  firm(nf*(nsim-nbr));
  matrix idt(nf*(nsim-nbr));
  matrix  ct(nf*(nsim-nbr)),  clt(nf*(nsim-nbr)),  cllt(nf*(nsim-nbr));
  matrix  ddt(nf*(nsim-nbr)),  ddlt(nf*(nsim-nbr)),  ddllt(nf*(nsim-nbr));
  matrix  dt(nf*(nsim-nbr)), dlt(nf*(nsim-nbr));
  matrix nwt(nf*(nsim-nbr));
	double du = 0, dd = 0, dul = 0, dull = 0, dulll = 0, ddl = 0,
  ddll = 0, ddlll = 0;
  int p_i;
  #pragma omp parallel for private(k_p, k_ll, k_l, c_l, c_p, r_l, b_ll, b_l, b_p, v_l, v_p, dul, du, ddl, dd, cf_l, z, i_p, b_low, wb, wz, r_p, cf_p, track, in_ll, in_l) 
  for(f = 0; f < nf; f++){
    //srand(dat.seedi + 3*f + 1 + 10000 * dat.resi);
    b_l = b[(int) floor(nb * randus(nf*nsim + f*3+1))];
    k_l = k[(int) floor(nk * randus(nf*nsim + f*3+2))];
    for(int t = 0; t < nsim; t++){
      z   = dat.zs(f * nsim + t, dat.resi);
      i_p = dat.is(f * nsim + t, dat.resi);
      
      b_low = 0;

      for(int i = 1; i < nb; i++){
        if(b_l > b[i]){
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
      
      wb = (b_l - b[b_low])/(b[b_low+1]-b[b_low]);
      //wz = (z - A[i_p])/(A[i_p + 1] - A[i_p]);
      //cout << wz << " " << i_p << " " << z << " " << A[i_p] << endl;
      
      if(z > A[i_p]){
        i_p = i_p + 1;
      }
      
  		k_p = prob.interpolate(kp, b_low, wb, i_p);

  		b_p = prob.interpolate(bp, b_low, wb, i_p);

  		c_p = prob.interpolate(cp, b_low, wb, i_p);

  		v_p = prob.interpolate(vold, b_low, wb, i_p);

      cf_p = 0;
      if(b_p < 0){
        cf_p = -b_p;
      }

      r_p = z;

      if(t>=nbr){
  			track++;

  			kpt(f*(nsim-nbr) + t-nbr) = k_p, kt(f*(nsim-nbr) + t-nbr)  = k_l, klt(f*(nsim-nbr) + t-nbr) = k_ll;

        vt(f*(nsim-nbr) + t-nbr)  = v_p, vlt(f*(nsim-nbr) + t-nbr)  = v_l;
        ct(f*(nsim-nbr) + t-nbr)  = c_p, clt(f*(nsim-nbr) + t-nbr)  = c_l;

        rt(f*(nsim-nbr) + t-nbr)  = r_p, rlt(f*(nsim-nbr) + t-nbr)  = r_l;

        bpt(f*(nsim-nbr) + t-nbr) = b_p, bt(f*(nsim-nbr) + t-nbr)  = b_l, blt(f*(nsim-nbr) + t-nbr) = b_ll;

        year(f*(nsim-nbr) + t-nbr) = (double) t, firm(f*(nsim-nbr) + t-nbr) = (double) f;

        it(f*(nsim-nbr) + t-nbr)  = k_p, ilt(f*(nsim-nbr) + t-nbr)  = in_l;
        dt(f*(nsim-nbr) + t-nbr)  = b_p;
        
        nwt(f*(nsim - nbr) + t-nbr) = (1 - tau_c) * r_p + (1 - delta) - b_l;
      }

      in_ll = in_l, in_l = k_p;
      k_ll = k_l, k_l=k_p;
      c_l=c_p;
      r_l = r_p;
      b_ll = b_l, b_l=b_p;
      v_l=v_p;
      dul = du;
      ddl = dd;
      cf_l = cf_p;
   
    }
  }

  #pragma omp barrier

  dat.add_variable("kpt", kpt);
  dat.add_variable("kt", kt);
  dat.add_variable("klt", klt);

  dat.add_variable("bpt", bpt);
  dat.add_variable("bt", bt);
  dat.add_variable("blt", blt);

  dat.add_variable("kt", kt);
  dat.add_variable("klt", klt);

  dat.add_variable("it", it);
  dat.add_variable("ilt", ilt);

  dat.add_variable("year", year);
  dat.add_variable("firm", firm);

  dat.add_variable("vt", vt);
  dat.add_variable("vlt", vlt);

  dat.add_variable("rt", rt);
  dat.add_variable("dt", dt);
  dat.add_variable("rlt", rlt);

  dat.add_variable("ct",  ct);
  dat.add_variable("clt", clt);
  dat.add_variable("nwt", nwt);
  dat.resi = dat.resi + 1;

  //return dat;
}
