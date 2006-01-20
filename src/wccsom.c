/* wccsom.c: for mapping patters using the WCC criterion as a distance
   measure.
   This version: Jan 19, 2006
   Author: Ron Wehrens

   Adapted from:
   *
   *  class/class.c by W. N. Venables and B. D. Ripley  Copyright (C) 1994-2002
   */

#include <R.h>

#define RANDIN  GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()

double wcc_crosscorr(double *, double *, int, double *, int);
double wcc_autocorr(double *, int, double *, int);
double wcc_corr(double *, double *, int);

void WCC_onlineSOM(double *data, double *codes, double *nhbrdist,
		   double *alphas, double *pradius, double *wghts,
		   double *dataAcors, double *Acors, double *changes,
		   Sint *pn, Sint *pp, Sint *pncodes, Sint *prlen,
		   Sint *ptrwdth)
{
  int n = *pn, p = *pp, ncodes = *pncodes, trwdth= *ptrwdth, rlen = *prlen; 
  int cd, i, j, k, l, nearest, fivers, twenties, iter;
  double dm, wccval, tmp, alpha, radius, decay;

  /* radius: fast exponential decay ending at zero; lateron we add
     1 so that at least one group of neighbours will be
     updated. Minimal value of 0.1 */
  if (*pradius > 1.1)
    decay = log(*pradius-1.0);
  else
    decay = 0.1;

  /* First calculate all autocors */
  for (cd = 0; cd < ncodes; cd++)
    Acors[cd] = wcc_autocorr(codes + cd*p, p, wghts, trwdth);
  for (i = 0; i < n; i++)
    dataAcors[i] = wcc_autocorr(data + i*p, p, wghts, trwdth);

  RANDIN;

  fivers = rlen*n / 20;
  twenties = rlen*n / 5;

  iter = 0;
  for (k = 0; k < rlen; k++) {

    alpha = alphas[1] + ((alphas[0] - alphas[1]) * (rlen - k) / rlen);
    radius = 1.0 + (*pradius-1.0) * exp(-decay * k / rlen);

    changes[k] = 0.0;

    for (l = 0; l < n; l++) {
      if (iter++ % twenties == 0)
	fprintf(stdout, "%d%%", 20 * iter / twenties);
      else if (iter % fivers == 0)
	fprintf(stdout, ".");
      fflush(stdout);

      /* pick a random data point */
      i = (int)(n * UNIF);
      
      /* find the nearest code 'near' */
      nearest = 0;
      dm = -1.0;
      for (cd=0; cd<ncodes; cd++) {
	wccval = wcc_crosscorr(data + i*p, 
			       codes + cd*p, 
			       p, wghts, trwdth)/(Acors[cd] * dataAcors[i]);
	if (wccval > dm) {
	  nearest = cd;
	  dm = wccval;
	}
      }
      
      /* update all codes within certain radius of 'nearest' */
      for (cd = 0; cd < ncodes; cd++) {
	if (nhbrdist[cd + ncodes*nearest] > radius) continue;
	for (j = 0; j < p; j++) {
	  tmp = alpha / (1.0 + nhbrdist[cd + ncodes*nearest]);
	  codes[cd*p + j] += tmp * (data[i*p + j] - codes[cd*p + j]);
	}
	
	Acors[cd] = wcc_autocorr(codes + cd*p, p, wghts, trwdth);

	if (cd == nearest) {
	  tmp = wcc_crosscorr(data + i*p,
			      codes + cd*p,
			      p, wghts, trwdth)/(Acors[cd] * dataAcors[i]);
	  changes[k] += tmp - dm;
	}

      }
    }
  }
  
  for (k = 0; k < rlen; k++)
    changes[k] /= n;

  fprintf(stdout, "\n");
  fflush(stdout);
  RANDOUT;
}


/* knn1 procedure using the wcc value as a distance */

void wccassign(double *data, double *dataAcors, 
	       double *codes, double *Acors, 
	       double *wghts, int *classif, double *bestwccs,
	       Sint *pn, Sint *pp, Sint *pncodes, Sint *ptrwdth)
{
  int n = *pn, p = *pp, ncodes = *pncodes, trwdth= *ptrwdth;
  int cd, i, nearest;
  double dm, wccval;

  for (i=0; i<n; i++) {
    nearest = -1;
    dm = 0.0;

    for (cd=0; cd<ncodes; cd++) {
      wccval = wcc_crosscorr(data + i*p, codes + cd*p, p, wghts, trwdth) / 
	(Acors[cd]*dataAcors[i]);

      if (wccval > dm) {
	nearest = cd;
	dm = wccval;
      }
    }

    classif[i] = nearest+1;
    bestwccs[i] = dm;
  }
}


double wcc_crosscorr(double *p1, double *p2, int np, 
		     double *wghts, int trwdth) {
  int i;
  double crosscov;

  crosscov = wcc_corr(p1, p2, np);
  
  for (i=1; i<trwdth; i++) {
    crosscov+=(wcc_corr(p1, p2+i, np-i)*wghts[i]);
    crosscov+=(wcc_corr(p1+i, p2, np-i)*wghts[i]);
  }
    
  return(crosscov);
}

double wcc_autocorr(double *p1, int np, double *wghts, int trwdth) {
  int i;
  double autocov;

  autocov = wcc_corr(p1, p1, np);

  for (i=1; i<trwdth; i++) 
    autocov += 2*(wcc_corr(p1, p1+i, np-i)*wghts[i]);

  return(sqrt(autocov));
}

double wcc_corr(double *p1, double *p2, int npoints)
{
  int i;
  double anum;

  anum = 0.0;
  for (i=0; i<npoints; i++)
    anum += p1[i]*p2[i];

  return(anum);
}

void wccdist(double *p1, double *p2, Sint *pnpoints, 
	     double *wghts, Sint *ptrwdth, double *WCC)
{
  int npoints = *pnpoints, trwdth=*ptrwdth;

  *WCC = wcc_crosscorr(p1, p2, npoints, wghts, trwdth);
}

void wacdist(double *p1, Sint *pnpoints, 
	     double *wghts, Sint *ptrwdth, double *ACC) 
{
  int npoints = *pnpoints, trwdth=*ptrwdth;

  *ACC = wcc_autocorr(p1, npoints, wghts, trwdth);
}


void wacdists(double *patterns, Sint *pnobj, Sint *pnpoints, 
	      double *wghts, Sint *ptrwdth, double *ACC) 
{
  int nobj = *pnobj, npoints = *pnpoints, trwdth=*ptrwdth;
  int i;

  for (i=0; i<nobj; i++)
    ACC[i] = wcc_autocorr(patterns + i*npoints, npoints, wghts, trwdth);
}


