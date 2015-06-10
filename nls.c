#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include "nls.h"

#include <k.h>

#define INDEX(i,j,rows) ((i) + (rows)*(j))

double* y_values;
/* double* x_values; not needed, keep orig q x values*/
int ncol_x, nrow_x;
double* th_values;
char* fcnName;
double* weights;
K xArgGlobal;

void fcn(int* m, int* n, double* x, double* fvec, int* iflag) {
     /* this is a C function which calls q to evaluate the model function */
     
     fflush(stdout);
     K thetaArg, ans;
     
     if(*iflag == 0) {
	  printf("current estimates:\n\t");
	  for(int ii=0; ii != *n; ++ii) {
	       printf("%6.9f ", x[ii]);
	  }
	  printf("\n\n");
     }
     
     /* create function arguments */
     char qCall[2000] = {'\0'};
     sprintf(qCall, "{[x;y] %s[;y] each x}", fcnName);
     
     thetaArg  = ktn(KF, *n);
     for(int jj=0; jj!= *n; ++jj)
	  kF(thetaArg)[jj] = x[jj];

     /* make call to q and get function values */
     ans = k(0, qCall, r1(xArgGlobal), thetaArg, (K)0);

     for(int ii=0; ii != *m; ++ii)
	  fvec[ii] = weights[ii] * (y_values[ii] - kF(ans)[ii]);

     fflush(stdout);
}

void decode_opts(K opts,
		 double* ftol, double* xtol, double* gtol, int* maxfev, 
		 double* epsfcn, double* factor,  int* nprint) {

     /* try to be very forgiving to the user */
     if( 99 == opts->t ) {
	  K optKeys = kK(opts)[0];
	  K optVals = kK(opts)[1];

	  for(int ii = 0; ii != optKeys->n; ++ii) {
	       char* optName = kS(optKeys)[ii];

	       if(0 == strcmp(optName, "nprint")) {
		    *nprint = (int)kF(optVals)[ii];
		    printf("nprint changed to %d\n", *nprint);
		    continue;
	       }

	       if(0 == strcmp(optName, "maxfev")) {
		    *maxfev = (int)kF(optVals)[ii];
		    printf("maxfev changed to %d\n", *maxfev);
		    continue;
	       }

	       if(0 == strcmp(optName, "ftol")) {
		    *ftol = kF(optVals)[ii];
		    printf("ftol changed to %e\n", *ftol);
		    continue;
	       }

	       if(0 == strcmp(optName, "xtol")) {
		    *xtol = kF(optVals)[ii];
		    printf("xtol changed to %e\n", *xtol);
		    continue;
	       }

	       if(0 == strcmp(optName, "gtol")) {
		    *gtol = kF(optVals)[ii];
		    printf("gtol changed to %e\n", *gtol);
		    continue;
	       }

	       if(0 == strcmp(optName, "epsfcn")) {
		    *epsfcn = kF(optVals)[ii];
		    printf("epsfcn changed to %e\n", *epsfcn);
		    continue;
	       }

	       if(0 == strcmp(optName, "factor")) {
		    *factor = kF(optVals)[ii];
		    printf("factor changed to %e\n", *factor);
		    continue;
	       }

	  }
     }
}

		 
K q_nls(K yArg, K fnName, K xArg, K thetaArg, K wtsArg, K opts) {
     K ans = (K)0;
     K resid, coef, dof, iinfo;
     xArgGlobal = xArg;

     int n_fcnName = fnName->n;
     fcnName = (char*)malloc((n_fcnName+1)*sizeof(char));
     strncpy((char*)fcnName,(const char*)kS(fnName), n_fcnName);
     fcnName[n_fcnName] = '\0';

     int nrow_y, npars;
     nrow_y = yArg->n;
     nrow_x = xArg->n;
     ncol_x = kK(xArg)[0]->n;
     npars = thetaArg->n;

     if(nrow_x != nrow_y)
	  krr(ss("dim mismatch between LHS and RHS"));

     y_values = (double*)malloc(nrow_y * sizeof(double));
     th_values =(double*)malloc(npars * sizeof(double));
     weights = (double*)malloc(nrow_y * sizeof(double));

     double *fvec /*residuals*/, *diag, *fjac, *qtf, *wa1, *wa2, *wa3, *wa4;
     int *ipvt;

     /* default values these are configurable on the q side */
     double ftol=1.0e-9, xtol=1.0e-9, gtol=0.0, epsfcn=0.0, factor=100.0;
     int maxfev=2000*(npars+1), mode=1, nprint=0, info=0, nfev, ldfjac=nrow_y;

     decode_opts(opts, &ftol, &xtol, &gtol, &maxfev, &epsfcn, &factor, &nprint);

     fvec = (double*)malloc(nrow_y * sizeof(double));
     diag = (double*)malloc(npars * sizeof(double));
     fjac = (double*)malloc(nrow_y * npars * sizeof(double));
     qtf = (double*)malloc(npars * sizeof(double));
     wa1 = (double*)malloc(npars * sizeof(double));
     wa2 = (double*)malloc(npars * sizeof(double));
     wa3 = (double*)malloc(npars * sizeof(double));
     wa4 = (double*)malloc(nrow_y * sizeof(double));
     ipvt = (int*)malloc(npars * sizeof(int));
     
     /* set initial */
     for(int ii = 0; ii != npars; ++ii)
	  th_values[ii] = kF(kK(thetaArg)[ii])[0];

     for(int ii = 0; ii != nrow_x; ++ii) {
	  y_values[ii] = kF(kK(yArg)[ii])[0];
	  weights[ii] = sqrt(kF(kK(wtsArg)[ii])[0]);
     }

     /* optimization call */
     lmdif_(fcn, &nrow_y, &npars, th_values, fvec, &ftol, &xtol, &gtol, &maxfev,
	    &epsfcn, diag, &mode, &factor, &nprint, &info, &nfev, fjac, &ldfjac,
	    ipvt, qtf, wa1, wa2, wa3, wa4);

     resid = ktn(KF,0);
     coef = ktn(KF,0);
     dof = ki(nrow_y - npars);
     iinfo = ki(info);

     for(int ii = 0; ii != nrow_y; ++ii) {
	  double r = fvec[ii] / weights[ii];
	  ja(&resid, &r);
     }

     for(int ii = 0; ii != npars; ++ii)
	  ja(&coef, &th_values[ii]);

     ans = ktn(0, 0);
     jk(&ans, iinfo);
     jk(&ans, coef);
     jk(&ans, resid);
     jk(&ans, dof);

     free(fcnName);
     free(y_values);
     free(th_values);
     free(weights);
     free(fvec);
     free(diag);
     free(fjac);
     free(qtf);
     free(wa1);
     free(wa2);
     free(wa3);
     free(wa4);
     free(ipvt);

     fflush(stdout);
     return(ans);
}

void __attribute((constructor)) init_function(void) {
     printf("nls minpack C library loaded\n");
}

void __attribute((destructor)) fini_function(void) {
     printf("minpack library unloaded\n");
}
