#ifndef NLS_H
#define NLS_H

/* minpack interface */
void lmdif_ (
	     void (*fcn)(int*, int* , double*, double*, int*),
	     int* m,
	     int* n,
	     double* x,
	     double* fvec,
	     double* ftol,
	     double* xtol,
	     double* gtol,
	     int* maxfev,
	     double* epsfcn,
	     double* diag,
	     int* mode,
	     double* factor,
	     int* nprint,
	     int* info,
	     int* nfev,
	     double* fjac,
	     int* ldfjac,
	     int* ipvt,
	     double* qtf,
	     double* wa1,
	     double* wa2,
	     double* wa3,
	     double* wa4
	     );

#endif

