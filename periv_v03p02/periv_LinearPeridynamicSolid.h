#ifndef PERIV_LINEARPERIDYNAMICSOLID_H
#define PERIV_LINEARPERIDYNAMICSOLID_H

void periv_LinearPeridynamicSolid(const double *coord, double *disp, int *pointfam, int *numfam, int *nodefam, int *fail,
	double *fncst, double *wvolume, double *thetai, double *pforce);

#endif