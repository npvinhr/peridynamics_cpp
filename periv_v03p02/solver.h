#ifndef SOLVER_H
#define SOLVER_H

void solver(double *coord, double *acc, double *vel, double *disp,
	int *pointfam, int *numfam, int *nodefam, int *fail,
	double *fncst, double *dmg, double *endtime, double *wvolume, double *thetai, double *edp, double *pforce, double *bforce,
	int totint, int totbottom, int tottop);

#endif