#ifndef PERIV_UTILITIES_H
#define PERIV_UTILITIES_H

//Influence function
double omega(double zeta);

void computeWeightedVolume(double *coord, int *pointfam, int *numfam, int *nodefam, int *fail, double *fncst, double *wvolume);

void computeDilatation(double *coord, double *disp, int *pointfam, int *numfam, int *nodefam, int *fail,
	double *fncst, double *wvolume, double *thetai);


#endif