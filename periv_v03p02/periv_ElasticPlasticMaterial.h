#ifndef PERIV_ELASTICPLASTICMATERIAL_H
#define PERIV_ELASTICPLASTICMATERIAL_H

void periv_ElasticPlasticMaterial(const double *coord, const double *disp, const int *pointfam, const int *numfam, const int *nodefam, int *fail,
	const double *fncst, const double *wvolume, const double *thetai, double *edp, double *pforce);

#endif