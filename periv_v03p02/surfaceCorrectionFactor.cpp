#include <cmath>
#include "inp.h"

void surfaceCorrectionFactor(double *coord, double *disp, int *pointfam, int *numfam, int *nodefam, double *stendens, double *fncst)
{
	using namespace inp;

	//idist : Initial distance
	double idist = 0.0e0;
	//cnode : Current material point
	int cnode = 0;
	//Length of deformed bond
	double nlength = 0.0e0;
	//fac: Volume correction factor
	double fac = 0.0e0;
	
	//Consider bond breaking or not? - No
	//Test loading condition on the actual structure

	//-------------------------------------------------------------------------------------------------//
	//Loading 1
	for (int i = 0; i < totnode; ++i)
	{
		disp[2*i+0] = 0.001e0 * coord[2*i+0];
		disp[2*i+1] = 0.0e0;
	}

	for (int i = 0; i < totnode; ++i)
	{
		stendens[2*i+0] = 0.0e0;
		for (int j = 0; j < numfam[i]; ++j)
		{
			//cnode = nodefam[pointfam[2*i+0] + j - 1][0];
			cnode = nodefam[pointfam[i] + j];
			idist = sqrt((coord[2*cnode+0] - coord[2*i+0])*(coord[2*cnode+0] - coord[2*i+0]) +
				(coord[2*cnode+1] - coord[2*i+1])*(coord[2*cnode+1] - coord[2*i+1]));
			nlength = sqrt((coord[2*cnode+0] + disp[2*cnode+0] - coord[2*i+0] - disp[2*i+0])*(coord[2*cnode+0] + disp[2*cnode+0] - coord[2*i+0] - disp[2*i+0]) +
				(coord[2*cnode+1] + disp[2*cnode+1] - coord[2*i+1] - disp[2*i+1])*(coord[2*cnode+1] + disp[2*cnode+1] - coord[2*i+1] - disp[2*i+1]));
			if (idist <= delta - radij)
				fac = 1.0e0;
			else if (idist <= delta + radij)
				fac = (delta + radij - idist) / (2.0e0*radij);
			else
				fac = 0.0e0;
			stendens[2*i+0] = stendens[2*i+0] + 0.5e0 * 0.5e0 * bc * ((nlength - idist) / idist)*((nlength - idist) / idist) * idist * vol * fac;
		}
		//Calculation of surface correction factor in x direction
		//by finding the ratio of the analytical strain energy density value
		//to the strain energy density value obtained from PD Theory
		fncst[2*i+0] = sedload1 / stendens[2*i+0];
	}
	//-------------------------------------------------------------------------------------------------//



	//-------------------------------------------------------------------------------------------------//
	//Loading 2
	for (int i = 0; i < totnode; ++i)
	{
		disp[2*i+0] = 0.0e0;
		disp[2*i+1] = 0.001e0 * coord[2*i+1];
	}

	for (int i = 0; i < totnode; ++i)
	{
		stendens[2*i+1] = 0.0e0;
		for (int j = 0; j < numfam[i]; ++j)
		{
			//cnode = nodefam[pointfam[2*i+0] + j - 1][0];
			cnode = nodefam[pointfam[i] + j];
			idist = sqrt((coord[2*cnode+0] - coord[2*i+0])*(coord[2*cnode+0] - coord[2*i+0]) +
				(coord[2*cnode+1] - coord[2*i+1])*(coord[2*cnode+1] - coord[2*i+1]));
			nlength = sqrt((coord[2*cnode+0] + disp[2*cnode+0] - coord[2*i+0] - disp[2*i+0])*(coord[2*cnode+0] + disp[2*cnode+0] - coord[2*i+0] - disp[2*i+0]) +
				(coord[2*cnode+1] + disp[2*cnode+1] - coord[2*i+1] - disp[2*i+1])*(coord[2*cnode+1] + disp[2*cnode+1] - coord[2*i+1] - disp[2*i+1]));
			if (idist <= delta - radij)
				fac = 1.0e0;
			else if (idist <= delta + radij)
				fac = (delta + radij - idist) / (2.0e0*radij);
			else
				fac = 0.0e0;
			stendens[2*i+1] = stendens[2*i+0] + 0.5e0 * 0.5e0 * bc * ((nlength - idist) / idist)*((nlength - idist) / idist) * idist * vol * fac;
		}
		//Calculation of surface correction factor in x direction
		//by finding the ratio of the analytical strain energy density value
		//to the strain energy density value obtained from PD Theory
		fncst[2*i+1] = sedload1 / stendens[2*i+0];
	}
	//-------------------------------------------------------------------------------------------------//
}