//#include <iostream>
//#include <fstream>
//#include <string>

#include <cmath>
#include "inp.h"
#include "periv_utilities.h"

void periv_LinearPeridynamicSolid(const double *coord, double *disp, int *pointfam, int *numfam, int *nodefam, int *fail,
	double *fncst, double *wvolume, double *thetai, double *pforce)
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
	//dforce1 : x component of the PD force between two material points
	double dforce1 = 0.0e0;
	//dforce1 : y component of the PD force between two material points
	double dforce2 = 0.0e0;
	//e_d : deviatoric part of the extension (e^d)
	double e_d = 0.0e0;
	//t_u : scalar force state
	double t_u = 0.0e0;
	//M1 : x component of the deformed direction vector state
	double M1 = 0.0e0;
	//M2 : y component of the deformed direction vector state
	double M2 = 0.0e0;

	for (int i = 0; i < totnode; ++i)
	{
		pforce[2 * i + 0] = 0.0e0;
		pforce[2 * i + 1] = 0.0e0;
	}

	for (int i = 0; i < totnode; ++i)
	{
		for (int j = 0; j < numfam[i]; ++j)
		{
			cnode = nodefam[pointfam[i] + j];
			if (fail[maxfam*i + j] == 1)
			{
				idist = sqrt((coord[2 * cnode + 0] - coord[2 * i + 0])*(coord[2 * cnode + 0] - coord[2 * i + 0]) +
					(coord[2 * cnode + 1] - coord[2 * i + 1])*(coord[2 * cnode + 1] - coord[2 * i + 1]));
				nlength = sqrt((coord[2 * cnode + 0] + disp[2 * cnode + 0] - coord[2 * i + 0] - disp[2 * i + 0])*(coord[2 * cnode + 0] + disp[2 * cnode + 0] - coord[2 * i + 0] - disp[2 * i + 0]) +
					(coord[2 * cnode + 1] + disp[2 * cnode + 1] - coord[2 * i + 1] - disp[2 * i + 1])*(coord[2 * cnode + 1] + disp[2 * cnode + 1] - coord[2 * i + 1] - disp[2 * i + 1]));
				e_d = (nlength - idist) - thetai[i]*idist/3.0e0;
				t_u = (3.0e0/wvolume[i]*bulkModulus*thetai[i]*idist + 15.0e0*shearModulus/wvolume[i]*e_d) * omega(idist);
				M1 = (coord[2 * cnode + 0] + disp[2 * cnode + 0] - coord[2 * i + 0] - disp[2 * i + 0])/nlength;
				M2 = (coord[2 * cnode + 1] + disp[2 * cnode + 1] - coord[2 * i + 1] - disp[2 * i + 1])/nlength;

				//Volume correction		//[Peridynamic Theory and Its Applications - 7.2]
				if (idist <= delta - radij)
					fac = 1.0e0;
				else if (idist <= delta + radij)
					fac = (delta + radij - idist) / (2.0e0*radij);
				else
					fac = 0.0e0;

				double theta;
				if (abs(coord[2 * cnode + 1] - coord[2 * i + 1]) <= 1.0e-10)
					theta = 0.0e0;
				else if (abs(coord[2 * cnode + 0] - coord[2 * i + 0]) <= 1.0e-10)
					theta = 90.0e0 * pi / 180.0e0;
				else
					theta = atan(abs(coord[2 * cnode + 1] - coord[2 * i + 1]) / abs(coord[2 * cnode + 0] - coord[2 * i + 0]));

				//Determination of the surface correction between two material points		//[Peridynamic Theory and Its Applications - 7.7]
				double scx, scy, scr;
				scx = (fncst[2 * i + 0] + fncst[2 * cnode + 0]) / 2.0e0;
				scy = (fncst[2 * i + 1] + fncst[2 * cnode + 1]) / 2.0e0;
				scr = 1.0e0 / (((cos(theta))*(cos(theta)) / (scx*scx)) + ((sin(theta))*(sin(theta)) / (scy*scy)));
				scr = sqrt(scr);

				//Calculation of the peridynamic force in x and y directions
				//acting on a material point i due to a material point j
				dforce1 = t_u * M1 * vol * scr * fac;
				dforce2 = t_u * M2 * vol * scr * fac;
			}
			else
			{
				dforce1 = 0.0e0;
				dforce2 = 0.0e0;
			}

			

			pforce[2 * i + 0]		= pforce[2 * i + 0]		+ dforce1;
			pforce[2 * i + 1]		= pforce[2 * i + 1]		+ dforce2;

			pforce[2 * cnode + 0]	= pforce[2 * cnode + 0] - dforce1;
			pforce[2 * cnode + 1]	= pforce[2 * cnode + 1] - dforce2;
		}
	}
}