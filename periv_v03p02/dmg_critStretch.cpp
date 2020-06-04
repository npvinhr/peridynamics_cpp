#include <cmath>
#include "inp.h"

void dmg_critStretch(double *coord, double *disp, int *pointfam, int *numfam, int *nodefam, int *fail, double *dmg)
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
	
	

	for (int i = 0; i < totnode; ++i)
	{
		double dmgpar1 = 0.0e0;
		double dmgpar2 = 0.0e0;
		for (int j = 0; j < numfam[i]; ++j)
		{
			cnode = nodefam[pointfam[i] + j];
			idist = sqrt((coord[2 * cnode + 0] - coord[2 * i + 0])*(coord[2 * cnode + 0] - coord[2 * i + 0]) +
				(coord[2 * cnode + 1] - coord[2 * i + 1])*(coord[2 * cnode + 1] - coord[2 * i + 1]));
			nlength = sqrt((coord[2 * cnode + 0] + disp[2 * cnode + 0] - coord[2 * i + 0] - disp[2 * i + 0])*(coord[2 * cnode + 0] + disp[2 * cnode + 0] - coord[2 * i + 0] - disp[2 * i + 0]) +
				(coord[2 * cnode + 1] + disp[2 * cnode + 1] - coord[2 * i + 1] - disp[2 * i + 1])*(coord[2 * cnode + 1] + disp[2 * cnode + 1] - coord[2 * i + 1] - disp[2 * i + 1]));

			//Volume correction		//[Peridynamic Theory and Its Applications - 7.2]
			if (idist <= delta - radij)
				fac = 1.0e0;
			else if (idist <= delta + radij)
				fac = (delta + radij - idist) / (2.0e0*radij);
			else
				fac = 0.0e0;

			//Definition of a no - fail zone
			if (abs((nlength - idist) / idist) > scr0)
				if (abs(coord[2 * i + 1]) <= (length / 4.0e0))
					fail[maxfam*i + j] = 0;

			dmgpar1 = dmgpar1 + fail[maxfam*i + j] * vol * fac;
			dmgpar2 = dmgpar2 + vol * fac;
		}
		dmg[i] = 1.0e0 - dmgpar1 / dmgpar2;
	}
	
}