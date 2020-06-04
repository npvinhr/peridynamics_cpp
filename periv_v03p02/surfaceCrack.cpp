#include <iostream>
#include <fstream>
#include <string>

#include <cmath>
#include "inp.h"

void surfaceCrack(double *coord, int *pointfam, int *numfam, int *nodefam, int *fail)
{
	//Definition of the crack surface
	//PD bonds penetrating through the crack surface are broken

	using namespace inp;

	// Center crack
	//cnode : Current material point
	int cnode = 0;
	for (int i = 0; i < totnode; ++i)
	{
		for (int j = 0; j < numfam[i]; ++j)
		{
			//cnode = nodefam[pointfam[2*i+0] + j - 1];
			cnode = nodefam[pointfam[i] + j];
			if ((coord[2*cnode+1] > 0.0e0) and (coord[2*i+1] < 0.0e0))
			{
				if ((abs(coord[2*i+0]) - (crlength / 2.0e0)) <= 1.0e-10)
					fail[maxfam*i+j] = 0;
				else if ((abs(coord[2*cnode+0]) - (crlength / 2.0e0)) <= 1.0e-10)
					fail[maxfam*i+j] = 0;
			}
			else if ((coord[2*i+1] > 0.0e0) and (coord[2*cnode+1] < 0.0e0))
			{
				if ((abs(coord[2*i+0]) - (crlength / 2.0e0)) <= 1.0e-10)
					fail[maxfam*i+j] = 0;
				else if ((abs(coord[2*cnode+0]) - (crlength / 2.0e0)) <= 1.0e-10)
					fail[maxfam*i+j] = 0;
			}
		}
	}

	//// Edge crack
	////cnode : Current material point
	//int cnode = 0;
	//double cracktip = 1.0e-10;
	////double cracktip = -0.0125;
	//for (int i = 0; i < totnode; ++i)
	//{
	//	for (int j = 0; j < numfam[i]; ++j)
	//	{
	//		cnode = nodefam[pointfam[i] + j];
	//		if ((coord[2 * cnode + 1] > 0.0e0) and (coord[2 * i + 1] < 0.0e0))
	//		{
	//			if (coord[2 * i + 0] <= cracktip)
	//				fail[maxfam*i + j] = 0;
	//			else if (coord[2 * cnode + 0] <= cracktip)
	//				fail[maxfam*i + j] = 0;
	//		}
	//		else if ((coord[2 * i + 1] > 0.0e0) and (coord[2 * cnode + 1] < 0.0e0))
	//		{
	//			if (coord[2 * i + 0] <= cracktip)
	//				fail[maxfam*i + j] = 0;
	//			else if (coord[2 * cnode + 0] <= cracktip)
	//				fail[maxfam*i + j] = 0;
	//		}
	//	}
	//}

	/*std::string filename = "\test\surface.txt";
	std::ofstream outf(filename);

	for (int i = 0; i < totnode; ++i)
		for (int j = 0; j < numfam[i]; ++j)
		{
			outf << fail[idc(i, 0, maxfam)] << std::endl;
		}
	outf.close();*/
}