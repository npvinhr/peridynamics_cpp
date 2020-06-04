#include <iostream>
#include <fstream>
#include <string>

#include <cmath>
#include "inp.h"

void matPoint(double *coord, int &totintOut, int &totbottomOut, int &tottopOut)		//totint, totbottom, tottot will be updated
{
	//Specification of the locations of material points
	//Material points of the internal region

	using namespace inp;

	//nnum : Material point number
	int nnum = 0;

	for (int i = 0; i < ndivy; ++i)
	{
		for (int j = 0; j < ndivx; ++j)
		{
			nnum = nnum + 1;
			//coord[nnum][0] = (-1.0e0 * length / 2.0e0) + (dx / 2.0e0) + (j - 1) * dx;
			//coord[nnum][1] = (-1.0e0 * width / 2.0e0)  + (dx / 2.0e0) + (i - 1) * dx;
			coord[2*(nnum-1)+0] = (-1.0e0 * length / 2.0e0) + (dx / 2.0e0) + j * dx;
			coord[2*(nnum-1)+1] = (-1.0e0 * width  / 2.0e0) + (dx / 2.0e0) + i * dx;
		}
	}

	totintOut = nnum;

	//Material points of the boundary region - bottom
	for (int i = 0; i < nbnd; ++i)
	{
		for (int j = 0; j < ndivx; ++j)
		{
			nnum = nnum + 1;
			//coord[nnum][0] = -1.0e0 / 2.0e0 * length + (dx / 2.0e0) + (j - 1) * dx;
			//coord[nnum][1] = -1.0e0 / 2.0e0 * width -  (dx / 2.0e0) - (i - 1) * dx;
			coord[2*(nnum-1)+0] = -1.0e0 / 2.0e0 * length + (dx / 2.0e0) + j * dx;
			coord[2*(nnum-1)+1] = -1.0e0 / 2.0e0 * width  - (dx / 2.0e0) - i * dx;
		}
	}

	totbottomOut = nnum;

	//Material points of the boundary region - top
	for (int i = 0; i < nbnd; ++i)
	{
		for (int j = 0; j < ndivx; ++j)
		{
			nnum = nnum + 1;
			//coord[nnum][0] = -1.0e0 / 2.0e0 * length + (dx / 2.0e0) + (j - 1) * dx;
			//coord[nnum][1] =  1.0e0 / 2.0e0 * width +  (dx / 2.0e0) + (i - 1) * dx;
			coord[2*(nnum-1)+0] = -1.0e0 / 2.0e0 * length + (dx / 2.0e0) + j * dx;
			coord[2*(nnum-1)+1] =  1.0e0 / 2.0e0 * width  + (dx / 2.0e0) + i * dx;
		}
	}

	tottopOut = nnum;

	std::string filename = "matPoint.txt";
	std::ofstream outf(filename);

	for (int i = 0; i < totnode; ++i)
		outf << "P     " << coord[2 * i + 0] << "     " << coord[2 * i + 1] << "     " << "0.00000"  << std::endl;
	outf.close();
}