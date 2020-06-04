#include <iostream>
#include <fstream>
#include <string>

#include <cmath>
#include "inp.h"

void horizon(double *coord, int *pointfam, int *numfam, int *nodefam)
{
	//Determination of material points inside the horizon of each material point

	using namespace inp;

	//idist : Initial distance
	double idist = 0.0e0;

	for (int i = 0; i < totnode; ++i)
	{
		char istring[10];
		sprintf_s(istring, "%06d", i);
		std::string istr(istring);
		std::string istep = "i=" + istr;
		//istep = "i" + istr;
		std::cout << istep << std::endl;

		if (i == 0)
			pointfam[i] = 0;
		else
			pointfam[i] = pointfam[i - 1] + numfam[i - 1];
		for (int j = 0; j < totnode; ++j)
		{
			idist = sqrt((coord[2*j+0] - coord[2*i+0])*(coord[2*j+0] - coord[2*i+0]) + (coord[2*j+1] - coord[2*i+1])*(coord[2*j+1] - coord[2*i+1]));
			if (i != j)
			{
				if (idist < delta)
				{
					numfam[i] = numfam[i] + 1;
					nodefam[pointfam[i] + numfam[i] - 1] = j;
				}
			}
		}
	}

	std::string filename = "horizon.txt";
	std::ofstream outf(filename);

	for (int i = 0; i < 100; ++i)
		outf << nodefam[i] << std::endl;
	outf.close();
}