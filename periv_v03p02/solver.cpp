#include <iostream>
#include <fstream>
#include <string>

#include <cmath>
#include "inp.h"
#include "periv_PropotypeMicroelastic.h"
#include "dmg_critStretch.h"
#include "periv_utilities.h"
#include "periv_LinearPeridynamicSolid.h"
#include "periv_ElasticPlasticMaterial.h"



void solver(double *coord, double *acc, double *vel, double *disp,
	int *pointfam, int *numfam, int *nodefam, int *fail,
	double *fncst, double *dmg, double *endtime, double *wvolume, double *thetai, double *edp, double *pforce, double *bforce,
	int totint, int totbottom, int tottop)
{
	using namespace inp;

	//ctime : Current time
	double ctime = 0.0e0;
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

	//double	*pforcetest = new double[nt];
	
	//Initialization of displacements and velocities
	for (int i = 0; i < totnode; ++i)
	{
		vel[2 * i + 0] = 0.0e0;
		vel[2 * i + 1] = 0.0e0;
		disp[2 * i + 0] = 0.0e0;
		disp[2 * i + 1] = 0.0e0;
	}


	//For calculation of peridynamics forces - LPS
	computeWeightedVolume(coord, pointfam, numfam, nodefam, fail, fncst, wvolume);



	//Time integration
	for (int tt = 0; tt < nt + 1; ++tt)
	{
		std::cout << "tt = " << tt << std::endl;
		ctime = tt * dt;

		double srate = 10.0e0;

		//Application of boundary conditions at the top and bottom edges
		for (int i = totint; i < totbottom; ++i)
		{
			vel [2 * i + 1] = -srate;
			disp[2 * i + 1] = -srate * tt * dt;
		}
		for (int i = totbottom; i < tottop; ++i)
		{
			vel [2 * i + 1] = srate;
			disp[2 * i + 1] = srate * tt * dt;
		}
		
		

		

		


		//Calculation of peridynamics forces - PMB (with dmg_critStretch)
		//periv_PropotypeMicroelastic(coord, disp, pointfam, numfam, nodefam, fail, fncst, pforce);
		

		//Calculation of peridynamics forces - LPS
		/*computeWeightedVolume(coord, pointfam, numfam, nodefam, fncst, wvolume);*/		//Put outside of time integration loop
		computeDilatation(coord, disp, pointfam, numfam, nodefam, fail, fncst, wvolume, thetai);
		periv_LinearPeridynamicSolid(coord, disp, pointfam, numfam, nodefam, fail,	fncst, wvolume, thetai, pforce);


		//Calculation of peridynamics forces - Elastic-Plastic
		//computeDilatation(coord, disp, pointfam, numfam, nodefam, fail, fncst, wvolume, thetai);
		//periv_ElasticPlasticMaterial(coord, disp, pointfam, numfam, nodefam, fail, fncst, wvolume, thetai, edp, pforce);


		//Bond damage and failure
		dmg_critStretch(coord, disp, pointfam, numfam, nodefam, fail, dmg);


		



		for (int i = 0; i < totint; ++i)
		{
			//Calculation of acceleration of material point i
			acc[2 * i + 0] = (pforce[2 * i + 0] + bforce[2 * i + 0]) / dens;
			acc[2 * i + 1] = (pforce[2 * i + 1] + bforce[2 * i + 1]) / dens;
			//Calculation of velocity of material point i
			//by integrating the acceleration of material point i
			vel[2 * i + 0] = vel[2 * i + 0] + acc[2 * i + 0] * dt;
			vel[2 * i + 1] = vel[2 * i + 1] + acc[2 * i + 1] * dt;
			//Calculation of displacement of material point i
			//by integrating the velocity of material point i
			disp[2 * i + 0] = disp[2 * i + 0] + vel[2 * i + 0] * dt;
			disp[2 * i + 1] = disp[2 * i + 1] + vel[2 * i + 1] * dt;
		}
		for (int i = totint; i < totbottom; ++i)
		{
			acc[2 * i + 0] = (pforce[2 * i + 0] + bforce[2 * i + 0]) / dens;
			vel[2 * i + 0] = vel[2 * i + 0] + acc[2 * i + 0] * dt;
			disp[2 * i + 0] = disp[2 * i + 0] + vel[2 * i + 0] * dt;
		}
		for (int i = totbottom; i < tottop; ++i)
		{
			acc[2 * i + 0] = (pforce[2 * i + 0] + bforce[2 * i + 0]) / dens;
			vel[2 * i + 0] = vel[2 * i + 0] + acc[2 * i + 0] * dt;
			disp[2 * i + 0] = disp[2 * i + 0] + vel[2 * i + 0] * dt;
		}

		


		endtime[tt] = ctime;

		//pforcetest[tt] = pforce[2*(totnode/2+ndivx/2)+1];
		/*pforcetest[tt] = 0.0;
		for (int ilay = ndivx * ndivy / 2; ilay < ndivx*(ndivy / 2 + 1); ++ilay)
		{
			pforcetest[tt] = pforcetest[tt] + pforce[2 * ilay + 1];
		}*/


		//Output results to file - extended xyz format for Ovito
		std::string filename;
		if ((tt % 10) == 0)
		{
			using namespace std;
			//filename = "peri" + std::to_string(tt);

			char ttstring[10];
			sprintf_s(ttstring, "%06d", tt);
			std::string ttstr(ttstring);
			std::string filename = "periv" + ttstr + ".xyz";
			std::cout << filename << std::endl;

			ofstream outf(filename);

			outf << totint << endl;
			outf << "Properties=species:S:1:pos:R:3:dmg:R:1" << endl;
			for (int i = 0; i < totint; ++i)
				//outf << coord[2*i+0] << "     " << coord[2*i+1] << "     " << disp[2*i+0] << "     " << disp[2*i+1] << "     " << dmg[2*i+0] << endl;
				//Extended xyz format
				outf << "P     " << coord[2 * i + 0] + disp[2 * i + 0] << "     " << coord[2 * i + 1] + disp[2 * i + 1] << "     " << "0.00000" << "     " << dmg[i] << endl;
			outf.close();
		}
	}

	
	//std::string filename = "pforcetest.txt";
	//std::ofstream outf(filename);
	////outf << totint << std::endl;
	////outf << "Properties=species:S:1:pos:R:3:dmg:R:1" << endl;
	//for (int tt = 0; tt < nt; tt=tt+10)
	//	outf << tt << "       " << pforcetest[tt] << std::endl;
	//outf.close();
}