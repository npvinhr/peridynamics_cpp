//#include <iostream>
//#include <fstream>
//#include <string>

#include <cmath>
#include "inp.h"
#include "periv_utilities.h"

void periv_ElasticPlasticMaterial(const double *coord, const double *disp, const int *pointfam, const int *numfam, const int *nodefam, int *fail,
	const double *fncst, const double *wvolume, const double *thetai, double *edp, double *pforce)
{
	using namespace inp;

	//idist : Initial distance
	double idist = 0.0;
	//cnode : Current material point
	int cnode = 0;
	//Length of deformed bond
	double nlength = 0.0;
	//fac: Volume correction factor
	double fac = 0.0;
	//dforce1 : x component of the PD force between two material points
	double dforce1 = 0.0;
	//dforce1 : y component of the PD force between two material points
	double dforce2 = 0.0;
	//ed : deviatoric part of the extension (e^d)
	double ed = 0.0;
	//t_u : scalar force state
	double t_u = 0.0;
	//M1 : x component of the deformed direction vector state
	double M1 = 0.0;
	//M2 : y component of the deformed direction vector state
	double M2 = 0.0;
	//OMEGAP : unit influence function
	double OMEGAP = 1.0;
	//double thickness = 1.0;


	double yieldValue(1.123e1);
	//2D yield value (uniaxial stress)
	//yieldValue = 225.0 / 3.0 * yieldStress * yieldStress / 8 / pi / thick / pow(delta, 4);
	//3D yield value (uniaxial stress)
	//yieldValue = 25.0 * yieldStress * yieldStress / 8 / pi / pow(delta, 5);


	for (int i = 0; i < totnode; ++i)
	{
		pforce[2 * i + 0] = 0.0;
		pforce[2 * i + 1] = 0.0;
	}

	double tdNorm	= 0.0;
	double tdTrial	= 0.0;
	double edpN		= 0.0;
	for (int i = 0; i < totnode; ++i)
	{
		double alpha = 15.0 * shearModulus / wvolume[i];
		
		//Compute Deviatoric Force State Norm
		tdNorm = 0.0;
		for (int j = 0; j < numfam[i]; ++j)
		{
			cnode = nodefam[pointfam[i] + j];
			if (fail[maxfam*i + j] == 1)
			{
				idist = sqrt((coord[2 * cnode + 0] - coord[2 * i + 0])*(coord[2 * cnode + 0] - coord[2 * i + 0]) +
					(coord[2 * cnode + 1] - coord[2 * i + 1])*(coord[2 * cnode + 1] - coord[2 * i + 1]));
				nlength = sqrt((coord[2 * cnode + 0] + disp[2 * cnode + 0] - coord[2 * i + 0] - disp[2 * i + 0])*(coord[2 * cnode + 0] + disp[2 * cnode + 0] - coord[2 * i + 0] - disp[2 * i + 0]) +
					(coord[2 * cnode + 1] + disp[2 * cnode + 1] - coord[2 * i + 1] - disp[2 * i + 1])*(coord[2 * cnode + 1] + disp[2 * cnode + 1] - coord[2 * i + 1] - disp[2 * i + 1]));
				ed = (nlength - idist); //- thetai[i] * idist / 3.0;
				double test = thetai[i];
				edpN = edp[j];
				tdTrial = alpha * OMEGAP * (ed - edpN);

				//Volume correction		//[Peridynamic Theory and Its Applications - 7.2]
				if (idist <= delta - radij)
					fac = 1.0;
				else if (idist <= delta + radij)
					fac = (delta + radij - idist) / (2.0e0*radij);
				else
					fac = 0.0;

				double theta;
				if (abs(coord[2 * cnode + 1] - coord[2 * i + 1]) <= 1.0e-10)
					theta = 0.0;
				else if (abs(coord[2 * cnode + 0] - coord[2 * i + 0]) <= 1.0e-10)
					theta = 90.0 * pi / 180.0;
				else
					theta = atan(abs(coord[2 * cnode + 1] - coord[2 * i + 1]) / abs(coord[2 * cnode + 0] - coord[2 * i + 0]));

				//Determination of the surface correction between two material points		//[Peridynamic Theory and Its Applications - 7.7]
				double scx, scy, scr;
				scx = (fncst[2 * i + 0] + fncst[2 * cnode + 0]) / 2.0;
				scy = (fncst[2 * i + 1] + fncst[2 * cnode + 1]) / 2.0;
				scr = 1.0 / (((cos(theta))*(cos(theta)) / (scx*scx)) + ((sin(theta))*(sin(theta)) / (scy*scy)));
				scr = sqrt(scr);

				tdNorm = tdNorm + tdTrial * tdTrial * vol * scr * fac;
				double p = 1;
			}
		}	
		tdNorm = sqrt(tdNorm);
		//Compute Deviatoric Force State Norm
		
		
		//Yield function
		double f;
		f = tdNorm * tdNorm / 2.0 - yieldValue;
		double p = 1;

		//Plastic condition
		double deltaLambda;
		if (f > 0)		
		{
			deltaLambda = (tdNorm / sqrt(2.0*yieldValue) - 1.0) / alpha;
		}
		
		
		//Peridynamic force state
		for (int j = 0; j < numfam[i]; ++j)
		{
			cnode = nodefam[pointfam[i] + j];
			if (fail[maxfam*i + j] == 1)
			{
				idist = sqrt((coord[2 * cnode + 0] - coord[2 * i + 0])*(coord[2 * cnode + 0] - coord[2 * i + 0]) +
					(coord[2 * cnode + 1] - coord[2 * i + 1])*(coord[2 * cnode + 1] - coord[2 * i + 1]));
				nlength = sqrt((coord[2 * cnode + 0] + disp[2 * cnode + 0] - coord[2 * i + 0] - disp[2 * i + 0])*(coord[2 * cnode + 0] + disp[2 * cnode + 0] - coord[2 * i + 0] - disp[2 * i + 0]) +
					(coord[2 * cnode + 1] + disp[2 * cnode + 1] - coord[2 * i + 1] - disp[2 * i + 1])*(coord[2 * cnode + 1] + disp[2 * cnode + 1] - coord[2 * i + 1] - disp[2 * i + 1]));
				ed = (nlength - idist); //-thetai[i] * idist / 3.0;
				edpN = edp[j];
				tdTrial = alpha * OMEGAP * (ed - edpN);

				//Deviatoric part of force state
				double td;
				if (f > 0)		//Plastic condition
				{
					td		= sqrt(2.0*yieldValue) * tdTrial / tdNorm;
					//double testf = td * td / 2.0 - yieldValue;
					edp[j]	= edpN + td * deltaLambda;
				}
				else
				{
					td		= tdTrial;
					edp[j]	= edpN;
				}

				//Isotropic part of force state
				double ti;
				ti = 3.0 / wvolume[i] * bulkModulus*thetai[i] * OMEGAP * idist;

				//Force state
				double t_u;
				t_u = ti + td;
				M1 = (coord[2 * cnode + 0] + disp[2 * cnode + 0] - coord[2 * i + 0] - disp[2 * i + 0]) / nlength;
				M2 = (coord[2 * cnode + 1] + disp[2 * cnode + 1] - coord[2 * i + 1] - disp[2 * i + 1]) / nlength;

				//Volume correction		//[Peridynamic Theory and Its Applications - 7.2]
				if (idist <= delta - radij)
					fac = 1.0;
				else if (idist <= delta + radij)
					fac = (delta + radij - idist) / (2.0*radij);
				else
					fac = 0.0;

				double theta;
				if (abs(coord[2 * cnode + 1] - coord[2 * i + 1]) <= 1.0e-10)
					theta = 0.0;
				else if (abs(coord[2 * cnode + 0] - coord[2 * i + 0]) <= 1.0e-10)
					theta = 90.0 * pi / 180.0;
				else
					theta = atan(abs(coord[2 * cnode + 1] - coord[2 * i + 1]) / abs(coord[2 * cnode + 0] - coord[2 * i + 0]));

				//Determination of the surface correction between two material points		//[Peridynamic Theory and Its Applications - 7.7]
				double scx, scy, scr;
				scx = (fncst[2 * i + 0] + fncst[2 * cnode + 0]) / 2.0;
				scy = (fncst[2 * i + 1] + fncst[2 * cnode + 1]) / 2.0;
				scr = 1.0 / (((cos(theta))*(cos(theta)) / (scx*scx)) + ((sin(theta))*(sin(theta)) / (scy*scy)));
				scr = sqrt(scr);

				//Calculation of the peridynamic force in x and y directions
				//acting on a material point i due to a material point j
				dforce1 = t_u * M1 * vol * scr * fac;
				dforce2 = t_u * M2 * vol * scr * fac;
				//dforce1 = td * M1 * vol * scr * fac;
				//dforce2 = td * M2 * vol * scr * fac;
			}
			else
			{
				dforce1 = 0.0;
				dforce2 = 0.0;
			}

			

			pforce[2 * i + 0]		= pforce[2 * i + 0]		+ dforce1;
			pforce[2 * i + 1]		= pforce[2 * i + 1]		+ dforce2;

			pforce[2 * cnode + 0]	= pforce[2 * cnode + 0] - dforce1;
			pforce[2 * cnode + 1]	= pforce[2 * cnode + 1] - dforce2;
		}
	}
}