#include <iostream>
#include <fstream>
#include <string>

#include "inp.h"
#include "matPoint.h"
#include "horizon.h"
#include "surfaceCrack.h"
#include "surfaceCorrectionFactor.h"
#include "solver.h"



int main()
{
	using namespace inp;

	std::cout << "Total number of node: " << totnode << std::endl;

	double	*coord		= new double[totnode * 2];
	double	*acc		= new double[totnode * 2];
	double	*vel		= new double[totnode * 2];
	double	*disp		= new double[totnode * 2];

	int		*pointfam	= new int[totnode];
	int		*numfam		= new int[totnode];
	int		*nodefam	= new int[nodefammax];
	int		*fail		= new int[totnode*maxfam];
	
	double	*stendens	= new double[totnode * 2];
	double	*fncst		= new double[totnode * 2];
	double	*pforce		= new double[totnode * 2];
	double	*bforce		= new double[totnode * 2];
	double	*dmg		= new double[totnode];
	
	double	*wvolume	= new double[totnode];
	double	*thetai		= new double[totnode];
	
	double	*enddisp	= new double[nt];
	double	*endtime	= new double[nt];
	double	*pforceold	= new double[totnode * 2];
	double	*velhalfold = new double[totnode * 2];
	double	*velhalf	= new double[totnode * 2];
	double	*massvec	= new double[totnode * 2];

	double	*edp		= new double[nodefammax];	//edp: e^dp in plasticity model


//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//
//INITIALIZATION//
	for (int i = 0; i < totnode; ++i)
	{
		//coord: Material point locations, 0 : x - coord, 1 : y - coord
		coord[2*i+0] = 0.0e0;
		coord[2*i+1] = 0.0e0;
		//numfam : Number of family members of each material point
		numfam[i] = 0;
		//pointfam : index array to find the family members in nodefam array
		pointfam[i] = 0;
		//pforce : total peridynamic force acting on a material point, 0 : x - coord, 1 : y - coord
		pforce[2*i+0] = 0.0e0;
		pforce[2*i+1] = 0.0e0;
		//pforceold : total peridynamic force acting on a material point in the previous time step
		//0 : x - coord, 1 : y - coord
		pforceold[2*i+0] = 0.0e0;
		pforceold[2*i+1] = 0.0e0;
		//bforce : body load acting on a material point, 0 : x - coord, 1 : y - coord
		bforce[2*i+0] = 0.0e0;
		bforce[2*i+1] = 0.0e0;
		//stendens : strain energy of a material point, 0 : loading 1, 1 : loading 2
		stendens[2*i+0] = 0.0e0;
		stendens[2*i+1] = 0.0e0;
		//fncst : surface correction factors of a material point, 0 : loading 1, 1 : loading 2
		fncst[2*i+0] = 1.0e0;
		fncst[2*i+1] = 1.0e0;
		//disp : displacement of a material point, 0 : x - coord, 1 : y - coord
		disp[2*i+0] = 0.0e0;
		disp[2*i+1] = 0.0e0;
		//vel : velocity of a material point, 0 : x - coord, 1 : y - coord
		vel[2*i+0] = 0.0e0;
		vel[2*i+1] = 0.0e0;
		velhalfold[2*i+0] = 0.0e0;
		velhalfold[2*i+1] = 0.0e0;
		velhalf[2*i+0] = 0.0e0;
		velhalf[2*i+1] = 0.0e0;
		//acc: acceleration of a material point, 0 : x - coord, 1 : y - coord
		acc[2*i+0] = 0.0e0;
		acc[2*i+1] = 0.0e0;
		//massvec : massvector for adaptive dynamic relaxation, 0 : x - coord, 1 : y - coord
		massvec[2*i+0] = 0.0e0;
		massvec[2*i+1] = 0.0e0;
		//fail : Failure array
		for (int j = 0; j < maxfam; ++j)
		{
			fail[maxfam*i+j] = 0;
		}
		//dmg: Damage of a material point
		dmg[i] = 0.0e0;
		wvolume[i] = 0.0e0;
		thetai[i] = 0.0e0;
	}

	for (int i = 0; i < nodefammax; ++i)
	{
		//nodefam: array containing family members of all material points
		nodefam[i] = 0;
		//edp: e^dp
		edp[i] = 0.0e0;
	}
	for (int i = 0; i < nt; ++i)
	{
		enddisp[i] = 0.0e0;
		endtime[i] = 0.0e0;
	}
	//Initialization of displacements and velocities
	for (int i = 0; i < totnode; ++i)
	{
		vel	[2*i+0] = 0.0e0;
		vel	[2*i+1] = 0.0e0;
		disp[2*i+0] = 0.0e0;
		disp[2*i+1] = 0.0e0;
	}
	//Initialization of fail flag array
	//1 means no failure, 0 means failure of the PD bond
	for (int i = 0; i < totnode; ++i)
	{
		for (int j = 0; j < maxfam; ++j)
		{
			fail[maxfam*i+j] = 1;
		}
	}
//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//



	




//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//
	int totint(0), totbottom(0), tottop(0);
	matPoint(coord, totint, totbottom, tottop);		//totint, totbottom, tottot will be updated
//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//	

		
	

//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//	
	horizon(coord, pointfam, numfam, nodefam);
//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//	

	

//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//	
	surfaceCrack(coord, pointfam, numfam, nodefam, fail);
//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//	

	

//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//	
	surfaceCorrectionFactor(coord, disp, pointfam, numfam, nodefam, stendens, fncst);
//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//	

	

//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//	
	solver(coord, acc, vel, disp,
		pointfam, numfam, nodefam, fail,
		fncst, dmg, endtime, wvolume, thetai, edp, pforce, bforce,
		totint, totbottom, tottop);
//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//	





	//std::string filename = "\results\periv.txt";
	//std::ofstream outf(filename);

	//outf << totint << std::endl;
	//outf << "Properties=species:S:1:pos:R:3:dmg:R:1" << std::endl;
	//for (int i = 1; i < totint + 1; ++i)
	//	//outf << coord[2*i+0] << "     " << coord[2*i+1] << "     " << disp[2*i+0] << "     " << disp[2*i+1] << "     " << dmg[2*i+0] << endl;
	//	//Extended xyz format
	//	outf << "P     " << coord[2*i+0] + disp[2*i+0] << "     " << coord[2*i+1] + disp[2*i+1] << "     " << "0.00000" << "     " << dmg[i] << std::endl;
	//outf.close();




	std::cout << "Total number of node: " << totnode << std::endl;

	std::cin.clear(); // reset any error flags
	std::cin.ignore(32767, '\n'); // ignore any characters in the input buffer until we find an enter character
	std::cin.get(); // get one more char from the user
	return 0;
}