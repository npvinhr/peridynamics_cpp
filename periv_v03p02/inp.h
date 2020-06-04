#include <math.h>

#ifndef INP_H
#define INP_H

namespace inp
{
	const double	pi(3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862);

	// Center crack
	const double	length	(0.05e0);			//length: Total length of the plate - x direction
	const double	width	(0.05e0);			//width : Total width of the plate - y direction
	const int		ndivx	(50);		//ndivx: Number of divisions in x direction - except boundary region
	const int		ndivy	(50);		//ndivy: Number of divisions in y direction - except boundary region
	// Edge crack
	//const double	length(0.1e0);			//length: Total length of the plate - x direction
	//const double	width(0.04e0);			//width : Total width of the plate - y direction
	//const int		ndivx	(100*5/2);	//ndivx: Number of divisions in x direction - except boundary region
	//const int		ndivy	(100);		//ndivy: Number of divisions in y direction - except boundary region

	//const double	length	(0.001e0*5.0);			//length: Total length of the plate - x direction
	//const double	width	(0.001e0*15.0);			//width : Total width of the plate - y direction
	//const int		ndivx	(5);	//ndivx: Number of divisions in x direction - except boundary region
	//const int		ndivy	(15);		//ndivy: Number of divisions in y direction - except boundary region
	
	const int		nbnd	(3);		//nbnd: Number of divisions in the boundary region
	const int		totnode	(ndivx*(ndivy + 2 * nbnd));	//totnode: Total number of material points
	//const int		totnode(2);	//totnode: Total number of material points

	const double	dx		(length / ndivx);	//dx : Spacing between material points
	//const double	dx(1.0);	//dx : Spacing between material points
	const double	delta	(3.015 * dx);		//delta : Horizon
	const double	thick	(dx);				//thick : Thickness of the plate						
	const double	area	(dx * dx);			//area : Cross - sectional area
	const double	vol		(area * dx);		//vol : Volume of a material point
	const double	radij	(dx / 2.0e0);		//radij : Material point radius
	
	const double	bulkModulus	(163.0e9);		//bulkModulus : bulk modulus
	const double	shearModulus(79.0e9);		//shearModulus : shear modulus
	const double	yieldStress	(250.0e6*1e3);		//yieldStress: yield stress

	const double	dens	(8000.0e0);			//dens : Density
	const double	emod	(192.0e9);			//emod : Elastic modulus
	const double	pratio	(1.0e0 / 3.0e0);	//pratio12 = Poisson's ratio						
	const double	crlength(0.01e0);			//crlength : Crack length
	const double	scr0	(0.04472e0);		//scr0 : Critical stretch
									
	const double	bc		(9.0e0 * emod / (pi * thick * (delta*delta*delta)));	//bc : Bond constant
	const double	sedload1(9.0e0 / 16.0e0 * emod * 1.0e-6);				//sedload1 : Strain energy density for the first loading
	const double	sedload2(9.0e0 / 16.0e0 * emod * 1.0e-6);				//sedload2 : Strain energy density for the second loading
	
	const double	dt		(0.8e0 * sqrt(2.0e0*dens*dx / (pi*delta*delta * dx*bc)));	//dt : Time interval
	const int		nt		(1000);			//nt: Total number of time step
	const double	totime	(nt * dt);			//totime : Total time
	const int		maxfam	(100);				//maxfam: Maximum number of material points inside a horizon of a material point
	const int		nodefammax(10000000);
}
#endif
