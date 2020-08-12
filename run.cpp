/// \reference Wang, M., Wang, J., Pan, N., & Chen, S. (2007). 
/// Mesoscopic predictions of the effective thermal conductivity for microscale random porous media. 
/// Physical Review E - Statistical, Nonlinear, and Soft Matter Physics, 75(3), 1¨C10.
/// https://doi.org/10.1103/PhysRevE.75.036702
/// 

#include <time.h>
#include <iostream>
#include <cstring>

#include "Porous2D.h"
#include "Porous3D.h"

using namespace std;

void Generate2D();
void Generate3D();

int main()
{
	clock_t startTime, endTime;
	startTime = clock();

	Generate3D();

	endTime = clock();
	cout << "Totle Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}

void Generate2D()
{
	const int M = 200;
	const int N = 200;
	const double phi = 0.5;			/// target porosity
	const double p_cd = 0.005;		/// core distribution probability
	const double z = 0.05;     /// Probability in orthogonal direction
	const double f = 0.0125;   /// Probability in oblique direction
	
	int* s = new int[M * N];
	memset(s, 0, sizeof(int)*M*N);		/// initialize to 0

	Porous2D porous(M, N, phi, p_cd, z, f);
	porous.Generation(s);

	porous.output2tecplot(M, N, s, "test.plt"); 

	delete[] s;
}

void Generate3D()
{
	const int Nx = 50;
	const int Ny = 50;
	const int Nz = 50;
	const double phi = 0.5;			/// target porosity
	const double p_cd = 0.005;		/// core distribution probability
	const double p_surface = 0.2;     /// probability to grow on surface
	const double p_edge = 0.0167;	  /// probability to grow on edge, default value (p_surface / 12.0)
	const double p_point = 0.0021;	/// probability to grow on point, default value (p_surface / 96.0)
	
	int* s = new int[Nx * Ny * Nz];
	memset(s, 0, sizeof(int)*Nx*Ny*Nz);		/// initialize to 0

	Porous3D porous(Nx, Ny, Nz, phi, p_cd, p_surface, p_edge, p_point);
	porous.Generation(s);

	porous.output2tecplot("test.plt"); 

	delete[] s;
}