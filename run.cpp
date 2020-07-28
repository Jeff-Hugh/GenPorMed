/// \reference Wang, M., Wang, J., Pan, N., & Chen, S. (2007). 
/// Mesoscopic predictions of the effective thermal conductivity for microscale random porous media. 
/// Physical Review E - Statistical, Nonlinear, and Soft Matter Physics, 75(3), 1¨C10.
/// https://doi.org/10.1103/PhysRevE.75.036702
/// 

#include <string>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <assert.h>

#include "Porous2D.h"

using namespace std;

int main()
{
	clock_t startTime, endTime;
	startTime = clock();

	const int M = 200;
	const int N = 200;
	const double phi = 0.3;			/// target porosity
	const double p_cd = 0.02;		/// core distribution probability
	const double z = 0.05;     /// Probability in orthogonal direction
	const double f = 0.0125;   /// Probability in oblique direction
	int* s = new int[M * N]{ 0 };

	Porous2D porous(M, N, phi, p_cd, z, f);
	porous.generation(s);

	porous.output2tecplot(M, N, s, "test.plt"); 

	delete[] s;
	endTime = clock();
	cout << "Totle Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
	return 0;
}