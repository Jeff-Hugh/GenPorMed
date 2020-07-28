#ifndef POROUS3D_H
#define POROUS3D_H

#include<string>
#include<vector>

using namespace std;

class Porous3D
{
public:
	const int NX;         /// grid number in x direction
	const int NY;         /// grid number in y direction
	const int NZ;         /// grid number in y direction
	int* Solid;      /// the return value
	const double phi;	/// target porosity
	const double p_cd;		/// core distribution probability
	const double p_surface;     /// probability to grow on surface
	const double p_edge;	/// probability to grow on edge, default value (p_surface / 12.0)
	const double p_point;	/// probability to grow on point, default value (p_surface / 96.0)
	
	/// Constructor
	Porous3D(int nx, int ny, int nz, int *s, double por, double p, double p_s, double p_e, double p_p);
	~Porous3D();

	void output2tecplot(std::string filename);
	bool smooth();
	bool fixhole();
	double CalcPor();
	int NearSolid(int i, int j, int k);

private:
	int Grow_Times;         /// for iteration

	void grow();
	bool is_interior(const int i, const int j, const int k);
	int CellIndex(int x, int y, int z);

};

#endif
