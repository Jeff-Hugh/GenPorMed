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
	double phi;	/// target porosity
	const double p_cd;		/// core distribution probability
	const double p_surface;     /// probability to grow on surface
	const double p_edge;	/// probability to grow on edge, default value (p_surface / 12.0)
	const double p_point;	/// probability to grow on point, default value (p_surface / 96.0)
	
	/// Constructor
	Porous3D(int nx, int ny, int nz, double por, double p, double p_s, double p_e, double p_p);
	~Porous3D();

	void Generation(int* s);
	void output2tecplot(std::string filename);
	double CalcPor();
	int NearSolid(int i, int j, int k);
	void DeleteDeadZone();

private:
	int Grow_Times;         /// for iteration

	void grow();
	int CellIndex(int x, int y, int z);
	void InfectFluid(int i, int j, int k);
	int DifferentNumber(int i, int j, int k, int core_id);

};

#endif
