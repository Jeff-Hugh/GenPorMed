#include <iostream>
#include <stdlib.h> 
#include <time.h>
#include <string>
#include <fstream>
#include <vector>
#include <omp.h>
#include <assert.h>

#include "Porous3D.h"

using namespace std;

Porous3D::Porous3D(int nx, int ny, int nz, double por = 0.6, double p = 0.05,
	double p_s = 0.2, double p_e = 0.0167, double p_p = 0.0021):NX(nx), NY(ny), NZ(nz), phi(por), p_cd(p), p_surface(p_s), p_edge(p_e), p_point(p_p)
{
}

void Porous3D::Generation(int* s)
{
	double phi_p;
	/// Generate random seed
	srand((unsigned)time(NULL));

	/// Initialization
	this->Solid = new int[NX*NY*NZ]{};

	// every core has different id
	int core_id = 1;

	/// Generate growth core
	for (int i = 0; i < this->NX; i++)
		for (int j = 0; j < this->NY; j++)
			for (int k = 0; k < this->NZ; k++)
			{
				const int cell = CellIndex(i, j, k);
				if (((rand() / double(RAND_MAX)) < p_cd) && (NearSolid(i,j,k) == 0)) this->Solid[cell] = core_id++;
			}
				
	///  Generate process
	this->Grow_Times = 0;
	do
	{
		grow();
		
		/// Calculate porosity
		phi_p = CalcPor();
		//std:: cout << "Porosity: " << phi_p << endl;
	} while (phi_p > phi);

// delete CoreID
#pragma omp parallel for
	for (int i = 0; i < NX * NY * NZ; i++)
	{
		if (Solid[i] > 0) Solid[i] = 1;
	}

	DeleteDeadZone();
	phi_p = CalcPor();

	phi = phi_p;
	std::cout << "Grow completed, Porosity: " << phi_p << endl;

	//copy Solid (inner) to s (outsider)
# pragma omp parallel for
	for (int i = 0; i < NX*NY*NZ; i++) s[i] = Solid[i];
		
}

void Porous3D::grow()
{
	Grow_Times += 1;
	// std::cout << "Grow Times: " << Grow_Times << endl;

# pragma omp parallel for 
	for (int i = 0; i < this->NX; i++)
		for (int j = 0; j < this->NY; j++)
			for (int k = 0; k < this->NZ; k++)
			{
				const int cell = CellIndex(i, j, k);
				if (this->Solid[cell] != 0)		/// if solid
				{
					int core_id = Solid[cell];
					for (int x = -1; x <= 1; x++)
						for (int y = -1; y <= 1; y++)
							for (int z = -1; z <= 1; z++)
							{
								if (!((i + x) < 0 || (j + y) < 0 || (k + z) < 0 || (i + x) > this->NX - 1 || (j + y) > this->NY - 1 || (k + z) > this->NZ - 1))
								{
									const int cell_n = CellIndex(i + x, j + y, k + z);
									switch (abs(x) + abs(y) + abs(z))
									{
									case 0:
										break;
									case 1:		// face neighbor
										// grow in X direction
										if ((x != 0) && (rand() / double(RAND_MAX)) < p_surface) 
										{
											int diff_num = DifferentNumber(i + x, j + y, k + z, core_id);
											if (diff_num == 0) Solid[cell_n] = core_id;
										}
											
										// grow in Y direction
										if ((y != 0) && (rand() / double(RAND_MAX)) < p_surface) 
										{
											int diff_num = DifferentNumber(i + x, j + y, k + z, core_id);
											if (diff_num == 0) Solid[cell_n] = core_id;
										}

										// grow in Z direction
										if ((z != 0) && (rand() / double(RAND_MAX)) < p_surface) 
										{
											int diff_num = DifferentNumber(i + x, j + y, k + z, core_id);
											if (diff_num == 0) Solid[cell_n] = core_id;
										}
										break;
									case 2:		// edge neighbor
										// if delete this case option, the boundary of the final porous media will be more smooth
										if ((rand() / double(RAND_MAX)) < p_edge) 
										{
											int diff_num = DifferentNumber(i + x, j + y, k + z, core_id);
											if (diff_num == 0) Solid[cell_n] = core_id;
										}
										break;
									case 3:		/// point neighbor
										// if delete this case option, the boundary of the final porous media will be more smooth
										if ((rand() / double(RAND_MAX)) < p_point)
										{
											int diff_num = DifferentNumber(i + x, j + y, k + z, core_id);
											if (diff_num == 0) Solid[cell_n] = core_id;
										}
										break;
									default:
										break;
									}
								}
							}
				}
			}

}

void Porous3D::output2tecplot(std::string filename)
{
	std::string out_buffer = "";
	ofstream outfile;
	std::string title = "3D porous media";
	outfile.open(filename);
	outfile << "TITLE = \"" << title << "\"" << endl;
	/// output solid position
	outfile << "VARIABLES = \"X\", \"Y\", \"Z\", \"value\" " << endl;
	outfile << "ZONE  t=\"solid\" I = " << this->NX << ", J = " << this->NY << ", K = " << this->NZ << ", F = point" << endl;
	double dx = 1.0 / this->NX;
	double dy = 1.0 / this->NY;
	double dz = 1.0 / this->NZ;
	for (int i = 0; i < this->NX; i++)
		for (int j = 0; j < this->NY; j++)
		{
			out_buffer = "";

			for (int k = 0; k < this->NZ; k++)
			{
				const int cell = CellIndex(i, j, k);
				out_buffer += std::to_string(i * dx) + "," + std::to_string(j * dy) + "," + std::to_string(k * dz) + "," + std::to_string(this->Solid[cell]) + "\n";
			}
			outfile << out_buffer;
		}

	outfile.close();
	std::cout << "Output completed" << endl;
}

double Porous3D::CalcPor()
{
	double t_temp = 0.0;
	int i;
#pragma omp parallel for private(i) reduction(+:t_temp)
	for (i = 0; i < NX * NY * NZ; i++)
	{
			t_temp += ((Solid[i]>0)? 1 : 0);
	}
			
	return (1 - t_temp / (NX * NY * NZ));
}

int Porous3D::CellIndex(int x, int y, int z)
{
	assert(x >= 0 && x < NX);
	assert(y >= 0 && y < NY);
	assert(z >= 0 && z < NZ);

	return x + y * NX + z * NX * NY;
}

int Porous3D::NearSolid(int i, int j, int k)
{
	int number = 0;
	for (int x = -1; x <= 1; x++)
		for (int y = -1; y <= 1; y++)
			for (int z = -1; z <= 1; z++)
			{
				if (!((i + x) < 0 || (j + y) < 0 || (k + z) < 0 || (i + x) > this->NX - 1 || (j + y) > this->NY - 1 || (k + z) > this->NZ - 1))
				{
					const int cell_n = CellIndex(i + x, j + y, k + z);
					number += Solid[cell_n];
				}
			}
	return number;
}

Porous3D::~Porous3D()
{
	delete[] this->Solid;
}

/// 0: fluid
/// 1: solid
/// 2: infected lattice
void Porous3D::InfectFluid(int i, int j, int k)
{
	for (int x = -1; x <= 1; x++)
		for (int y = -1; y <= 1; y++)
			for (int z = -1; z <= 1; z++)
			{
				if (!((i + x) < 0 || (j + y) < 0 || (k + z) < 0 || (i + x) > this->NX - 1 || (j + y) > this->NY - 1 || (k + z) > this->NZ - 1))
				{
					int cell = CellIndex(i + x, j + y, k + z);
					if (Solid[cell] == 0)
					{
						Solid[cell] = 2;
						InfectFluid(i + x, j + y, k + z);
					}
				}
			}
}

void Porous3D::DeleteDeadZone()
{
	/// infect fluid lattice from surface z = 0
	/// marked 2 as infected lattice
	/// 0: fluid
	/// 1: solid
	for (int i = 0; i < NX; i++)
		for (int j = 0; j < NY; j++)
		{
			int idx = CellIndex(i, j, 0);
			if (Solid[idx] == 0) Solid[idx] = 2;
			InfectFluid(i, j, 0);
		}

	/// modify fluid lattice in dead zone to solid lattice
#pragma omp parallel for
	for (int i = 0; i < NX * NY * NZ; i++)
		if (Solid[i] == 0) Solid[i] = 1;
	/// modify infected lattice to fluid lattice
#pragma omp parallel for
	for (int i = 0; i < NX * NY * NZ; i++)
		if (Solid[i] == 2) Solid[i] = 0;


	/// infect fluid lattice from surface z = NZ - 1
	/// marked 2 as infected lattice
	/// 0: fluid
	/// 1: solid
	for (int i = 0; i < NX; i++)
		for (int j = 0; j < NY; j++)
		{
			int idx = CellIndex(i, j, NZ - 1);
			if (Solid[idx] == 0) Solid[idx] = 2;
			InfectFluid(i, j, NZ - 1);
		}

	/// modify fluid lattice in dead zone to solid lattice
#pragma omp parallel for
	for (int i = 0; i < NX * NY * NZ; i++)
		if (Solid[i] == 0) Solid[i] = 1;
	/// modify infected lattice to fluid lattice
#pragma omp parallel for
	for (int i = 0; i < NX * NY * NZ; i++)
		if (Solid[i] == 2) Solid[i] = 0;

	
	/// infect fluid lattice from surface y = 0
	/// marked 2 as infected lattice
	/// 0: fluid
	/// 1: solid
	for (int i = 0; i < NX; i++)
		for (int k = 0; k < NZ; k++)
		{
			int idx = CellIndex(i, 0, k);
			if (Solid[idx] == 0) Solid[idx] = 2;
			InfectFluid(i, 0, k);
		}

	/// modify fluid lattice in dead zone to solid lattice
#pragma omp parallel for
	for (int i = 0; i < NX * NY * NZ; i++)
		if (Solid[i] == 0) Solid[i] = 1;
	/// modify infected lattice to fluid lattice
#pragma omp parallel for
	for (int i = 0; i < NX * NY * NZ; i++)
		if (Solid[i] == 2) Solid[i] = 0;

	/// infect fluid lattice from surface y = NY - 1
	/// marked 2 as infected lattice
	/// 0: fluid
	/// 1: solid
	for (int i = 0; i < NX; i++)
		for (int k = 0; k < NZ; k++)
		{
			int idx = CellIndex(i, NY - 1, k);
			if (Solid[idx] == 0) Solid[idx] = 2;
			InfectFluid(i, NY - 1, k);
		}

	/// modify fluid lattice in dead zone to solid lattice
#pragma omp parallel for
	for (int i = 0; i < NX * NY * NZ; i++)
		if (Solid[i] == 0) Solid[i] = 1;
	/// modify infected lattice to fluid lattice
#pragma omp parallel for
	for (int i = 0; i < NX * NY * NZ; i++)
		if (Solid[i] == 2) Solid[i] = 0;

	/// infect fluid lattice from surface x = 0
	/// marked 2 as infected lattice
	/// 0: fluid
	/// 1: solid
	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
		{
			int idx = CellIndex(0, j, k);
			if (Solid[idx] == 0) Solid[idx] = 2;
			InfectFluid(0, j, k);
		}

	/// modify fluid lattice in dead zone to solid lattice
#pragma omp parallel for
	for (int i = 0; i < NX * NY * NZ; i++)
		if (Solid[i] == 0) Solid[i] = 1;
	/// modify infected lattice to fluid lattice
#pragma omp parallel for
	for (int i = 0; i < NX * NY * NZ; i++)
		if (Solid[i] == 2) Solid[i] = 0;

	/// infect fluid lattice from surface x = NX - 1
	/// marked 2 as infected lattice
	/// 0: fluid
	/// 1: solid
	for (int j = 0; j < NY; j++)
		for (int k = 0; k < NZ; k++)
		{
			int idx = CellIndex(NX - 1, j, k);
			if (Solid[idx] == 0) Solid[idx] = 2;
			InfectFluid(NX - 1, j, k);
		}

	/// modify fluid lattice in dead zone to solid lattice
#pragma omp parallel for
	for (int i = 0; i < NX * NY * NZ; i++)
		if (Solid[i] == 0) Solid[i] = 1;
	/// modify infected lattice to fluid lattice
#pragma omp parallel for
	for (int i = 0; i < NX * NY * NZ; i++)
		if (Solid[i] == 2) Solid[i] = 0;
}

int Porous3D::DifferentNumber(int i, int j, int k, int core_id)
{
	int Num = 0;
	for (int x = -1; x <= 1; x++)
		for (int y = -1; y <= 1; y++)
			for (int z = -1; z <= 1; z++)
			{
				if (!((i + x) < 0 || (j + y) < 0 || (k + z) < 0 || (i + x) > this->NX - 1 || (j + y) > this->NY - 1 || (k + z) > this->NZ - 1))
				{
					int cell = CellIndex(i + x, j + y, k + z);
					if ((Solid[cell] != core_id) && (Solid[cell] != 0)) Num++;
				}
			}

	return Num;
}