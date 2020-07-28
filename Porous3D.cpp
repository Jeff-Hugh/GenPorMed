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

bool Porous3D::is_interior(const int i, const int j, const int k)
{
	const int cell = CellIndex(i, j, k);
	if (this->Solid[cell] == 1)
	{
		int sum_i = 0;
		for (int x = -1; x <= 1; x++)
			for (int y = -1; y <= 1; y++)
				for (int z = -1; z <= 1; z++)
				{
					if (!((i + x) < 0 || (j + y) < 0 || (k + z) < 0 || 
						(i + x) > this->NX - 1 || (j + y) > this->NY - 1 || 
						(k + z) > this->NZ - 1))
					{
						const int cell_n = CellIndex(i + x, j + y, k + z);
						sum_i =+ this->Solid[cell_n];
						if (this->Solid[cell_n] == 0) return false;
					}
					else
					{
						sum_i = +1;
					}

				}

		if (sum_i == 27)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}

Porous3D::Porous3D(int nx, int ny, int nz, int *s, double por = 0.6, double p = 0.05,
	double p_s = 0.2, double p_e = 0.0167, double p_p = 0.0021):NX(nx), NY(ny), NZ(nz), phi(por), p_cd(p), p_surface(p_s), p_edge(p_e), p_point(p_p)
{
	double phi_p;
	/// Generate random seed
	srand((unsigned)time(NULL));

	/// Initialization
	this->Solid = new int[NX*NY*NZ]{};

	/// Generate growth core
	for (int i = 0; i < this->NX; i++)
		for (int j = 0; j < this->NY; j++)
			for (int k = 0; k < this->NZ; k++)
			{
				const int cell = CellIndex(i, j, k);
				if ((rand() / double(RAND_MAX)) < p_cd) this->Solid[cell] = 1;
			}
				
	///  Generate process
	this->Grow_Times = 0;
	do
	{
		grow();
		
		if (Grow_Times % 5 == 0)
		{
			bool FixSmooth = false;
			do
			{
				FixSmooth = fixhole();
			} while (!FixSmooth);
		}

		if (Grow_Times % 5 == 0)
		{
			bool FixSmooth = false;
			do
			{
				FixSmooth = smooth();
			} while (!FixSmooth);
		}

		/// Calculate porosity
		phi_p = CalcPor();
		//std:: cout << "Porosity: " << phi_p << endl;
	} while (phi_p > phi);

# pragma omp parallel for
	for (int i = 0; i < this->NX; i++)
		for (int j = 0; j < this->NY; j++)
			for (int k = 0; k < this->NZ; k++)
			{
				const int cell = CellIndex(i, j, k);
				s[cell] = this->Solid[cell];
			}
			
	std::cout << "Grow completed, Porosity: " << phi_p << endl;
}

void Porous3D::grow()
{
	Grow_Times += 1;
	std::cout << "Grow Times: " << Grow_Times << endl;

	int* Solid_p = new int[NX * NY * NZ]{};
# pragma omp parallel for
	for (int i = 0; i < this->NX; i++)
		for (int j = 0; j < this->NY; j++)
			for (int k = 0; k < this->NZ; k++)
			{
				const int cell = CellIndex(i, j, k);
				Solid_p[cell] = Solid[cell];
			}
	// cout << "Copy compeleted" << endl;

# pragma omp parallel for 
	for (int i = 0; i < this->NX; i++)
		for (int j = 0; j < this->NY; j++)
			for (int k = 0; k < this->NZ; k++)
			{
				const int cell = CellIndex(i, j, k);
				if (this->Solid[cell] == 1)
				{
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
										Solid_p[cell_n] = 1;
										break;
									case 1:		// face neighbor
										if ((x != 0) && (rand() / double(RAND_MAX)) < p_surface) 
											Solid_p[cell_n] = 1;
										if ((y != 0) && (rand() / double(RAND_MAX)) < p_surface)
											Solid_p[cell_n] = 1;
										if ((z != 0) && (rand() / double(RAND_MAX)) < p_surface / 8)
											Solid_p[cell_n] = 1;
										break;
									case 2:		// edge neighbor
										if ((rand() / double(RAND_MAX)) < p_edge) Solid_p[cell_n] = 1;
										break;
									case 3:		/// point neighbor
										if ((rand() / double(RAND_MAX)) < p_point) Solid_p[cell_n] = 1;
										break;
									default:
										break;
									}
								}
							}
				}
			}

# pragma omp parallel for
	for (int i = 0; i < this->NX; i++)
		for (int j = 0; j < this->NY; j++)
			for (int k = 0; k < this->NZ; k++)
			{
				const int cell = CellIndex(i, j, k);
				int s_temp = this->Solid[cell] + Solid_p[cell];
				if (s_temp > 0) this->Solid[cell] = 1;
			}
	
	delete[] Solid_p;
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
	{
		out_buffer = "";

		for (int j = 0; j < this->NY; j++)
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
	int i, j, k;
#pragma omp parallel for private(i,j,k) reduction(+:t_temp)
	for (i = 0; i < NX; i++)
		for (j = 0; j < NY; j++)
			for (k = 0; k < NZ; k++)
			{
				const int cell = CellIndex(i, j, k);
				t_temp += Solid[cell];
			}
			
	return (1 - t_temp / (NX * NY * NZ));
}

bool Porous3D::fixhole()
{
	return true;
}

bool Porous3D::smooth()
{
	return true;
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
	
}