/*
\file GSGS.cpp
\reference Wang, M., Wang, J., Pan, N., & Chen, S. (2007).
Mesoscopic predictions of the effective thermal conductivity for microscale random porous media.
Physical Review E - Statistical, Nonlinear, and Soft Matter Physics, 75(3), 1¨C10.
https://doi.org/10.1103/PhysRevE.75.036702

zero for fluid
one for solid
*/

#include <iostream>
#include <stdlib.h> 
#include <time.h>
#include <string>
#include <fstream>
#include <vector>
#include <assert.h>
#include <omp.h>

#include "Porous2D.h"
#include "gnuplot-iostream.h"

using namespace std;

Porous2D::Porous2D(int nx, int ny, double por = 0.5, double p = 0.05, double p_z = 0.05, 
	double p_f = 0.0125):NX(nx), NY(ny), phi(por), p_cd(p), z(p_z), f(p_f)
{
}

void Porous2D::Generation(int* s)
{
	double phi_p;
	/// Generate random seed
	srand((unsigned)time(NULL));
	// std::cout << "Generate random seed completed" << endl;

	/// Initialization to zero fluid
	this->Solid = new int[NX * NY]{};

	// every core has different id
	int core_id = 1;

	/// Generate growth core
	for (int i = 0; i < this->NX; i++)
		for (int j = 0; j < this->NY; j++)
		{
			const int cell = this->CellIndex(i, j, this->NX, this->NY);
			if ( ((rand() / double(RAND_MAX)) < p_cd) && (NearSolid(i,j) == 0) )
			{
				this->Solid[cell] = core_id++;
			}
		}
	///  Generate process
	this->Grow_Times = 0;
	do
	{
		grow_alone();
		phi_p = CalcPor();
	}while (phi_p > phi);
	

	// delete CoreID
#pragma omp parallel for
	for (int i = 0; i < NX * NY; i++)
	{
		if (Solid[i] >= 1) Solid[i] = 1;
	}

	DeleteDeadZone();
	phi_p = CalcPor();

	phi = phi_p;
	std::cout << "Grow completed. Porosoty: " << phi << endl;


	//copy Solid (inner) to s (outsider)
#pragma omp parallel for
	for (int i = 0; i < NX * NY; i++) s[i] = Solid[i];

	delete[] this->Solid;
}

Porous2D::~Porous2D()
{
}

void Porous2D::output2tecplot(const int Total_M, const int Total_N, int* S, string filename)
{
	ofstream outfile;
	string title = "2D porous media";
	outfile.open(filename);

	outfile << "TITLE = \"" << title << "\"" << endl;

	/// output solid position
	outfile << "VARIABLES = \"X\", \"Y\",\"value\" " << endl;
	outfile << "ZONE  t = \"solid\" " << " I = " << Total_M << ", J = " << Total_N << " F = point" << endl;

	for (int j = 0; j < Total_N; j++)
		for (int i = 0; i < Total_M; i++)
		{
			const int cell = CellIndex(i, j, Total_M, Total_N);
			outfile << i << "," << j << ", " << S[cell] << endl;
		}

	outfile.close();
	std::cout << "Output completed for tecplot." << endl;
}


/// Grow in eight directions.
///
///*****    6    2    5   *****
///*****                  *****
///*****    3    C    1   *****
///*****                  *****
///*****    7    4    8   *****
///

void Porous2D::grow_alone()
{
	this->Grow_Times += 1;
	//std::cout << "Grow Times: " << this->Grow_Times << endl;

	// grow
	for (int i = 0; i < this->NX; i++)
	{
		for (int j = 0; j < this->NY; j++)
		{
			const int cell = this->CellIndex(i, j, this->NX, this->NY);

			if (i == 0)
			{
				if (j == 0)
				{
					if (this->Solid[cell] != 0)
					{
						grow_d1(i, j, Solid, Solid[cell]);
						grow_d2(i, j, Solid, Solid[cell]);
						//grow_d5(i, j, Solid, Solid[cell]);
					}
				}
				else if (j == this->NY - 1)
				{
					if (this->Solid[cell] != 0)
					{
						grow_d1(i, j, Solid, Solid[cell]);
						grow_d4(i, j, Solid, Solid[cell]);
						//grow_d8(i, j, Solid, Solid[cell]);
					}
				}
				else
				{
					if (this->Solid[cell] != 0)
					{
						grow_d1(i, j, Solid, Solid[cell]);
						grow_d2(i, j, Solid, Solid[cell]);
						grow_d4(i, j, Solid, Solid[cell]);
						//grow_d5(i, j, Solid, Solid[cell]);
						//grow_d8(i, j, Solid, Solid[cell]);
					}
				}
			}
			else if (i == this->NX - 1)
			{
				if (j == 0)
				{
					if (this->Solid[cell] != 0)
					{
						grow_d2(i, j, Solid, Solid[cell]);
						grow_d3(i, j, Solid, Solid[cell]);
						//grow_d6(i, j, Solid, Solid[cell]);
					}
				}
				else if (j == this->NY - 1)
				{
					if (this->Solid[cell] != 0)
					{
						grow_d3(i, j, Solid, Solid[cell]);
						grow_d4(i, j, Solid, Solid[cell]);
						//grow_d7(i, j, Solid, Solid[cell]);
					}
				}
				else
				{
					if (this->Solid[cell] != 0)
					{
						grow_d2(i, j, Solid, Solid[cell]);
						grow_d3(i, j, Solid, Solid[cell]);
						grow_d4(i, j, Solid, Solid[cell]);
						//grow_d6(i, j, Solid, Solid[cell]);
						//grow_d7(i, j, Solid, Solid[cell]);
					}
				}
			}
			else
			{
				if (j == 0)
				{
					if (this->Solid[cell] != 0)
					{
						grow_d1(i, j, Solid, Solid[cell]);
						grow_d2(i, j, Solid, Solid[cell]);
						grow_d3(i, j, Solid, Solid[cell]);
						//grow_d5(i, j, Solid, Solid[cell]);
						//grow_d6(i, j, Solid, Solid[cell]);
					}
				}
				else if (j == this->NY - 1)
				{
					if (this->Solid[cell] != 0)
					{
						grow_d1(i, j, Solid, Solid[cell]);
						grow_d3(i, j, Solid, Solid[cell]);
						grow_d4(i, j, Solid, Solid[cell]);
						//grow_d7(i, j, Solid, Solid[cell]);
						//grow_d8(i, j, Solid, Solid[cell]);
					}
				}
				else
				{
					if (this->Solid[cell] != 0)
					{
						grow_d1(i, j, Solid, Solid[cell]);
						grow_d2(i, j, Solid, Solid[cell]);
						grow_d3(i, j, Solid, Solid[cell]);
						grow_d4(i, j, Solid, Solid[cell]);
						//grow_d5(i, j, Solid, Solid[cell]);
						//grow_d6(i, j, Solid, Solid[cell]);
						//grow_d7(i, j, Solid, Solid[cell]);
						//grow_d8(i, j, Solid, Solid[cell]);
					}
				}
			}
		}
	}
		
}

bool Porous2D::fixhole()
{
	int n = 0;

#pragma omp parallel for
	for (int i = 0; i < this->NX; i++)
	{
		for (int j = 0; j < this->NY; j++)
		{
			const int cell = this->CellIndex(i, j, this->NX, this->NY);
			if (this->Solid[cell] == 0)
			{
				if (i == 0)
				{
					if (j == 0)
					{
						if (NearSolid(i, j) == 3) this->Solid[cell] = 1, n++;
					}
					else if (j == NY - 1)
					{
						if (NearSolid(i, j) == 3) this->Solid[cell] = 1, n++;
					}
					else
					{
						if (NearSolid(i, j) >= 4) this->Solid[cell] = 1, n++;
					}
				}
				else if (i == NX - 1)
				{
					if (j == 0)
					{
						if (NearSolid(i, j) == 3) this->Solid[cell] = 1, n++;
					}
					else if (j == NY - 1)
					{
						if (NearSolid(i, j) == 3) this->Solid[cell] = 1, n++;
					}
					else
					{
						if (NearSolid(i, j) >= 4) this->Solid[cell] = 1, n++;
					}
				}
				else
				{
					if (j == 0)
					{
						if (NearSolid(i, j) >= 4) this->Solid[cell] = 1, n++;
					}
					else if (j == NY - 1)
					{
						if (NearSolid(i, j) >= 4) this->Solid[cell] = 1, n++;
					}
					else
					{
						if (NearSolid(i, j) >= 7) this->Solid[cell] = 1, n++;
					}
				}
			}
		}
	}
	if (n == 0)
		return true;	/// fix over
	else
		return false;
}

int Porous2D::NearSolid(int i, int j)
{
	if (i == 0)
	{
		if (j == 0)
		{
			int tmp = 0;
			tmp += (Solid[this->CellIndex(i+1, j, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i+1, j+1, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i, j+1, this->NX, this->NY)] > 0) ? 1 : 0;
			return   tmp;
		}
		else if (j == NY - 1)
		{
			int tmp = 0;
			tmp += (Solid[this->CellIndex(i+1, j, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i+1, j-1, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i, j-1, this->NX, this->NY)] > 0) ? 1 : 0;
			return tmp;
		}
		else
		{
			int tmp = 0;
			tmp += (Solid[this->CellIndex(i+1, j, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i+1, j-1, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i, j-1, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i+1, j+1, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i, j+1, this->NX, this->NY)] > 0) ? 1 : 0;
			return tmp;
		}
	}
	else if (i == NX - 1)
	{
		if (j == 0)
		{
			int tmp = 0;
			tmp += (Solid[this->CellIndex(i-1, j, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i-1, j+1, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i, j+1, this->NX, this->NY)] > 0) ? 1 : 0;
			return tmp;
		}
		else if (j == NY - 1)
		{
			int tmp = 0;
			tmp += (Solid[this->CellIndex(i-1, j, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i-1, j-1, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i, j-1, this->NX, this->NY)] > 0) ? 1 : 0;
			return tmp;
		}
		else
		{
			int tmp = 0;
			tmp += (Solid[this->CellIndex(i-1, j, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i-1, j+1, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i, j+1, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i-1, j-1, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i, j-1, this->NX, this->NY)] > 0) ? 1 : 0;
			return tmp;
		}
	}
	else
	{
		if (j == 0)
		{
			int tmp = 0;
			tmp += (Solid[this->CellIndex(i-1, j, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i-1, j+1, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i, j+1, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i+1, j, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i+1, j+1, this->NX, this->NY)] > 0) ? 1 : 0;
			return  tmp;
		}
		else if (j == NY - 1)
		{
			int tmp = 0;
			tmp += (Solid[this->CellIndex(i-1, j, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i-1, j-1, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i, j-1, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i+1, j, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i+1, j-1, this->NX, this->NY)] > 0) ? 1 : 0;
			return  tmp;
		}
		else
		{
			int tmp = 0;
			tmp += (Solid[this->CellIndex(i-1, j, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i+1, j, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i-1, j+1, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i, j+1, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i+1, j+1, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i-1, j-1, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i, j-1, this->NX, this->NY)] > 0) ? 1 : 0;
			tmp += (Solid[this->CellIndex(i+1, j-1, this->NX, this->NY)] > 0) ? 1 : 0;
			return  tmp;
		}
	}
}

bool Porous2D::smooth()
{
	int n = 0;
#pragma omp parallel for
	for (int i = 0; i < this->NX; i++)
	{
		for (int j = 0; j < this->NY; j++)
		{
			const int cell = this->CellIndex(i, j, this->NX, this->NY);
			if (this->Solid[cell] == 1)
			{
				if (i == 0)
				{
					if (j == 0)
					{
						if (NearSolid(i, j) == 0) this->Solid[cell] = 0, n++;
					}
					else if (j == NY - 1)
					{
						if (NearSolid(i, j) == 0) this->Solid[cell] = 0, n++;
					}
					else
					{
						if (NearSolid(i, j) <= 1) this->Solid[cell] = 0, n++;
					}
				}
				else if (i == NX - 1)
				{
					if (j == 0)
					{
						if (NearSolid(i, j) == 0) this->Solid[cell] = 0, n++;
					}
					else if (j == NY - 1)
					{
						if (NearSolid(i, j) == 0) this->Solid[cell] = 0, n++;
					}
					else
					{
						if (NearSolid(i, j) <= 1) this->Solid[cell] = 0, n++;
					}
				}
				else
				{
					if (j == 0)
					{
						if (NearSolid(i, j) <= 1) this->Solid[cell] = 0, n++;
					}
					else if (j == NY - 1)
					{
						if (NearSolid(i, j) <= 1) this->Solid[cell] = 0, n++;
					}
					else
					{
						if (NearSolid(i, j) <= 2) this->Solid[cell] = 0, n++;
					}
				}
			}
		}
	}
	if (n == 0)
		return true;	/// smooth over
	else
		return false;
}

double Porous2D::CalcPor()
{
	double t_temp = 0;
	int i;
#pragma omp parallel for private(i) reduction(+:t_temp)
	for (i = 0; i < NX * NY; i++)
	{
			t_temp += ((Solid[i]>0)? 1 : 0);
	}
	
	return (1 - t_temp / (NX * NY));
}

int Porous2D::CellIndex(int x, int y, const int MM, const int NN)
{
	assert(x >= 0 && x < MM);
	assert(y >= 0 && y < NN);

	return x + y * MM;
}

void Porous2D::offset(const int Total_M, const int Total_N, const int offset_x, const int offset_y, int* S, const int M, const int N, int* Solid)
{
	if (Total_M <= (offset_x + M) || Total_N <= (offset_y + N))
	{
		std::cout << "No need to offset purous aera." << endl;
	}
	else
	{
#pragma omp parallel for
		for (int i = 0; i < M; i++)
			for (int j = 0; j < N; j++)
			{
				const int cell = CellIndex(i, j, M, N);
				const int Total_cell = CellIndex(offset_x + i, offset_y + j, Total_M, Total_N);
				if (Solid[cell] == 1)
					S[Total_cell] = 1;
				else
					S[Total_cell] = 0;
			}
	}
}

void Porous2D::output2png(const int Total_M, const int Total_N, int* S, std::string filename)
{
	Gnuplot gp;

	// create temp txt file for output
	ofstream outfile;
	string  temp_name = filename + ".txt";
	outfile.open(temp_name);
	for (int j = 0; j < Total_N; j++)
	{
		for (int i = 0; i < Total_M; i++)
		{
			const int cell = CellIndex(i, j, Total_M, Total_N);
			outfile << S[cell] << " ";
		}
		outfile << "\n";
	}
	outfile.close();

	gp << "set term png size 800,800\n";
	gp << "set output \"" << filename << "\"\n";
	gp << "set pm3d map\n";
	gp << "unset colorbox\n";
	gp << "unset xtics\n";
	gp << "unset ytics\n";
	gp << "splot '" << temp_name << "' matrix\n";

	//char* temp = temp_name.data();
	//if (remove(temp) != 0) std::cout << "Failure to remove temp file.";

#ifdef _WIN32
	// For Windows, prompt for a keystroke before the Gnuplot object goes out of scope so that
	// the gnuplot window doesn't get closed.
	//std::cout << "Press enter to exit." << std::endl;
	//std::cin.get();
#endif
	std::cout << "Output png file completed." << endl;
}

void Porous2D::grow_d1(int i, int j, int* s, int CoreID)
{
	if (((rand() / double(RAND_MAX)) < p_d1) && (i + 1 < this->NX))
	{
		if (s[CellIndex(i + 1, j, NX, NY)] == 0)
		{
			if (i + 1 == NX - 1)
				// if it is next to the boundary, set to solid 
				s[CellIndex(i + 1, j, NX, NY)] = CoreID;
			else
			{
				// check if exist solid that do not equal to CoreID around the cell(i+1, j)
				if (j == NY - 1)
				{
					if ((s[CellIndex(i + 2, j, NX, NY)] > 0) && (s[CellIndex(i + 2, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 2, j - 1, NX, NY)] > 0) && (s[CellIndex(i + 2, j - 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 1, j - 1, NX, NY)] > 0) && (s[CellIndex(i + 1, j - 1, NX, NY)] != CoreID))
						s[CellIndex(i + 1, j, NX, NY)] = 0;
					else
						s[CellIndex(i + 1, j, NX, NY)] = CoreID;
					
				}
				else if (j == 0)
				{
					if ((s[CellIndex(i + 2, j, NX, NY)] > 0) && (s[CellIndex(i + 2, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 2, j + 1, NX, NY)] > 0) && (s[CellIndex(i + 2, j + 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 1, j + 1, NX, NY)] > 0) && (s[CellIndex(i + 1, j + 1, NX, NY)] != CoreID))
						s[CellIndex(i + 1, j, NX, NY)] = 0;
					else
						s[CellIndex(i + 1, j, NX, NY)] = CoreID;
				}
				else
				{
					if ((s[CellIndex(i + 2, j, NX, NY)] > 0) && (s[CellIndex(i + 2, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 2, j - 1, NX, NY)] > 0) && (s[CellIndex(i + 2, j - 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 1, j - 1, NX, NY)] > 0) && (s[CellIndex(i + 1, j - 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 2, j + 1, NX, NY)] > 0) && (s[CellIndex(i + 2, j + 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 1, j + 1, NX, NY)] > 0) && (s[CellIndex(i + 1, j + 1, NX, NY)] != CoreID))
						s[CellIndex(i + 1, j, NX, NY)] = 0;
					else
						s[CellIndex(i + 1, j, NX, NY)] = CoreID;
				}
			}	
		}
	}
}

void Porous2D::grow_d2(int i, int j, int* s, int CoreID)
{
	if (((rand() / double(RAND_MAX)) < p_d2) && (j + 1 < this->NY))
	{
		if (s[CellIndex(i, j + 1, NX, NY)] == 0)
		{
			if (j + 2 == this->NY)
				// if it is next to the boundary, set to solid
				s[this->CellIndex(i, j + 1, this->NX, this->NY)] = CoreID;
			else
			{
				// check if exist solid that do not equal to CoreID around the cell(i+1, j)
				if (i == NX - 1)
				{
					if ((s[CellIndex(i, j + 2, NX, NY)] > 0) && (s[CellIndex(i, j + 2, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 1, j + 2, NX, NY)] > 0) && (s[CellIndex(i - 1, j + 2, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 1, j + 1, NX, NY)] > 0) && (s[CellIndex(i - 1, j + 1, NX, NY)] != CoreID))
						s[CellIndex(i, j + 1, NX, NY)] = 0;
					else
						s[CellIndex(i, j + 1, NX, NY)] = CoreID;
				}
				else if (i == 0)
				{
					if ((s[CellIndex(i, j + 2, NX, NY)] > 0) && (s[CellIndex(i, j + 2, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 1, j + 2, NX, NY)] > 0) && (s[CellIndex(i + 1, j + 2, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 1, j + 1, NX, NY)] > 0) && (s[CellIndex(i + 1, j + 1, NX, NY)] != CoreID))
						s[CellIndex(i, j + 1, NX, NY)] = 0;
					else
						s[CellIndex(i, j + 1, NX, NY)] = CoreID;
				}
				else
				{
					if ((s[CellIndex(i, j + 2, NX, NY)] > 0) && (s[CellIndex(i, j + 2, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 1, j + 2, NX, NY)] > 0) && (s[CellIndex(i - 1, j + 2, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 1, j + 1, NX, NY)] > 0) && (s[CellIndex(i - 1, j + 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 1, j + 2, NX, NY)] > 0) && (s[CellIndex(i + 1, j + 2, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 1, j + 1, NX, NY)] > 0) && (s[CellIndex(i + 1, j + 1, NX, NY)] != CoreID))
						s[CellIndex(i, j + 1, NX, NY)] = 0;
					else
						s[CellIndex(i, j + 1, NX, NY)] = CoreID;
				}
			}
		}
	}
}

void Porous2D::grow_d3(int i, int j, int* s, int CoreID)
{
	if (((rand() / double(RAND_MAX)) < p_d3) && (i > 0))
	{
		if (s[CellIndex(i - 1, j, NX, NY)] == 0)
		{
			if (i == 1)
				// if it is next to the boundary, set to solid 
				s[CellIndex(i - 1, j, NX, NY)] = CoreID;
			else
			{
				// check if exist solid that do not equal to CoreID around the cell(i+1, j)
				if (j == NY - 1)
				{
					if ((s[CellIndex(i - 2, j, NX, NY)] > 0) && (s[CellIndex(i - 2, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 2, j - 1, NX, NY)] > 0) && (s[CellIndex(i - 2, j - 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 1, j - 1, NX, NY)] > 0) && (s[CellIndex(i - 1, j - 1, NX, NY)] != CoreID))
						s[CellIndex(i - 1, j, NX, NY)] = 0;
					else
						s[CellIndex(i - 1, j, NX, NY)] = CoreID;
				}
				else if (j == 0)
				{
					if ((s[CellIndex(i - 2, j, NX, NY)] > 0) && (s[CellIndex(i - 2, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 2, j + 1, NX, NY)] > 0) && (s[CellIndex(i - 2, j + 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 1, j + 1, NX, NY)] > 0) && (s[CellIndex(i - 1, j + 1, NX, NY)] != CoreID))
						s[CellIndex(i - 1, j, NX, NY)] = 0;
					else
						s[CellIndex(i - 1, j, NX, NY)] = CoreID;
				}
				else
				{
					if ((s[CellIndex(i - 2, j, NX, NY)] > 0) && (s[CellIndex(i - 2, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 2, j - 1, NX, NY)] > 0) && (s[CellIndex(i - 2, j - 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 1, j - 1, NX, NY)] > 0) && (s[CellIndex(i - 1, j - 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 2, j + 1, NX, NY)] > 0) && (s[CellIndex(i - 2, j + 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 1, j + 1, NX, NY)] > 0) && (s[CellIndex(i - 1, j + 1, NX, NY)] != CoreID))
						s[CellIndex(i - 1, j, NX, NY)] = 0;
					else
						s[CellIndex(i - 1, j, NX, NY)] = CoreID;
				}
			}
		}
	}
}

void Porous2D::grow_d4(int i, int j, int* s, int CoreID)
{
	if (((rand() / double(RAND_MAX)) < p_d4) && (j > 0))
	{
		if (s[CellIndex(i, j - 1, NX, NY)] == 0)
		{
			if (j == 1)
				// if it is next to the boundary, set to solid
				s[this->CellIndex(i, j - 1, this->NX, this->NY)] = CoreID;
			else
			{
				// check if exist solid that do not equal to CoreID around the cell(i+1, j)
				if (i == NX - 1)
				{
					if ((s[CellIndex(i, j - 2, NX, NY)] > 0) && (s[CellIndex(i, j - 2, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 1, j - 2, NX, NY)] > 0) && (s[CellIndex(i - 1, j - 2, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 1, j - 1, NX, NY)] > 0) && (s[CellIndex(i - 1, j - 1, NX, NY)] != CoreID))
						s[CellIndex(i, j - 1, NX, NY)] = 0;
					else
						s[CellIndex(i, j - 1, NX, NY)] = CoreID;
				}
				else if (i == 0)
				{
					if ((s[CellIndex(i, j - 2, NX, NY)] > 0) && (s[CellIndex(i, j - 2, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 1, j - 2, NX, NY)] > 0) && (s[CellIndex(i + 1, j - 2, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 1, j - 1, NX, NY)] > 0) && (s[CellIndex(i + 1, j - 1, NX, NY)] != CoreID))
						s[CellIndex(i, j - 1, NX, NY)] = 0;
					else
						s[CellIndex(i, j - 1, NX, NY)] = CoreID;
				}
				else
				{
					if ((s[CellIndex(i, j +- 2, NX, NY)] > 0) && (s[CellIndex(i, j - 2, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 1, j - 2, NX, NY)] > 0) && (s[CellIndex(i - 1, j - 2, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 1, j - 1, NX, NY)] > 0) && (s[CellIndex(i - 1, j - 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 1, j - 2, NX, NY)] > 0) && (s[CellIndex(i + 1, j - 2, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 1, j - 1, NX, NY)] > 0) && (s[CellIndex(i + 1, j - 1, NX, NY)] != CoreID))
						s[CellIndex(i, j - 1, NX, NY)] = 0;
					else
						s[CellIndex(i, j - 1, NX, NY)] = CoreID;
				}
			}
		}
	}
}

void Porous2D::grow_d5(int i, int j, int* s, int CoreID)
{
	if (((rand() / double(RAND_MAX)) < p_d5) && (j + 1 < this->NY) && (i + 1 < this->NX))
	{
		if (s[CellIndex(i + 1, j + 1, NX, NY)] == 0)
		{
			if ((i + 1 == NX - 1) && (j + 1 == NY - 1))
				// if it is next to the boundary, set to solid 
				s[CellIndex(i + 1, j + 1, NX, NY)] = CoreID;
			else
			{
				// check if exist solid that do not equal to CoreID around the cell(i+1, j)
				if ((j + 1 == NY - 1) && (i + 1 < NX - 1))
				{
					if ((s[CellIndex(i, j + 1, NX, NY)] > 0) && (s[CellIndex(i, j + 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 1, j, NX, NY)] > 0) && (s[CellIndex(i + 1, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 2, j, NX, NY)] > 0) && (s[CellIndex(i + 2, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 2, j + 1, NX, NY)] > 0) && (s[CellIndex(i + 2, j + 1, NX, NY)] != CoreID))
						s[CellIndex(i + 1, j + 1, NX, NY)]  = 0;
					else
						s[CellIndex(i + 1, j + 1, NX, NY)]  = CoreID;
					
				}
				else if ((j + 1 < NY - 1) && (i + 1 == NX - 1))
				{
					if ((s[CellIndex(i, j + 1, NX, NY)] > 0) && (s[CellIndex(i, j + 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 1, j, NX, NY)] > 0) && (s[CellIndex(i + 1, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i, j + 2, NX, NY)] > 0) && (s[CellIndex(i, j + 2, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 1, j + 2, NX, NY)] > 0) && (s[CellIndex(i + 1, j + 2, NX, NY)] != CoreID))
						s[CellIndex(i + 1, j + 1, NX, NY)]  = 0;
					else
						s[CellIndex(i + 1, j + 1, NX, NY)]  = CoreID;
				}
				else
				{
					if ((s[CellIndex(i, j + 1, NX, NY)] > 0) && (s[CellIndex(i, j + 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 1, j, NX, NY)] > 0) && (s[CellIndex(i + 1, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 2, j, NX, NY)] > 0) && (s[CellIndex(i + 2, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 2, j + 1, NX, NY)] > 0) && (s[CellIndex(i + 2, j + 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i, j + 2, NX, NY)] > 0) && (s[CellIndex(i, j + 2, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 1, j + 2, NX, NY)] > 0) && (s[CellIndex(i + 1, j + 2, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 2, j + 2, NX, NY)] > 0) && (s[CellIndex(i + 2, j + 2, NX, NY)] != CoreID) )
						s[CellIndex(i + 1, j + 1, NX, NY)]  = 0;
					else
						s[CellIndex(i + 1, j + 1, NX, NY)]  = CoreID;
				}
			}	
		}
	}
}

void Porous2D::grow_d6(int i, int j, int* s, int CoreID)
{
	if (((rand() / double(RAND_MAX)) < p_d6) && (j + 1 < this->NY) && (i > 0))
	{
		if (s[CellIndex(i - 1, j + 1, NX, NY)] == 0)
		{
			if ((i - 1 == 0) && (j + 1 == NY - 1))
				// if it is next to the boundary, set to solid 
				s[CellIndex(i - 1, j + 1, NX, NY)] = CoreID;
			else
			{
				// check if exist solid that do not equal to CoreID around the cell(i+1, j)
				if ((j + 1 == NY - 1) && (i - 1 > 0)) 
				{
					if ((s[CellIndex(i, j + 1, NX, NY)] > 0) && (s[CellIndex(i, j + 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 1, j, NX, NY)] > 0) && (s[CellIndex(i - 1, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 2, j, NX, NY)] > 0) && (s[CellIndex(i - 2, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 2, j + 1, NX, NY)] > 0) && (s[CellIndex(i - 2, j + 1, NX, NY)] != CoreID))
						s[CellIndex(i - 1, j + 1, NX, NY)] = 0;
					else
						s[CellIndex(i - 1, j + 1, NX, NY)] = CoreID;
					
				}
				else if ((j + 1 < NY - 1) && (i == 1))
				{
					if ((s[CellIndex(i, j + 1, NX, NY)] > 0) && (s[CellIndex(i, j + 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 1, j, NX, NY)] > 0) && (s[CellIndex(i - 1, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i, j + 2, NX, NY)] > 0) && (s[CellIndex(i, j + 2, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 1, j + 2, NX, NY)] > 0) && (s[CellIndex(i - 1, j + 2, NX, NY)] != CoreID))
						s[CellIndex(i - 1, j + 1, NX, NY)] = 0;
					else
						s[CellIndex(i - 1, j + 1, NX, NY)] = CoreID;
				}
				else
				{
					if ((s[CellIndex(i, j + 1, NX, NY)] > 0) && (s[CellIndex(i, j + 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 1, j, NX, NY)] > 0) && (s[CellIndex(i - 1, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 2, j, NX, NY)] > 0) && (s[CellIndex(i - 2, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 2, j + 1, NX, NY)] > 0) && (s[CellIndex(i - 2, j + 1, NX, NY)] != CoreID) 
						|| (s[CellIndex(i, j + 2, NX, NY)] > 0) && (s[CellIndex(i, j + 2, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 1, j + 2, NX, NY)] > 0) && (s[CellIndex(i - 1, j + 2, NX, NY)] != CoreID) 
						|| (s[CellIndex(i - 2, j + 2, NX, NY)] > 0) && (s[CellIndex(i - 2, j + 2, NX, NY)] != CoreID))
						s[CellIndex(i - 1, j + 1, NX, NY)] = 0;
					else
						s[CellIndex(i - 1, j + 1, NX, NY)] = CoreID;
				}
			}	
		}
	}
}

void Porous2D::grow_d7(int i, int j, int* s, int CoreID)
{
	if (((rand() / double(RAND_MAX)) < p_d7) && (j > 0) && (i > 0))
	{
		if (s[CellIndex(i - 1, j - 1, NX, NY)] == 0)
		{
			if ((i - 1 == 0) && (j - 1 == 0))
				// if it is next to the boundary, set to solid 
				s[CellIndex(i - 1, j - 1, NX, NY)] = CoreID;
			else
			{
				// check if exist solid that do not equal to CoreID around the cell(i+1, j)
				if ((j - 1 == 0) && (i - 1 > 0)) 
				{
					if ((s[CellIndex(i, j - 1, NX, NY)] > 0) && (s[CellIndex(i, j - 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 1, j, NX, NY)] > 0) && (s[CellIndex(i - 1, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 2, j, NX, NY)] > 0) && (s[CellIndex(i - 2, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 2, j - 1, NX, NY)] > 0) && (s[CellIndex(i - 2, j - 1, NX, NY)] != CoreID))
						s[CellIndex(i - 1, j - 1, NX, NY)] = 0;
					else
						s[CellIndex(i - 1, j - 1, NX, NY)] = CoreID;
					
				}
				else if ((j - 1 > 0) && (i == 1))
				{
					if ((s[CellIndex(i, j - 1, NX, NY)] > 0) && (s[CellIndex(i, j - 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 1, j, NX, NY)] > 0) && (s[CellIndex(i - 1, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i, j - 2, NX, NY)] > 0) && (s[CellIndex(i, j - 2, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 1, j - 2, NX, NY)] > 0) && (s[CellIndex(i - 1, j - 2, NX, NY)] != CoreID))
						s[CellIndex(i - 1, j - 1, NX, NY)] = 0;
					else
						s[CellIndex(i - 1, j - 1, NX, NY)] = CoreID;
				}
				else
				{
					if ((s[CellIndex(i, j - 1, NX, NY)] > 0) && (s[CellIndex(i, j - 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 1, j, NX, NY)] > 0) && (s[CellIndex(i - 1, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 2, j, NX, NY)] > 0) && (s[CellIndex(i - 2, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 2, j - 1, NX, NY)] > 0) && (s[CellIndex(i - 2, j - 1, NX, NY)] != CoreID) 
						|| (s[CellIndex(i, j - 2, NX, NY)] > 0) && (s[CellIndex(i, j - 2, NX, NY)] != CoreID)
						|| (s[CellIndex(i - 1, j - 2, NX, NY)] > 0) && (s[CellIndex(i - 1, j - 2, NX, NY)] != CoreID) 
						|| (s[CellIndex(i - 2, j - 2, NX, NY)] > 0) && (s[CellIndex(i - 2, j - 2, NX, NY)] != CoreID))
						s[CellIndex(i - 1, j - 1, NX, NY)] = 0;
					else
						s[CellIndex(i - 1, j - 1, NX, NY)] = CoreID;
				}
			}	
		}
	}
}

void Porous2D::grow_d8(int i, int j, int* s, int CoreID)
{
	if (((rand() / double(RAND_MAX)) < p_d8) && (j > 0) && (i + 1 < this->NX))
	{
		if (s[CellIndex(i + 1, j - 1, NX, NY)] == 0)
		{
			if ((i + 1 == NX - 1) && (j - 1 == 0))
				// if it is next to the boundary, set to solid 
				s[CellIndex(i + 1, j - 1, NX, NY)] = CoreID;
			else
			{
				// check if exist solid that do not equal to CoreID around the cell(i+1, j)
				if ((j - 1 == 0) && (i + 1 < NX - 1)) 
				{
					if ((s[CellIndex(i, j - 1, NX, NY)] > 0) && (s[CellIndex(i, j - 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 1, j, NX, NY)] > 0) && (s[CellIndex(i + 1, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 2, j, NX, NY)] > 0) && (s[CellIndex(i + 2, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 2, j - 1, NX, NY)] > 0) && (s[CellIndex(i + 2, j - 1, NX, NY)] != CoreID))
						s[CellIndex(i + 1, j - 1, NX, NY)] = 0;
					else
						s[CellIndex(i + 1, j - 1, NX, NY)] = CoreID;
					
				}
				else if ((j - 1 > 0) && (i + 1 == NX - 1))
				{
					if ((s[CellIndex(i, j - 1, NX, NY)] > 0) && (s[CellIndex(i, j - 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 1, j, NX, NY)] > 0) && (s[CellIndex(i + 1, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i, j - 2, NX, NY)] > 0) && (s[CellIndex(i, j - 2, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 1, j - 2, NX, NY)] > 0) && (s[CellIndex(i + 1, j - 2, NX, NY)] != CoreID))
						s[CellIndex(i + 1, j - 1, NX, NY)] = 0;
					else
						s[CellIndex(i + 1, j - 1, NX, NY)] = CoreID;
				}
				else
				{
					if ((s[CellIndex(i, j - 1, NX, NY)] > 0) && (s[CellIndex(i, j - 1, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 1, j, NX, NY)] > 0) && (s[CellIndex(i + 1, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 2, j, NX, NY)] > 0) && (s[CellIndex(i + 2, j, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 2, j - 1, NX, NY)] > 0) && (s[CellIndex(i + 2, j - 1, NX, NY)] != CoreID) 
						|| (s[CellIndex(i, j - 2, NX, NY)] > 0) && (s[CellIndex(i, j - 2, NX, NY)] != CoreID)
						|| (s[CellIndex(i + 1, j - 2, NX, NY)] > 0) && (s[CellIndex(i + 1, j - 2, NX, NY)] != CoreID) 
						|| (s[CellIndex(i + 2, j - 2, NX, NY)] > 0) && (s[CellIndex(i + 2, j - 2, NX, NY)] != CoreID))
						s[CellIndex(i + 1, j - 1, NX, NY)] = 0;
					else
						s[CellIndex(i + 1, j - 1, NX, NY)] = CoreID;
				}
			}	
		}
	}
}

void Porous2D::DeleteDeadZone()
{
	/// infect fluid lattice from edge y = 0
	/// marked 2 as infected lattice
	for (int i = 0; i < NX; i++)
	{
		int idx = CellIndex(i, 0, NX, NY);
		if (Solid[idx] == 0) Solid[idx] = 2;
		InfectFluid(i, 0);
	}

#pragma omp parallel for
	for (int i = 0; i < NX * NY; i++)
		if (Solid[i] == 0) Solid[i] = 1;

#pragma omp parallel for
	for (int i = 0; i < NX * NY; i++)
		if (Solid[i] == 2) Solid[i] = 0;
//////////////////////////////////////////////////////////////
	/// infect fluid lattice from edge x = 0
	/// marked 2 as infected lattice
	for (int j = 0; j < NY; j++)
	{
		int idx = CellIndex(0, j, NX, NY);
		if (Solid[idx] == 0) Solid[idx] = 2;
		InfectFluid(0, j);
	}

#pragma omp parallel for
	for (int i = 0; i < NX * NY; i++)
		if (Solid[i] == 0) Solid[i] = 1;

#pragma omp parallel for
	for (int i = 0; i < NX * NY; i++)
		if (Solid[i] == 2) Solid[i] = 0;
//////////////////////////////////////////////////////////////
	/// infect fluid lattice from y = Ny - 1
	/// marked 2 as infected lattice
	for (int i = 0; i < NX; i++)
	{
		int idx = CellIndex(i, NY - 1, NX, NY);
		if (Solid[idx] == 0) Solid[idx] = 2;
		InfectFluid(i, NY - 1);
	}

#pragma omp parallel for
	for (int i = 0; i < NX * NY; i++)
		if (Solid[i] == 0) Solid[i] = 1;

#pragma omp parallel for
	for (int i = 0; i < NX * NY; i++)
		if (Solid[i] == 2) Solid[i] = 0;
//////////////////////////////////////////////////////////////
		/// infect fluid lattice from x = NX - 1
		/// marked 2 as infected lattice
	for (int j = 0; j < NY; j++)
	{
		int idx = CellIndex(NX - 1, j, NX, NY);
		if (Solid[idx] == 0) Solid[idx] = 2;
		InfectFluid(NX - 1, j);
	}

#pragma omp parallel for
	for (int i = 0; i < NX * NY; i++)
		if (Solid[i] == 0) Solid[i] = 1;

#pragma omp parallel for
	for (int i = 0; i < NX * NY; i++)
		if (Solid[i] == 2) Solid[i] = 0;
}

/// 0: fluid
/// 1: solid
/// 2: infected lattice
void Porous2D::InfectFluid(int i, int j)
{
	if (i == 0)
	{
		if (j == 0)
		{
			int idx_up = CellIndex(i + 1, j, NX, NY);
			int idx_right = CellIndex(i, j + 1, NX, NY);
			if (Solid[idx_up] == 0)
			{
				Solid[idx_up] = 2;
				InfectFluid(i + 1, j);
			}
			if (Solid[idx_right] == 0)
			{
				Solid[idx_right] = 2;
				InfectFluid(i, j + 1);
			}
		}
		else if (j == NY - 1)
		{
			int idx_up = CellIndex(i + 1, j, NX, NY);
			int idx_left = CellIndex(i, j - 1, NX, NY);
			if (Solid[idx_up] == 0)
			{
				Solid[idx_up] = 2;
				InfectFluid(i + 1, j);
			}
			if (Solid[idx_left] == 0)
			{
				Solid[idx_left] = 2;
				InfectFluid(i, j - 1);
			}
		}
		else
		{
			int idx_up = CellIndex(i + 1, j, NX, NY);
			int idx_right = CellIndex(i, j + 1, NX, NY);
			int idx_left = CellIndex(i, j - 1, NX, NY);
			if (Solid[idx_up] == 0)
			{
				Solid[idx_up] = 2;
				InfectFluid(i + 1, j);
			}
			if (Solid[idx_right] == 0)
			{
				Solid[idx_right] = 2;
				InfectFluid(i, j + 1);
			}
			if (Solid[idx_left] == 0)
			{
				Solid[idx_left] = 2;
				InfectFluid(i, j - 1);
			}
		}
	}
	else if (i == NX - 1)
	{
		if (j == 0)
		{
			int idx_down = CellIndex(i - 1, j, NX, NY);
			int idx_right = CellIndex(i, j + 1, NX, NY);
			if (Solid[idx_down] == 0)
			{
				Solid[idx_down] = 2;
				InfectFluid(i - 1, j);
			}
			if (Solid[idx_right] == 0)
			{
				Solid[idx_right] = 2;
				InfectFluid(i, j + 1);
			}
		}
		else if (j == NY - 1)
		{
			int idx_down = CellIndex(i - 1, j, NX, NY);
			int idx_left = CellIndex(i, j - 1, NX, NY);
			if (Solid[idx_down] == 0)
			{
				Solid[idx_down] = 2;
				InfectFluid(i - 1, j);
			}
			if (Solid[idx_left] == 0)
			{
				Solid[idx_left] = 2;
				InfectFluid(i, j - 1);
			}
		}
		else
		{
			int idx_down = CellIndex(i - 1, j, NX, NY);
			int idx_right = CellIndex(i, j + 1, NX, NY);
			int idx_left = CellIndex(i, j - 1, NX, NY);
			if (Solid[idx_down] == 0)
			{
				Solid[idx_down] = 2;
				InfectFluid(i - 1, j);
			}
			if (Solid[idx_right] == 0)
			{
				Solid[idx_right] = 2;
				InfectFluid(i, j + 1);
			}
			if (Solid[idx_left] == 0)
			{
				Solid[idx_left] = 2;
				InfectFluid(i, j - 1);
			}
		}
	}
	else
	{
		if (j == 0)
		{
			int idx_down = CellIndex(i - 1, j, NX, NY);
			int idx_right = CellIndex(i, j + 1, NX, NY);
			int idx_up = CellIndex(i + 1, j, NX, NY);
			if (Solid[idx_up] == 0)
			{
				Solid[idx_up] = 2;
				InfectFluid(i + 1, j);
			}
			if (Solid[idx_down] == 0)
			{
				Solid[idx_down] = 2;
				InfectFluid(i - 1, j);
			}
			if (Solid[idx_right] == 0)
			{
				Solid[idx_right] = 2;
				InfectFluid(i, j + 1);
			}
		}
		else if (j == NY - 1)
		{
			int idx_down = CellIndex(i - 1, j, NX, NY);
			int idx_left = CellIndex(i, j - 1, NX, NY);
			int idx_up = CellIndex(i + 1, j, NX, NY);
			if (Solid[idx_up] == 0)
			{
				Solid[idx_up] = 2;
				InfectFluid(i + 1, j);
			}
			if (Solid[idx_down] == 0)
			{
				Solid[idx_down] = 2;
				InfectFluid(i - 1, j);
			}
			if (Solid[idx_left] == 0)
			{
				Solid[idx_left] = 2;
				InfectFluid(i, j - 1);
			}
		}
		else
		{
			int idx_down = CellIndex(i - 1, j, NX, NY);
			int idx_left = CellIndex(i, j - 1, NX, NY);
			int idx_up = CellIndex(i + 1, j, NX, NY);
			int idx_right = CellIndex(i, j + 1, NX, NY);
			if (Solid[idx_up] == 0)
			{
				Solid[idx_up] = 2;
				InfectFluid(i + 1, j);
			}
			if (Solid[idx_down] == 0)
			{
				Solid[idx_down] = 2;
				InfectFluid(i - 1, j);
			}
			if (Solid[idx_left] == 0)
			{
				Solid[idx_left] = 2;
				InfectFluid(i, j - 1);
			}
			if (Solid[idx_right] == 0)
			{
				Solid[idx_right] = 2;
				InfectFluid(i, j + 1);
			}

		}
	}
}
