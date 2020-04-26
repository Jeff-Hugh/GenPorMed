/// \file GSGS.cpp
/// \reference Wang, M., Wang, J., Pan, N., & Chen, S. (2007). 
/// Mesoscopic predictions of the effective thermal conductivity for microscale random porous media. 
/// Physical Review E - Statistical, Nonlinear, and Soft Matter Physics, 75(3), 1¨C10.
/// https://doi.org/10.1103/PhysRevE.75.036702
/// 
/// 

#include <iostream>
#include <stdlib.h> 
#include <time.h>
#include <string>
#include <fstream>
#include <assert.h>
#include <omp.h>

using namespace std;

const float p_cd = 0.05;		/// core distribution probability
const float z = 0.05, f = 0.0125;
const float p_d1 = z, p_d2 = z, p_d3 = z, p_d4 = z;
const float p_d5 = f, p_d6 = f, p_d7 = f, p_d8 = f;

int Grow_Times = 0;

int CellIndex(int x, int y, const int MM, const int NN)
{
	assert(x >= 0 && x < MM);
	assert(y >= 0 && y < NN);

	return x + y * MM;
}


/// \brief Grow in eight directions.
///
/// Grow directions
///*****    6    2    5   *****
///*****                  *****
///*****    3    C    1   *****
///*****                  *****
///*****    7    4    8   *****
///
void grow(const int M, const int N, int * Solid)
{
	Grow_Times += 1;
	cout << "Grow Times: " << Grow_Times << endl;

	int* Solid_p = new int[M * N]{};
#pragma omp parallel for
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
		{
			int cell = CellIndex(i,j,M,N);
			Solid_p[cell] = Solid[cell];
		}
			
	cout << "Copy compeleted" << endl;

	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			const int cell = CellIndex(i,j,M,N);
			if (i == 0)
			{
				if (j == 0)
				{
					if (Solid[cell] == 1)
					{
						if ((rand() / double(RAND_MAX)) < p_d1) Solid_p[CellIndex(i+1,j,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d2) Solid_p[CellIndex(i,j+1,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d5) Solid_p[CellIndex(i+1,j+1,M,N)] = 1;
					}
				}
				else if (j == N - 1)
				{
					if (Solid[cell] == 1)
					{
						if ((rand() / double(RAND_MAX)) < p_d1) Solid_p[CellIndex(i+1,j,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d4) Solid_p[CellIndex(i,j-1,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d8) Solid_p[CellIndex(i+1,j-1,M,N)] = 1;
					}
				}
				else
				{
					if (Solid[cell]== 1)
					{
						if ((rand() / double(RAND_MAX)) < p_d1) Solid_p[CellIndex(i+1,j,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d2) Solid_p[CellIndex(i,j+1,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d4) Solid_p[CellIndex(i,j-1,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d5) Solid_p[CellIndex(i+1,j+1,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d8) Solid_p[CellIndex(i+1,j-1,M,N)] = 1;
					}
				}
			}
			else if (i == N - 1)
			{
				if (j == 0)
				{
					if (Solid[cell] == 1)
					{
						if ((rand() / double(RAND_MAX)) < p_d2) Solid_p[CellIndex(i,j+1,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d3) Solid_p[CellIndex(i-1,j,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d6) Solid_p[CellIndex(i-1,j+1,M,N)] = 1;
					}
				}
				else if (j == N - 1)
				{
					if (Solid[cell] == 1)
					{
						if ((rand() / double(RAND_MAX)) < p_d3) Solid_p[CellIndex(i-1,j,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d4) Solid_p[CellIndex(i,j-1,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d7) Solid_p[CellIndex(i-1,j-1,M,N)] = 1;
					}
				}
				else
				{
					if (Solid[cell] == 1)
					{
						if ((rand() / double(RAND_MAX)) < p_d2) Solid_p[CellIndex(i,j+1,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d3) Solid_p[CellIndex(i-1,j,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d4) Solid_p[CellIndex(i,j-1,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d6) Solid_p[CellIndex(i-1,j+1,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d7) Solid_p[CellIndex(i-1,j-1,M,N)] = 1;
					}
				}
			}
			else
			{
				if (j == 0)
				{
					if (Solid[cell] == 1)
					{
						if ((rand() / double(RAND_MAX)) < p_d1) Solid_p[CellIndex(i+1,j,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d2) Solid_p[CellIndex(i,j+1,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d3) Solid_p[CellIndex(i-1,j,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d5) Solid_p[CellIndex(i+1,j+1,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d6) Solid_p[CellIndex(i-1,j+1,M,N)] = 1;
					}
				}
				else if (j == N - 1)
				{
					if (Solid[cell] == 1)
					{
						if ((rand() / double(RAND_MAX)) < p_d1) Solid_p[CellIndex(i+1,j,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d3) Solid_p[CellIndex(i-1,j,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d4) Solid_p[CellIndex(i,j-1,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d7) Solid_p[CellIndex(i-1,j-1,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d8) Solid_p[CellIndex(i+1,j-1,M,N)] = 1;
					}
				}
				else
				{
					if (Solid[cell] == 1)
					{
						if ((rand() / double(RAND_MAX)) < p_d1) Solid_p[CellIndex(i+1,j,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d2) Solid_p[CellIndex(i,j+1,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d3) Solid_p[CellIndex(i-1,j,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d4) Solid_p[CellIndex(i,j-1,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d5) Solid_p[CellIndex(i+1,j+1,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d6) Solid_p[CellIndex(i-1,j+1,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d7) Solid_p[CellIndex(i-1,j-1,M,N)] = 1;
						if ((rand() / double(RAND_MAX)) < p_d8) Solid_p[CellIndex(i+1,j-1,M,N)] = 1;
					}
				}
			}
		}
	}
	cout << "Neighbor generated" << endl;

#pragma omp parallel for
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
		{
			const int cell = CellIndex(i,j,M,N);
			int s_temp = Solid[cell] + Solid_p[cell];
			if (s_temp > 0) Solid[cell] = 1;
		}
	cout << "One time grow compeleted" << endl;
}

void output2tecplot(const int Total_M, const int Total_N, int *S, string filename)
{
	ofstream outfile;
	string title = "2D porous media";
	outfile.open(filename);

	outfile << "TITLE = \"" << title << "\"" << endl;

	/// output solid position
	outfile << "VARIABLES = \"X\", \"Y\",\"value\" " << endl;
	outfile << "ZONE  t = \"solid\" " << " I = " << Total_M << ", J = " << Total_N << " F = point" << endl;
	outfile  << endl;

	for (int j = 0; j < Total_N; j++)
		for (int i = 0; i < Total_M; i++)
		{
			const int cell = CellIndex(i,j,Total_M, Total_N);
			outfile <<  i << "," << j << ", " << S[cell] << endl;
		}
	
	outfile.close();
	cout << "Output completed" << endl;
}

void output2gnuplot(const int Total_M, const int Total_N, int *S, string filename)
{
	string out_buffer = "";
	ofstream outfile;

	outfile.open(filename);
	for (int i = 0; i < Total_M; i++)
		for (int j = 0; j < Total_N; j++)
		{
			const int cell = CellIndex(i,j,Total_M, Total_N);
			out_buffer += std::to_string(S[cell]) + " ";
		}
		out_buffer += "\n";
	outfile << out_buffer;
	outfile.close();
	cout << "Output completed" << endl;
}

void offset(const int Total_M, const int Total_N, const int offset_x, const int offset_y, int *S, const int M, const int N, int * Solid)
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
				const int cell = CellIndex(i,j,M,N);
				const int Total_cell = CellIndex(offset_x + i, offset_y + j, Total_M, Total_N);
				if (Solid[cell] == 1)
					S[Total_cell] = 1;
				else
					S[Total_cell] = 0;
			}	
	}
}


int main()
{
	const int M = 100;
	const int N = 100;
	const float phi = 0.7;			/// target porosity
	float phi_p;

	/// initialize to 0, all fluid
	int* Solid = new int[M * N]{};

	/// Generate random seed
	srand((unsigned)time(NULL));
	cout << "Generate random seed completed" <<endl;

	/// Generate growth core
	#pragma omp paralle for
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
		{
			const int cell = CellIndex(i,j,M,N);
			if ((rand() / double(RAND_MAX)) < p_cd)
			{
				Solid[cell] = 1;
			}
		}
			
	cout << "Generate growth core completed" << endl;

	///  Generate process
	do
	{
		grow(M, N, Solid);

		/// Calculate porosity
		float t_temp = 0;
		#pragma omp paralle for
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				const int cell = CellIndex(i,j,M,N);
				t_temp += Solid[cell];
			}

		phi_p = 1 - t_temp / (N * N);
		cout << phi_p << endl;

	} while (phi_p > phi);
	cout << "Grow completed" << endl;

	/// offset at a large area
	const int Total_M = 500;
	const int Total_N = 300;
	const int offset_x = 100;
	const int offset_y = 100;
	int * S = new int[Total_M * Total_N]{};
	offset(Total_M, Total_N, offset_x, offset_y, S, M, N, Solid);

	/// Output result
	// for tecplot 360
	output2tecplot(Total_M, Total_N, S, "offset_porous_2D.plt");
	output2tecplot(M, N, Solid, "porous_2D.plt");

	// for gnuplot
	//output2gnuplot(Total_M, Total_N, S);

	delete[] S;
	delete[] Solid;
	return 0;
}