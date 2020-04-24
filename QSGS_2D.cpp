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

using namespace std;

const float phi = 0.6;			/// target porosity
const float p_cd = 0.05;		/// core distribution probability
const float z = 0.05, f = 0.0125;
const float p_d1 = z, p_d2 = z, p_d3 = z, p_d4 = z;
const float p_d5 = f, p_d6 = f, p_d7 = f, p_d8 = f;

const int N = 300;
const int M = 300;
int Grow_Times = 0;

int CellIndex(int x, int y, const int M, const int N)
{
	assert(x >= 0 && x < M);
	assert(y >= 0 && y < N);

	return y + x * M;
}

/// No solid in the beginning.
int* Solid = new int[M * N]{};

/// \brief Grow in eight directions.
///
/// Grow directions
///*****    6    2    5   *****
///*****                  *****
///*****    3    C    1   *****
///*****                  *****
///*****    7    4    8   *****
///
void grow()
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

void output2tecplot(const int Total_M, const int Total_N, int *S)
{
	string out_buffer = "";

	ofstream outfile;
	string filename = "porous_2D.plt";
	string title = "2D porous media";
	outfile.open(filename);

	outfile << "TITLE = \"" << title << "\"" << endl;

	/// output solid position
	outfile << "VARIABLES = \"X\", \"Y\", \"Z\",\"value\" " << endl;
	outfile << "ZONE  t=\"solid\" I = " << Total_M << ", J = " << Total_N << ", K = " << 1 << ", F = point" << endl;
	for (int i = 0; i < Total_M; i++)
		for (int j = 0; j < Total_N; j++)
		{
			const int cell = CellIndex(i,j,Total_M, Total_N);
			out_buffer += std::to_string(i) + "," + std::to_string(j) + ", 0, " + std::to_string(S[cell]) + "\n";
		}
	outfile << out_buffer;
	outfile.close();
	cout << "Output completed" << endl;
}

void output2gnuplot(const int Total_M, const int Total_N, int *S)
{
	string out_buffer = "";
	ofstream outfile;
	string filename = "porous.out";

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

	system("pause");
}

void offset(const int Total_M, const int Total_N, const int offset_x, const int offset_y, int *S)
{
	if (Total_M <= (offset_x + N) || Total_N <= (offset_y + N))
	{
		std::cout << "No need to offset purous aera." << endl;
	}
	else
	{
		#pragma omp parallel for
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				const int cell = CellIndex(i,j,M,N);
				const int Total_cell = CellIndex(offset_x + i, offset_y + j, Total_M, Total_N);
				if (Solid[cell] == 1)	S[Total_cell] = 1;
			}	
	}
}


int main()
{
	float phi_p;

	/// Generate random seed
	srand((unsigned)time(NULL));
	cout << "Generate random seed completed" <<endl;

	/// Initialization
	init();
	cout << "Initialization completed" << endl;

	/// Generate growth core
	#pragma omp paralle for
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
		{
			const int cell = CellIndex(i,j);
			if ((rand() / double(RAND_MAX)) < p_cd)
			{
				Solid[cell] = 1;
			}
		}
			
	cout << "Generate growth core completed" << endl;

	///  Generate process
	do
	{
		grow();

		/// Calculate porosity
		float t_temp = 0;
		#pragma omp paralle for
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				const int cell =CellIndex(i,j);
				t_temp += Solid[i][j];
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
	int* S = new int[Total_M*Total_N]{};
	offset(Total_M, Total_N, offset_x, offset_y, S);

	/// Output result
	// for tecplot 360
	output2tecplot(Total_M, Total_N, S);

	// for gnuplot
	//output2gnuplot(Total_M, Total_N, S);

	system("pause");

	return 0;
}