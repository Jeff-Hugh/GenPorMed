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

const float phi = 0.1;
const float p_cd = 0.05;
const float z = 0.05, f = 0.0125;
const float p_d1 = z, p_d2 = z, p_d3 = z, p_d4 = z;
const float p_d5 = f, p_d6 = f, p_d7 = f, p_d8 = f;

const int N = 300;
int Solid[N][N];
int Grow_Times = 0;

/// No solid in the beginning.
void init()
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			Solid[i][j] = 0;
		}
	}
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
void grow()
{
	Grow_Times += 1;
	cout << "Grow Times: " << Grow_Times << endl;

	int Solid_p[N][N];
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			Solid_p[i][j] = Solid[i][j];
		}
	}
	cout << "Copy compeleted" << endl;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (i == 0)
			{
				if (j == 0)
				{
					if (Solid[i][j] == 1)
					{
						if ((rand() / double(RAND_MAX)) < p_d1) Solid_p[i + 1][j] = 1;
						if ((rand() / double(RAND_MAX)) < p_d2) Solid_p[i][j + 1] = 1;
						if ((rand() / double(RAND_MAX)) < p_d5) Solid_p[i + 1][j + 1] = 1;
					}
				}
				else if (j == N - 1)
				{
					if (Solid[i][j] == 1)
					{
						if ((rand() / double(RAND_MAX)) < p_d1) Solid_p[i + 1][j] = 1;
						if ((rand() / double(RAND_MAX)) < p_d4) Solid_p[i][j - 1] = 1;
						if ((rand() / double(RAND_MAX)) < p_d8) Solid_p[i + 1][j - 1] = 1;
					}
				}
				else
				{
					if (Solid[i][j] == 1)
					{
						if ((rand() / double(RAND_MAX)) < p_d1) Solid_p[i + 1][j] = 1;
						if ((rand() / double(RAND_MAX)) < p_d2) Solid_p[i][j + 1] = 1;
						if ((rand() / double(RAND_MAX)) < p_d4) Solid_p[i][j - 1] = 1;
						if ((rand() / double(RAND_MAX)) < p_d5) Solid_p[i + 1][j + 1] = 1;
						if ((rand() / double(RAND_MAX)) < p_d8) Solid_p[i + 1][j - 1] = 1;
					}
				}
			}
			else if (i == N - 1)
			{
				if (j == 0)
				{
					if (Solid[i][j] == 1)
					{
						if ((rand() / double(RAND_MAX)) < p_d2) Solid_p[i][j + 1] = 1;
						if ((rand() / double(RAND_MAX)) < p_d3) Solid_p[i - 1][j] = 1;
						if ((rand() / double(RAND_MAX)) < p_d6) Solid_p[i - 1][j + 1] = 1;
					}
				}
				else if (j == N - 1)
				{
					if (Solid[i][j] == 1)
					{
						if ((rand() / double(RAND_MAX)) < p_d3) Solid_p[i - 1][j] = 1;
						if ((rand() / double(RAND_MAX)) < p_d4) Solid_p[i][j - 1] = 1;
						if ((rand() / double(RAND_MAX)) < p_d7) Solid_p[i - 1][j - 1] = 1;
					}
				}
				else
				{
					if (Solid[i][j] == 1)
					{
						if ((rand() / double(RAND_MAX)) < p_d2) Solid_p[i][j + 1] = 1;
						if ((rand() / double(RAND_MAX)) < p_d3) Solid_p[i - 1][j] = 1;
						if ((rand() / double(RAND_MAX)) < p_d4) Solid_p[i][j - 1] = 1;
						if ((rand() / double(RAND_MAX)) < p_d6) Solid_p[i - 1][j + 1] = 1;
						if ((rand() / double(RAND_MAX)) < p_d7) Solid_p[i - 1][j - 1] = 1;
					}
				}
			}
			else
			{
				if (j == 0)
				{
					if (Solid[i][j] == 1)
					{
						if ((rand() / double(RAND_MAX)) < p_d1) Solid_p[i + 1][j] = 1;
						if ((rand() / double(RAND_MAX)) < p_d2) Solid_p[i][j + 1] = 1;
						if ((rand() / double(RAND_MAX)) < p_d3) Solid_p[i - 1][j] = 1;
						if ((rand() / double(RAND_MAX)) < p_d5) Solid_p[i + 1][j + 1] = 1;
						if ((rand() / double(RAND_MAX)) < p_d6) Solid_p[i - 1][j + 1] = 1;
					}
				}
				else if (j == N - 1)
				{
					if (Solid[i][j] == 1)
					{
						if ((rand() / double(RAND_MAX)) < p_d1) Solid_p[i + 1][j] = 1;
						if ((rand() / double(RAND_MAX)) < p_d3) Solid_p[i - 1][j] = 1;
						if ((rand() / double(RAND_MAX)) < p_d4) Solid_p[i][j - 1] = 1;
						if ((rand() / double(RAND_MAX)) < p_d7) Solid_p[i - 1][j - 1] = 1;
						if ((rand() / double(RAND_MAX)) < p_d8) Solid_p[i + 1][j - 1] = 1;
					}
				}
				else
				{
					if (Solid[i][j] == 1)
					{
						if ((rand() / double(RAND_MAX)) < p_d1) Solid_p[i + 1][j] = 1;
						if ((rand() / double(RAND_MAX)) < p_d2) Solid_p[i][j + 1] = 1;
						if ((rand() / double(RAND_MAX)) < p_d3) Solid_p[i - 1][j] = 1;
						if ((rand() / double(RAND_MAX)) < p_d4) Solid_p[i][j - 1] = 1;
						if ((rand() / double(RAND_MAX)) < p_d5) Solid_p[i + 1][j + 1] = 1;
						if ((rand() / double(RAND_MAX)) < p_d6) Solid_p[i - 1][j + 1] = 1;
						if ((rand() / double(RAND_MAX)) < p_d7) Solid_p[i - 1][j - 1] = 1;
						if ((rand() / double(RAND_MAX)) < p_d8) Solid_p[i + 1][j - 1] = 1;
					}
				}
			}
		}
	}
	cout << "Neighbor produced" << endl;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			Solid[i][j] = Solid_p[i][j];
		}
	}
	cout << "One time grow compeleted" << endl;
}


int main()
{
	float phi_p;

	/// Produce random seed
	srand((unsigned)time(NULL));
	cout << "Produce random seed completed" <<endl;

	/// Initialization
	init();
	cout << "Initialization completed" << endl;

	/// Produce growth core
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if ((rand() / double(RAND_MAX)) < p_cd)
			{
				Solid[i][j] = 1;
			}
		}
	}
	cout << "Produce growth core completed" << endl;

	///  Produce process
	do
	{
		grow();

		/// Calculate porosity
		float t_temp = 0;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				t_temp += Solid[i][j];
			}
		}
		phi_p = t_temp / (N * N);
		cout << phi_p << endl;

	} while (phi_p < phi);
	cout << "Grow completed" << endl;

	/// Output result
	ofstream outfile;
	string filename;
	filename = "porous.out";
	outfile.open(filename);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			outfile << Solid[i][j] << " ";
		}
		outfile << endl;
	}
	outfile.close();
	cout << "Output result completed" << endl;

	system("pause");

	return 0;
}