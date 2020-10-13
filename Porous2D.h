#ifndef POROUS2D_H
#define POROUS2D_H

#include <string>
class Porous2D
{
    public:
		const int NX;         /// grid number in x direction
        const int NY;         /// grid number in y direction
        double phi;	/// target porosity
        const double p_cd;		/// core distribution probability
		const double z;     /// Probability in orthogonal direction
        const double f;   /// Probability in oblique direction

        /// Constructor
        Porous2D(int nx, int ny, double por, double p, double p_z, double p_f);
		~Porous2D();

		void ReadTecplot2D(const int Total_M, const int Total_N, int* S, std::string filename);
		void Generation(int* s);
        void output2tecplot(const int Total_M, const int Total_N, int* S, std::string filename);
		void output2png(const int Total_M, const int Total_N, int* S, std::string filename);
		bool smooth();
		bool fixhole();
		void DeleteDeadZone();
		double CalcPor();
		int NearSolid(int i, int j);
		void offset(const int Total_M, const int Total_N, const int offset_x,
			const int offset_y, int* S, const int M, const int N, int* Solid);
        
    private:
        double p_d1 = z;
		double p_d2 = z;
		double p_d3 = z;
		double p_d4 = z;
        double p_d5 = f;
		double p_d6 = f;
		double p_d7 = f;
		double p_d8 = f;
        int * Solid;      /// the return value
        int Grow_Times;         /// for iteration

		void grow_alone();
		void grow_d1(int i, int j, int* s, int CoreID);
		void grow_d2(int i, int j, int* s, int CoreID);
		void grow_d3(int i, int j, int* s, int CoreID);
		void grow_d4(int i, int j, int* s, int CoreID);
		void grow_d5(int i, int j, int* s, int CoreID);
		void grow_d6(int i, int j, int* s, int CoreID);
		void grow_d7(int i, int j, int* s, int CoreID);
		void grow_d8(int i, int j, int* s, int CoreID);

		int CellIndex(int x, int y, const int MM, const int NN);
		void InfectFluid(int i, int j);
};

#endif
