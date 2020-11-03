#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

#include "4th_RK.hpp"
#include "misc.hpp"

using namespace std;

void Func( double t, double x[], double f[], void *params ){
	double *p = (double *)params;
	double G = p[0];
	double M = p[1];
	double r2 = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
	double r3 = pow(r2,1.5);
	f[0] = x[3];
	f[1] = x[4];
	f[2] = x[5];
	f[3] = -G*M/r3*x[0];
	f[4] = -G*M/r3*x[1];;
	f[5] = -G*M/r3*x[2];;

	// for( int i=0; i<6; ++i )
	// 	cout << "f[" << i << "] = " << f[i] << "x[" << i << "] = " << x[i] << endl;
	// exit(0);
}

int main()
{
	double p[2] = {1,1};
	int dim=6;

	double t=0;
	double t1=10;
	double x[6] = {1,0.1,0.1,0.1,0.5,0.3};

	ODE_System sys;
	sys.Init( Func, p, dim);

	ODE_RK4th_Workspace rk4th;
	rk4th.Init(dim);


	int NMAX = 1000;
	for( int i=1; i<=NMAX; ++i ){
		double ti = i*t1 / NMAX;
		ODE_Solver_RK4th( &rk4th, &sys, t, ti, x, 1e-6, 1e-6);
	// exit(0);
	
	//	contral output format
		cout.precision(10);
		cout.width(12);
		cout.unsetf(ios::left);
		cout << t << "\t";
		for( int j=0; j<6; ++j )
			cout << x[j] << " ";
		
		cout << endl;
	}

	return(0);
}


