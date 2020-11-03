#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>

#include "4th_RK.hpp"

using namespace std;


ODE_System::ODE_System(){
//	do nothing
}

ODE_System::~ODE_System(){
//	do nothing
}

//	Initialize the ODE system
void ODE_System::Init( void func(double t, double *x, double *f, void *params),
                       void *params,
                       int ndim ) {
    Ndim	= ndim;
    Func	= func;
    Params	= params;
}


//	Constructor
ODE_RK4th_Workspace::ODE_RK4th_Workspace(){
	k1 = NULL;
	k2 = NULL;
	k3 = NULL;
	k4 = NULL;
}

void ODE_RK4th_Workspace::Init(int ndim){
	Ndim = ndim;
	k1 = new double[ndim];
	k2 = new double[ndim];
	k3 = new double[ndim];
	k4 = new double[ndim];
}

//	Destructor
ODE_RK4th_Workspace::~ODE_RK4th_Workspace(){
	if( k1 != NULL)
		delete[] k1;
	if( k2 != NULL )
		delete[] k2;
	if( k3 != NULL )
		delete[] k3;
	if( k4 != NULL )
		delete[] k4;
}

void ODE_Solver_RK4th_Evalute_Ks( 	ODE_RK4th_Workspace *rk4th,
									double t,
									double x[],
									double h,
									void func(double t, double *x, double *f, void *params),
									void *params
								 ){
	int ndim = rk4th->Ndim;
	double f[ndim];
	double xtmp[ndim];

	for( int i=0; i<ndim; ++i ){
		xtmp[i] = x[i];
	}

//	Evaluate k1[]
	func(t, xtmp, f, params );
	for( int i=0; i<ndim; ++i ){
		rk4th->k1[i] = f[i];
		xtmp[i] = x[i] + 0.5*h*rk4th->k1[i];
	}

//	Evaluate k2[]
	func(t+0.5*h, xtmp, f, params );
	for( int i=0; i<ndim; ++i ){
		rk4th->k2[i] = f[i];
		xtmp[i] = x[i] + 0.5*h*rk4th->k2[i];
	}

//	Evaluate k3[]
	func(t+0.5*h, xtmp, f, params );
	for( int i=0; i<ndim; ++i ){
		rk4th->k3[i] = f[i];
		xtmp[i] = x[i] + h*rk4th->k3[i];
	}

//	Evaluate k4[]
	func(t+h, xtmp, f, params );
	for( int i=0; i<ndim; ++i ){
		rk4th->k4[i] = f[i];
	}
 }

//	===========================================================================
//	Evolve the ODE system from x to x + h.
void Evolve_OneStep_RK4th(  ODE_RK4th_Workspace *rk4th,
							ODE_System *ode_sys,
							double t,					//	initial position
                            double x[],					//	initial values
                            double h,					//	step_size
                            double eps_rel,				//	relative error
                            double eps_abs 				//	absolute error
                         ) {

	int ndim = ode_sys->Ndim;
	double x_new1[ndim];
	double x_new2[ndim];

	for(int i=0; i<ndim; ++i ){
		x_new1[i] = x_new2[i] = 0;
	}

//	step1: evolve the ODE system from x to x+h by one step
	ODE_Solver_RK4th_Evalute_Ks(rk4th, t, x, h,
								ode_sys->Func,
								ode_sys->Params);

	for( int i=0; i<ndim; ++i ){
		x_new1[i] = x[i] + (rk4th->k1[i] + 2*rk4th->k2[i] + 2*rk4th->k3[i] + rk4th->k4[i])*h/6;
	}

	int iter=1;
	while( iter <= _ITERATE_MAX_NUM_ ){
		
		int JMAX = pow(2,iter);
		double h_tmp = h / JMAX;
		
		for( int i=0; i<ndim; ++i ){
			x_new2[i] = x[i];
		}

		// get x_new2 from x_old via 'JMAX' steps
		for( int j=0; j<JMAX; ++j ){
			double t_tmp = t + j*h_tmp;
			ODE_Solver_RK4th_Evalute_Ks(rk4th, t_tmp, x_new2, h_tmp, ode_sys->Func, ode_sys->Params);
			for( int k=0; k<ndim; ++k){
				x_new2[k] += (rk4th->k1[k] + 2*rk4th->k2[k] + 2*rk4th->k3[k] + rk4th->k4[k])*h_tmp/6;
			}
		}

		// Evaluate errors
		double eps=0;
		for( int i=0; i<ndim; ++i){
			double eps_tmp = abs(x_new2[i]-x_new1[i]);
			if( eps < eps_tmp )
				eps = eps_tmp;
		}

		// Compare errors
		if( eps <= eps_abs ){
			for( int i=0; i<ndim; ++i ){
				x[i] = x_new2[i];
			}
			break;
		} else { // If eps is not small enough, copy x_new2 into x_new1
			for( int i=0; i<ndim; ++i ){
				x_new1[i] = x_new2[i];
			}
		}

		iter++;
	}


	if( iter >= _ITERATE_MAX_NUM_ ){
	// 	print warning message
		std::cout << "==> WARING: iterate number exceeds maximum number: "
				  << _ITERATE_MAX_NUM_ << std::endl;
		std::cout << "==> consider to choose another eps_abs." << std::endl;
	}
}

void ODE_Solver_RK4th( 	ODE_RK4th_Workspace *rk4th,
                        ODE_System *sys,
						double& t0,
                        double t1,
                        double x0[],
						double eps_rel,
						double eps_abs ){
	double h = t1-(t0);
	Evolve_OneStep_RK4th( rk4th, sys, t0, x0, h, eps_rel, eps_abs );
	t0 = t1;
}
