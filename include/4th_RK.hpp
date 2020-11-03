#include <vector>

#ifndef _4TH_RK_HPP_
#define _4TH_RK_HPP_

//	maximum number of iterations
#ifndef _ITERATE_MAX_NUM_
	#define _ITERATE_MAX_NUM_	20	// 2^10 = 1024
#endif

//	===========================================================================
// 	Define a structure to store the right hand side of the ODEs and parameters
//  that are needed to perform the calculations of the right hand side.
struct ODE_System{
//	dimension of the ODEs
    int Ndim;
//	right hand side the ODEs
    void (*Func)(double t, double *x, double *f, void *params);
//	parameters to be passed to the above "void func(****)"
    void *Params;

	ODE_System();
	~ODE_System();
	void Init(  void func(double t, double *x, double *f, void *params),
				void *params,
				int ndim);
};


//	===========================================================================
//	Define a structure to store vaiables needed by the 4th-RK algrithm.
struct ODE_RK4th_Workspace{
	int Ndim;
	double *k1;
	double *k2;
	double *k3;
	double *k4;

	ODE_RK4th_Workspace();
	~ODE_RK4th_Workspace();

	void Init(int ndim);
};

void ODE_Solver_RK4th_Evalute_Ks( 	ODE_RK4th_Workspace *rk4th,
									double t,
									double x[],
									double h,
									void func(double t, double *x, double *f, void *params),
									void *params
								 );

//	===========================================================================
//	Evolve the ODE system from x to x + h.
void Evolve_OneStep_RK4th(  ODE_RK4th_Workspace *rk4th,
							ODE_System *ode_sys,
							double t,					//	initial position
                            double x[],					//	initial values
                            double h,					//	step_size
                            double eps_rel,				//	relative error
                            double eps_abs );			//	absolute error

void ODE_Solver_RK4th( 	ODE_RK4th_Workspace *rk4th,
                        ODE_System *sys,
						double &t0,
                        double t1,
                        double x0[],
						double eps_rel,
						double eps_abs );

#endif