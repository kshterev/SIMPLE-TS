/*
 ============================================================================
 Name        : SIMPLE_TS_MPI
 Author      : Kiril S. Shterev
 Version     : v.1.1
 
 Copyright   : All rights reserved. The source code is freely distributed for non-commercial use.
 Non-commercial use: Developers or distributors can compile, use all code or any part of the code, redistribute, sell results calculated using code or part of it and not-only, but except commercial use.
 Commercial use    : It is consider any use of the code or part of it in any way related to software, which is sold. In this case has to be contacted to Kiril Shterev to negotiate terms and conditions of use. 
 In any usage of the code, the derivatives has to include the following statement: "This software contains source code provided by Kiril Shterev."
 
 In any case, that is used algorithm SIMPLE-TS has to be cited the main paper presented the algorithm [1] and any other related paper presented further development of the algorithm. The list of papers related to the algorithm are on web site contains source code of the algorithm or on the web page of the Kiril Shterev:
 http://www.imbm.bas.bg/index.php/en_US/pressure-based-finite-volume-method
 http://www.imbm.bas.bg/index.php/en_US/kiril-stoyanov-shterev
 
 No Support    : Kiril Shterev has no obligation to support or to continue providing or updating any of Materials.
 No Warranties : This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 
 Description : Computational Fluid Dynamic using method unsteady SIMPLE-TS, MPI 1.2, C++
 
 
 References:
 1. K. Shterev and S. Stefanov, Pressure based finite volume method for calculation of compressible viscous gas flows, Journal of Computational Physics 229 (2010) pp. 461-480,  doi:10.1016/j.jcp.2009.09.042
 ============================================================================
*/



//If is defined DEBUG_PROGRAM the program will write a lot of arrays to files
//to make debug easier.
//#define DEBUG_PROGRAM


//If is defined COMPILE_MPI the MPI will be compiled,
//otherwise no MPI will be compiled;).
//Just make or not #define COMPILE_MPI as comment.
#define COMPILE_MPI

#ifdef COMPILE_MPI

//Include 4 MPI_Barrier in loop2 and loop3 to stabilize the iterative process,
//when are used more processes (for 150 give error if are not included)
//#define USE_MPI_Barrier_in_iterative_process_to_be_more_stable
#define USE_MPI_Barrier_in_iterative_process_to_be_more_stable_loop2
//#define USE_MPI_Barrier_in_iterative_process_to_be_more_stable_loop3

#include <mpi.h>
#include "DomainDecomposition2D.h"

#endif //End of COMPILE_MPI


//Use fixed number of iteration in loop 3 (internal loop for calculation pressure and temperature)
//For unsteady pressure driven fluid flow in a 2D channel the results show that 2 or 4 iterations
//are good choice for good accuracy for a computational time
#define FIXED_NUMBER_OF_ITERATIONS_LOOP_3


// The solving of p and Temper have to be soved using Jacoby iterative method,
// when is used domain decomposition method (parralel - MPI).
#define Calculate_p_and_Temper_using_Jakoby_method


#ifdef Calculate_p_and_Temper_using_Jakoby_method
//Substitute P with p_old and TEMPER with Temper_old, when is used Jakoby iterative method
#define TEMPER Temper_old
#define P p_old
#else
//Substitute P with p and TEMPER with Temper, when is used Gauss-Seidel iterative method
#define TEMPER Temper
#define P p
#endif

#ifndef COMPILE_MPI
/*
Method of calculation p and Temper
Calculate p and Temper using Jakoby method
If variable Calculate_p_and_Temper_using_Jakoby_method is defined p and Temper
are calculated using method of Jacoby (p[ij] = a1p * p_old[i_1j] + ....).
If for loop for p and Temper is made one iteration, when is used Jacoby method, the
iteration process in some cases is not convergent. To make process converge can be made
some iterations, for example N_I_p_c = 5.

If variable Calculate_p_and_Temper_using_Jakoby_method is NOT defined p and Temper
are calculated using method of Gouse-Zaidel (p[ij] = a1p * p[i_1j] + ....)
*/

//#define Calculate_p_and_Temper_using_Jakoby_method

#endif //End of COMPILE_MPI


//From practical expiriance, when is used Jakoby method of calculation p and Temper
//have to make minimum 2 iteration, to garantee convergentbility and stablelity of this loop.
#ifdef Calculate_p_and_Temper_using_Jakoby_method
	const int MinimalNumberOfIterationInLoopFor_p_and_Temper = 2;
#else
	const int MinimalNumberOfIterationInLoopFor_p_and_Temper = 1;
#endif



#include "Objects_for_Solving.h"
#include "WriteReadData.h"
#include "stringConv.h"

#include "iostream"
#include <time.h>


#define _USE_MATH_DEFINES
#include <math.h>

#include <stdio.h>
#include <cstring>


//define all function which MUST be inline

//return maximum value
#define maximum(a, b) (((a) > (b)) ? (a) : (b))

//return maximum value
#define minimum(a, b) (((a) < (b)) ? (a) : (b))

//upwind shame - return rho * v
#define upwind(r0, r1, v1) (((v1) > 0) ? ((r0) * (v1)) : ((r1) * (v1)))

//upwind shame - return only rho
#define upwind_rho(r0, r1, v1) (((v1) == 0) ? (0.5 * ((r0) + (r1))) : (((v1) > 0) ? (r0) : (r1)))

//linear interpolation between two given values and steps in mesh
#define Li(fi_1, fi_2, h_1, h_2) (((fi_1) * (h_2) + (fi_2) * (h_1)) / ((h_1) + (h_2)))

//average harmonic between two given values and steps in mesh
#define Hi(fi_1, fi_2, h_1, h_2) (((h_1) + (h_2)) * (fi_1) * (fi_2) / ((h_1) * (fi_2) + (h_2) * (fi_1)))

//Solve index for points i+1, i-1, j+1, j-1
//This MACROS is NO need to be changed to be used in MPI.
#define ij_function(node1) \
	ij = (node1);\
	i_1 = i - 1 * (0 < i);\
	i_2 = i - 1 * (0 < i) - 1 * (1 < i);\
	i1 = i + 1 * (i < (Nx - 1));\
	i_1j = ij - 1 * (0 < i);\
	i1j = ij + 1 * (i < (Nx - 1));\
	ij_1 = ij - Nx * (0 < j);\
	ij1 = ij + Nx * (j < (Ny - 1));\
	i_1j_1 = ij - 1 * (0 < i) - Nx * (0 < j);\
	i_1j1 = ij - 1 * (0 < i) + Nx * (j < (Ny - 1));\
	i1j_1 = ij + 1 * (i < (Nx - 1)) - Nx * (0 < j);\
	i1j1 = ij + 1 * (i < (Nx - 1)) + Nx * (j < (Ny - 1));\
	i_2j = i_2 + j * Nx;



//Solve index for points i+1, i-1, j+1, j-1,
//for loop which start from i = 0, j = 0 and are to the i = Nx - 1, j = Ny -1.
//This MACROS is NEEDED to be changed to be used in MPI.
#ifndef COMPILE_MPI
#define ij_function_i_and_j_to_the_boundaries(node1) \
	ij = (node1);\
	i_1 = i - 1 * (0 < i);\
	i_2 = i - 1 * (0 < i) - 1 * (1 < i);\
	i1 = i + 1 * (i < (Nx - 1));\
	i_1j = ij - 1 * (0 < i);\
	i1j = ij + 1 * (i < (Nx - 1));\
	ij_1 = ij - Nx * (0 < j);\
	ij1 = ij + Nx * (j < (Ny - 1));\
	i_1j_1 = ij - 1 * (0 < i) - Nx * (0 < j);\
	i_1j1 = ij - 1 * (0 < i) + Nx * (j < (Ny - 1));\
	i1j_1 = ij + 1 * (i < (Nx - 1)) - Nx * (0 < j);\
	i1j1 = ij + 1 * (i < (Nx - 1)) + Nx * (j < (Ny - 1));\
	i_2j = i_2 + j * Nx;

#endif //End of COMPILE_MPI


#ifdef COMPILE_MPI
#define ij_function_i_and_j_to_the_boundaries(node1) \
	ij = (node1);\
	i_1 = i - 1 * (0 < i);\
	i_2 = i - 1 * (0 < i) - 1 * (1 < i);\
	i1 = i + 1 * (i < (Nx - 1));\
	i_1j = ij - 1 * (0 < i);\
	i1j = ij + 1 * (i < (Nx - 1));\
	ij_1 = ij - Nx * (0 < j);\
	ij1 = ij + Nx * (j < (Ny - 1));\
	i_1j_1 = ij - 1 * (0 < i) - Nx * (0 < j);\
	i_1j1 = ij - 1 * (0 < i) + Nx * (j < (Ny - 1));\
	i1j_1 = ij + 1 * (i < (Nx - 1)) - Nx * (0 < j);\
	i1j1 = ij + 1 * (i < (Nx - 1)) + Nx * (j < (Ny - 1));\
	i_2j = i_2 + j * Nx;

#endif //End of COMPILE_MPI



//Solve index for points i+1, i-1, i-2, j+1, j-1  for Periodic Boundary Conditions
//This MACROS have to be checked to be used in MPI.
#define ij_function_PeriodicBC(ipBC, jpBC) \
	if((ipBC) == 0) {i_1 = Nx - 1;} else {i_1 = (ipBC) - 1;}\
	if((ipBC) == 1) {i_2 = Nx - 1;} else if((ipBC) == 0) {i_2 = Nx - 2;} else {i_2 = (ipBC) - 2;}\
	if((ipBC) == (Nx - 1)) {i1 = 0;} else {i1 = (ipBC) + 1;}\
	ij = (ipBC) + (jpBC) * Nx;\
	i_1j = i_1 + (jpBC) * Nx;\
	i1j = i1 + (jpBC) * Nx;\
	ij_1 = (ipBC) + (jpBC - 1 * (0 < j)) * Nx;\
	ij1 = (ipBC) + (jpBC + 1 * (j < (Ny - 1))) * Nx;\
	i_1j_1 = i_1 + (jpBC - 1 * (0 < j)) * Nx;\
	i_1j1 = i_1 + (jpBC + 1 * (j < (Ny - 1))) * Nx;\
	i1j_1 = i1 + (jpBC - 1 * (0 < j)) * Nx;\
	i1j1 = i1 + (jpBC + 1 * (j < (Ny - 1))) * Nx;\
	i_2j = i_2 + (jpBC) * Nx;

//Solve index for points i+1, i-1, i-2, j+1, j-1  for Periodic Boundary Conditions
//for loop which start from i = 0, j = 0 and are to the i = Nx - 1, j = Ny -1.
//This MACROS have to be checked to be used in MPI.
#define ij_function_PeriodicBC_i_and_j_to_the_boundaries(ipBC, jpBC) \
	if((ipBC) == 0) {i_1 = Nx - 1;} else {i_1 = (ipBC) - 1;}\
	if((ipBC) == 1) {i_2 = Nx - 1;} else if((ipBC) == 0) {i_2 = Nx - 2;} else {i_2 = (ipBC) - 2;}\
	if((ipBC) == (Nx - 1)) {i1 = 0;} else {i1 = (ipBC) + 1;}\
	ij = (ipBC) + (jpBC) * Nx;\
	i_1j = i_1 + (jpBC) * Nx;\
	i1j = i1 + (jpBC) * Nx;\
	ij_1 = (ipBC) + (jpBC - 1 * (0 < j)) * Nx;\
	ij1 = (ipBC) + (jpBC + 1 * (j < (Ny - 1))) * Nx;\
	i_1j_1 = i_1 + (jpBC - 1 * (0 < j)) * Nx;\
	i_1j1 = i_1 + (jpBC + 1 * (j < (Ny - 1))) * Nx;\
	i1j_1 = i1 + (jpBC - 1 * (0 < j)) * Nx;\
	i1j1 = i1 + (jpBC + 1 * (j < (Ny - 1))) * Nx;\
	i_2j = i_2 + (jpBC) * Nx;


//Macros to change the values of two pointers to double arrays pointer_double_tmp.
//It is used when have to change the pointers p and p_old instead copy of arrays.
//Have to be defined pointer to double
#define Change_two_pointers_double(pointer_0, pointer_1) \
	pointer_double_tmp = pointer_1; \
	pointer_1 = pointer_0; \
	pointer_0 = pointer_double_tmp;

double * pointer_double_tmp;

using namespace std;

//const double pi = 3.14159265358979; //pi number to 15 singht

Objects_for_Solving * OforS; //Objects for Solving

//#ifndef SIMPLEST_ANL_compr_T_h_var_NoDim_Algorithm_H
//#define SIMPLEST_ANL_compr_T_h_var_NoDim_Algorithm_H

//#pragma once


	//Nimber of Dimention of rule vectors for begind and end of area - begin end = 2
	static const unsigned short ND = 2;


	//template <typename T_a, typename T_b>
	//inline double maximum(const T_a& a, const T_b& b)
	//{
	//	return (a > b ? a : b);
	//};

	//it, it_begin, Nt, ht - this is absolutely time in dimensionless values.
	//ht - step by time
	double it, it_begin, Nt, ht;
	unsigned int it_counter, Nx, Ny, Na, N_I, Iter, Iter_pr1, Iter_pr2;


	/*Criterial to continue iteration process.
	continue_iter == false -> stop iteration in this time step,
								because the convergence criterions are not satisfyed
	continue_iter == true -> continue iteration in this time step,
								because the convergence criterions are satisfyed

	continue_iter_p_c == false -> stop iteration in loop, where are calculated p and Temper,
								because the convergence criterions are not satisfyed
	continue_iter_p_c == true -> continue iteration in loop, where are calculated p and Temper,
								because the convergence criterions are satisfyed
	*/
	bool continue_iter, continue_iter_p_c;

	//countinue_time_step == true -> continue calculating next time step
	//countinue_time_step == false -> stop calculating time step, stop program
	bool countinue_time_step;


	void Write_it_ToFile(void);
	void Read_it_FromFile(void);

	//Change u_gas_nd in time of solve.
	//use_interpolation_for_u_BC_by_time == true; will use interpolation for u_BC
	//use_interpolation_for_u_BC_by_time == false; will NOT use interpolation for u_BC
	bool use_interpolation_for_u_BC_by_time;
	//u_gas_nd will be solved using linear interpolation.
	//u_BC_xb_tb - BC of u_gas_nd in time t_BC_b
	//u_BC_xb_te - BC of u_gas_nd in time t_BC_e
	double u_BC_xb_t_BC_b, u_BC_xb_t_BC_e;
	//Time for what is given u_gas_nd. Time is equivalent to nondimensionall time.
	//if(it < t_BC_b) then u_gas_nd = u_BC_xb_t_BC_b;
	//if(t_BC_e < it) then u_gas_nd = u_BC_xb_t_BC_e;
    double t_BC_b, t_BC_e;

	//Apply Boundary Conditions which are function of time.
	inline void BoundaryConditions_function_of_time(void);


	unsigned int i, j, ij, i_1j, i1j, ij_1, ij1, i_1j_1, i_1j1, i1j_1, i1j1, i_2j,
		i_1, i1, i_2,
		node, counter;


	//double rho_gas, mu_gas;

	//variables for mesh
	bool is_defined_array_for_mesh;
	double * x_f, * y_f, * hx, * hy, * x_v, * y_v;
	double h_tmp;

	//kind of mesh
	//                                                                          _____  __  ______
	//kind_of_mesh = Li_mesh - mesh is with linear interpolation step of kind:       \/  \/
	static const unsigned int Li_mesh = 1;
	//                                                                          _____  __  ______
	//kind_of_mesh = Par_mesh - mesh is with parabolik interpolation step:           \/  \/
	static const unsigned int Par_mesh = 2;
	//                                                                          ______    _______
	//kind_of_mesh = Li_flat_mesh - mesh is with linear interpolation step:           \__/
	static const unsigned int Li_mesh_flat = 3;
	//kind_of_mesh = Polynom3_flat_mesh - mesh is with polynom from 3 degree.
	//That mean that derivatiors in two given points are zero,                  ______    _______
	//smooth on the two ends of interpoletion:                                        \__/
	static const unsigned int Par_mesh_flat = 4;
	//                                                                          ______    _______
	//kind_of_mesh = Par_flat_mesh - mesh is with parabolik interpolation step:       \__/
	static const unsigned int Polynom3_flat_mesh = 5;
	unsigned int kind_of_mesh;

	//variable for mesh:
	double hxmin, hxmax, hxmid, xfmin, xfmax, xbmin, xbmax, xmid,
		axf3, axf2, axf1, axf0, axb3, axb2, axb1, axb0;

	double hxfront(double& x);
	double hxback(double& x);
	double hxfront_Polynom3(double& x);
	double hxback_Polynom3(double& x);
	double hx_all(double& x);
	double hx_all_par(double& x);
	double hx_all_par_flat(double& x);
	double hx_all_linear(double& x);
	double hx_all_linear_flat(double& x);
	double hx_all_Polynom3_flat_mesh(double& x);


	double hymin, hymax, hymid, ytmin, ytmax, ybmin, ybmax, ymid,
		ayt3, ayt2, ayt1, ayt0, ayb3, ayb2, ayb1, ayb0;

	double hytop(double& y);
	double hybottom(double& y);
	double hytop_Polynom3(double& y);
	double hybottom_Polynom3(double& y);
	double hy_all(double& y);
	double hy_all_par(double& y);
	double hy_all_par_flat(double& y);
	double hy_all_linear(double& y);
	double hy_all_linear_flat(double& y);
	double hy_all_Polynom3_flat_mesh(double& y);

	double x_b, x_e, y_b, y_e; //this are area boundaries


	//define rule vectors
	//The target of rule vectors is to make loop without any check.

	//define rule vectors.
	void define_rule_vectors(void);

	//is_rule_vectors_defined == true - rule vectors are already defined
	//is_rule_vectors_defined == false - rule vectors are NOT defined
	bool is_rule_vectors_defined;

	//vol_inf_fb_bool == true -> if in volume is fluid
	//vol_inf_fb_bool == false -> if in volume is body
	//The target is to eliminate the check in loop.
	//Example 1:
	//if(vol_inf_fb[ij] > fluid) fi = 0;
	//else fi = something;
	//equivalence: fi = vol_inf_fb_bool[ij] * something;

	//Example 2:
	//if(vol_inf_fb[ij] > fluid || vol_inf_fb[i1j] > fluid) fi = 0;
	//else fi = something;
	//equivalence: fi = vol_inf_fb_bool[ij] * vol_inf_fb_bool[i1j] * something;
	bool * vol_inf_fb_bool;

	//define rule vectors and variables for u
	//That is posible because the result from checks are the same EVERY LOOP!
	//Define vectors and vatriables for area where is no check
	//number of all elements of sloving for u of No Check - Na_nch_u
	unsigned int Na_nch_u;
	//vector contain ij of all elements of sloving for u of No Check - rule_vector_ij_nch_u
	unsigned int * rule_vector_ij_nch_u;
	//vector contain j of all elements of sloving for u of No Check - rule_vector_j_nch_u
	//becouse ij = i + j * Nx => i = ij - j * Nx;
	unsigned int * rule_vector_j_nch_u;


	//Define vectors and vatriables for area where is Boundary Conditions - here we have some checks
	//number of all elements of sloving for u of Boundary Conditions - Na_bc_u
	unsigned int Na_bc_u;
	//vector contain ij of all elements of sloving for u of Boundary Conditions - rule_vector_ij_bc_u
	unsigned int * rule_vector_ij_bc_u;
	//vector contain j of all elements of sloving for u of Boundary Conditions - rule_vector_j_bc_u
	//becouse ij = i + j * Nx => i = ij - j * Nx;
	unsigned int * rule_vector_j_bc_u;



	//define rule vectors and variables for v
	//That is posible because the result from checks are the same EVERY LOOP!
	//Define vectors and vatriables for area where is no check
	//number of all elements of sloving for v of No Check - Na_nch_v
	unsigned int Na_nch_v;
	//vector contain ij of all elements of sloving for v of No Check - rule_vector_ij_nch_v
	unsigned int * rule_vector_ij_nch_v;
	//vector contain j of all elements of sloving for v of No Check - rule_vector_j_nch_v
	//becouse ij = i + j * Nx => i = ij - j * Nx;
	unsigned int * rule_vector_j_nch_v;


	//Define vectors and vatriables for area where is Boundary Conditions - here we have some checks
	//number of all elements of sloving for v of Boundary Conditions - Na_bc_v
	unsigned int Na_bc_v;
	//vector contain ij of all elements of sloving for v of Boundary Conditions - rule_vector_ij_bc_v
	unsigned int * rule_vector_ij_bc_v;
	//vector contain j of all elements of sloving for v of Boundary Conditions - rule_vector_j_bc_v
	//becouse ij = i + j * Nx => i = ij - j * Nx;
	unsigned int * rule_vector_j_bc_v;



	//define rule vectors and variables for others - p, Temper, rho in center of volume
	//That is posible because the result from checks are the same EVERY LOOP!
	//Define vectors and vatriables for area where is no check
	//number of all elements of sloving for others of No Check - Na_nch_others
	unsigned int Na_nch_others;
	//vector contain ij of all elements of sloving for others of No Check - rule_vector_ij_nch_others
	unsigned int * rule_vector_ij_nch_others;
	//vector contain j of all elements of sloving for others of No Check - rule_vector_j_nch_others
	//becouse ij = i + j * Nx => i = ij - j * Nx;
	unsigned int * rule_vector_j_nch_others;


	//Define vectors and vatriables for area where is Boundary Conditions - here we have some checks
	//number of all elements of sloving for others of Boundary Conditions - Na_bc_others
	unsigned int Na_bc_others;
	//vector contain ij of all elements of sloving for others of Boundary Conditions - rule_vector_ij_bc_others
	unsigned int * rule_vector_ij_bc_others;
	//vector contain j of all elements of sloving for others of Boundary Conditions - rule_vector_j_bc_others
	//becouse ij = i + j * Nx => i = ij - j * Nx;
	unsigned int * rule_vector_j_bc_others;



	//sqrt_T_ff_V - this is diffusion coefficient (dynamic viscosity) for velocity equations (u and v).
	//sqrt_T_ff_V is interpolation of sqrt_T in points x_f[i], y_f[j].
	//double * sqrt_T_ff_V;

	//This is diffusion coefficient in control surface yf, from integration of diffusion terms of equation for u.
	double * Gamma_yf;

	//This is diffusion coefficient in control surface xf, from integration of diffusion terms of equation for v.
	double * Gamma_xf;

	//The fluxes in all area. They are the same for u and v.
	double * Fx, * Fy;

	//double au33, au323;
	//double au44, au424;

	//Data are stored in storage arrays.
	//Pointer u and u_pr change the array.
	//At this way is prevented copy of u_pr = u at the beginning of loop 1.
	double * storage_u_0, * storage_u_1,
		* u, //* u_old,
		* u_pr,
		* u_pseudo, u_gas,
		* D_ux, * D_uy,
		//D_u1, D_u2, D_u3, D_u4,
		//Fxu0, Fxu1, Fxu2, Fxu3, Fxu4,
		//Fyu3, Fyu4, Fyu23, Fyu24,
		F_xu1, F_xu2,
		au0, au1, au2, au3, au4, bu3,
		Scu, Spu;

	//This is temporary coefficient a:
	double a_tmp, a_tmp_ij, a_tmp_i1j, a_tmp_ij1;


	//Variables to apply velocity slip and temperature jump boundary conditions
	//Number of body in control volume, when apply boundary conditions
	int Number_of_body;


	//Coefficients to apply temperature jump boundary condition
	bool is_TemperatureJumpBC_body_i_1_tmp, is_TemperatureJumpBC_body_i1_tmp,
		is_TemperatureJumpBC_body_j_1_tmp, is_TemperatureJumpBC_body_j1_tmp;

	bool dTdn_on_wall_0_i_1_tmp, dTdn_on_wall_0_i1_tmp,
		dTdn_on_wall_0_j_1_tmp, dTdn_on_wall_0_j1_tmp;

	double T_body_i_1_tmp, T_body_i1_tmp, T_body_j_1_tmp, T_body_j1_tmp;
	double aTBC_i_1, aTBC_i1, aTBC_j_1, aTBC_j1;



	////contain information for u. If u > 0 is true, if u <= 0 is false.
	//bool * positive_u;

	double max_residual_in_u, max_residual_in_u_Iter_0;
	//double au_c0, au_c1, au_c2, au_c3, au_c4, bu_c;
	//double au_c33, au_c323, au_c44, au_c424;
	//double aut0, aut1, aut2, aut3, aut4, but;
	//double autce0, autce1, autce2, autce3, autce4, btuce;
	//double aFu_c1, aFu_c2;
	//double Fxu333, Fxu23323, Fxu444, Fxu24424;


	//double av11, av141;
	//double av22, av242;

	//Data are stored in storage arrays.
	//Pointer v and v_pr change the array.
	//At this way is prevented copy of v_pr = v at the beginning of loop 1.
	double * storage_v_0, * storage_v_1,
		* v, //* v_old,
		* v_pr,
		* v_pseudo, v_gas,
		* D_vx, * D_vy,
		//D_v1, D_v2, D_v3, D_v4,
		//* Fxv, * Fyv,
		//Fxv1, Fxv2, Fxv41, Fxv42,
		//Fyv0, Fyv1, Fyv2, Fyv3, Fyv4,
		F_yv3, F_yv4,
		av0, av1, av2, av3, av4, bv3,
		Scv, Spv;



	//contain information for v. If v > 0 is true, if v <= 0 is false.
	//bool * positive_v;

	double max_residual_in_v, max_residual_in_v_Iter_0;
	//double av_c0, av_c1, av_c2, av_c3, av_c4, bv_c;
	//double av_c11, av_c141, av_c22, av_c242;
	//double avt0, avt1, avt2, avt3, avt4, bvt;
	//double avtce0, avtce1, avtce2, avtce3, avtce4, btvce;
	//double aFv_c3, aFv_c4;
	//double Fyv111, Fyv41141, Fyv222, Fyv42242;


	//Data are stored in storage arrays.
	//Pointer p and p_pr change the array.
	//At this way is prevented copy of p_pr = p at the beginning of loop 1 and p_old = p in loop 3.
	double * storage_p_0, * storage_p_1,
		//* p_c, * p_c_old,
		* p, * p_pr,
		* d_p_c_12, * d_p_c_34,
		ap0,
		//ap1, ap2, ap3, ap4,
		* apx, * apy,
		bp0,
		MaxError_p_c;
	//double fabs_max_p;

#ifdef Calculate_p_and_Temper_using_Jakoby_method
	double  * storage_p_2, * p_old;
#endif //End of Calculate_p_and_Temper_using_Jakoby_method

	//double bp01, bp02, bp03, bp04;
	double * bpx, * bpy;
	double hyj_upwind_rho_1_0, hxi_upwind_rho_3_0;
	//double rho_0_p, rho_1_p, rho_2_p, rho_3_p, rho_4_p;

	//double max_residual_in_p_c;
	//double max_residual_in_p_c_Itre0;
	//double sum_ap_c_and_bpc;

	unsigned int Iter_p_c, N_I_p_c;

	//density - rho
	//Data are stored in storage arrays.
	//Pointer rho and rho_pr change the array.
	//At this way is prevented copy of rho_pr = rho at the beginning of loop 1.
	double * storage_rho_0, * storage_rho_1, * rho, * rho_pr;
	////rho_1[ij] = 1 / rho[ij];
	////Many coefficient must be zero in body and in this way we do not need to check.
	//double * rho_1;

	//rho_u[ij] = upwind_rho(rho[i_1j], rho[ij] ,u[ij])
	//rho_u[ij] - contain rho in point u[ij] (x_f[i], y_v[j])
	double * rho_u;

	//rho_v[ij] = upwind_rho(rho[ij_1], rho[ij] ,v[ij])
	//rho_v[ij] - contain rho in point v[ij] (x_v[i], y_f[j])
	double * rho_v;

	//define_OforS
	void define_OforS(void);
	//Number of matrixes of vol_inf
	unsigned int N_vol_inf;

	//vol_inf - matrixes witch give us information aboit point
	void define_vol_inf(void);

	//This variable contain information about fluid and body in point
	static const unsigned int fluid = 0; //if vol_inf_fb[ij] is 0 that mean that there is fluid for solving
	//Number of bodies circles and polyhedrons which can be counted in vil_inf is greater then fluid
	//Ncp_beg_b - Number of circle or polyhedron begin body
	static const unsigned int Ncp_beg_b = fluid + 1;

	//If vol_inf_fb[ij] == fluid - in ij poin is fluid
	//If vol_inf_fb[ij] == 3 - in ij poin is body number 3 - 1 = 2
	unsigned int * vol_inf_fb;


	//This variable contain information, where in domain is GIVEN p
	static const unsigned int no_gp = 0; //if vol_inf_p[ij] is 0 that mean that there is no given pressure
	//Number of bodies circles and polyhedrons which can be counted in vil_inf is greater then fluid
	//Ncp_beg_p - Number of circle or polyhedron bedin pressure
	static const unsigned int Ncp_beg_p = no_gp + 1;

	//If vol_inf_p[ij] == flase - in ij poin is NOT GIVEN p
	//If vol_inf_p[ij] == true - in ij poin is GIVEN p
	unsigned int * vol_inf_p;
	void ApplyPressureConditionsToPressure(void);


	//This variable contain information, where in domain is GIVEN V(u,v)
	static const unsigned int no_gV = 0; //if vol_inf_V[ij] is 0 that mean that there is no given Velocity
	//Number of bodies circles and polyhedrons which can be counted in vil_inf is greater then fluid
	//Ncp_beg_V - Number of circle or polyhedron bedin pressure
	static const unsigned int Ncp_beg_V = no_gV + 1;

	//If vol_inf_V[ij] == flase - in ij poin is NOT GIVEN V
	//If vol_inf_V[ij] == true - in ij poin is GIVEN V
	unsigned int * vol_inf_V;
	void ApplyVelocityConditionsToVelocities(void);


	//data for body
	bool ToImport_data_for_circles_from_file_b,
		ToImport_data_for_polyhedrons_from_file_b;

	//data for given pressure
	bool ToImport_data_for_circles_from_file_p,
		ToImport_data_for_polyhedrons_from_file_p;

	//data for given velocity
	bool ToImport_data_for_circles_from_file_V,
		ToImport_data_for_polyhedrons_from_file_V;

	//variable which contain where in matrix OforS is information for
	//body, pressure and Velocity
	static const int gb = 0;
	static const int gp = 1;
	static const int gV = 2;


	void Solve(void);


	void define_Nx_Ny(void);

	void WriteEnteredDataToFile(void);
	void ReadEnteredDataFromFile(void);

	bool ToReadSolvedDataFromFile;

	//ToReadWriteSolvedDataToBinaryFile == true --> read solved data from binary file
	//ToReadWriteSolvedDataToBinaryFile == false --> read solved data from text file
	bool ToReadSolvedDataFromBinaryFile;

	//ToWriteSolvedDataToBinaryFile == true --> write solved data to binary file
	//ToWriteSolvedDataToBinaryFile == false --> write solved data to text file
	bool ToWriteSolvedDataToBinaryFile;

	//ToWriteSolvedDataToNewFiles == true --> write solved data to new file.
	//Example for u: u.bin.100 - this is velocity u, solved in binary mode, it = 100
	//To write date to binary ot text mode depend from ToWriteSolvedDataToBinaryFile varieble.
	//ToWriteSolvedDataToNewFiles == false --> do not write solved data in new files.
	bool ToWriteSolvedDataToNewFiles;

	//This is true when we continue calculation readind data which are interpoleted
	//from croase mesh.
	bool ToContinueFromInterpolation;

	void WriteSolvedDataToFile(void);
	void ReadSolvedDataFromFile(void);

	void WriteTempDataToFile(void);
	//ToCheckForFileIsNow_Tmp == true - for first if file exist and the solve is begin will delete exist file
	//ToCheckForFileIsNow_Tmp == false - if solve is begin file will not be deleted
	bool ToCheckForFileIsNow_Tmp;

	void WriteIterationsInLoopsToFile(void);
	//ToCheckForFileIsNow_Iterations_in_loops == true - for first if file exist and the solve is begin will delete exist file
	//ToCheckForFileIsNow_Iterations_in_loops == false - if solve is begin file will not be deleted
	bool ToCheckForFileIsNow_Iterations_in_loops;


	//Write data in short format for fast reading.
	//The kind of format of Short Data:
	//if WriteDataInShortFormat_kind == 0; ==> to not write data in short format
	//if WriteDataInShortFormat_kind == 1; ==> write: Re Kn CD CD_p u_gas_nd q_xb q_xe p_BC_xb p_BC_xe
	int WriteDataInShortFormat_kind;


	//Convergents criterions
	double MaxError_Velocities;
	////Epsilon1 = 1e-4;
	//double Epsilon1;
	////Epsilon2 = 1e-6;
	//double Epsilon2;
	////Epsilon3 = 5 * 1e-5;
 //   double Epsilon3;

	//// fabs_rho_u_h0_max_Iter0 = fabs(gas.rho_gas * u[ij] / h[0]) maximum after Iter = 0
	//double fabs_rho_u_h0_max_Iter0;
	////Full velocity after Iter = 0. V = sqrt(u^2 + v^2)
	//double Vmax_Iter0;
	//// fabs_u_u_old_max_Iter0 = fabs(u[ij] - u_old[ij]) maximum after Iter = 0
	//double fabs_u_u_old_max_Iter0;
	//// fabs_v_v_old_max_Iter0 = fabs(v[ij] - v_old[ij]) maximum after Iter = 0
	//double fabs_v_v_old_max_Iter0;
	double u_new, v_new;
	double u_tmp, v_tmp, p_tmp, Temper_tmp;


	//Temporary variables used for applying of Velocity Slip BC
	double u_body_tmp, v_body_tmp, w_VelocitySlipBC_tmp, F_VelocitySlip_tmp;

	bool is_VelocitySlipBC_body_tmp;


	//Temporary variables used for applying of Temperature Jump BC
	double T_body_tmp, dTdn_on_wall_0_tmp, F_TemperatureJump_tmp,
		F_TemperatureJump_Kn_for_this_BC_tmp;

	bool is_TemperatureJumpBC_body_tmp;


	// fabs_rho_u_h0_max_Iter0 = fabs(gas.rho_gas * v[ij] / h[1]) maximum after Iter = 0
	double fabs_rho_v_h1_max_Iter0;
	// fabs_rho_u_h0_max = fabs(gas.rho_gas * u[ij] / h[0]) maximum after check
	double fabs_rho_u_h0_max;
	// fabs_rho_v_h1_max = fabs(gas.rho_gas * v[ij] / h[1]) maximum after check
	double fabs_rho_v_h1_max;

	//Full velocity after check. V = sqrt(u^2 + v^2)
	double Vmax;
	// fabs_u_u_old_max = fabs(u[ij] - u_old[ij]) maximum after check
	double fabs_u_u_old_max;
	// fabs_v_v_old_max = fabs(v[ij] - v_old[ij]) maximum after check
	double fabs_v_v_old_max;

	int Nt_save_solved_data, Nt_save_solved_data_counter;
    int Nt_save_temp_data, Nt_save_temp_data_counter;

	//Nimber of iterations, which are passed before check process for convergence
	int N_I_check, N_I_check_counter;






	void define_area_for_solving(void);

	//Varieble which rules area for solving for velocities, pressure and etc.
	//Nx_u[0] begin area for solving for u
	//Nx_u[1] end area for solving for u
	unsigned int * Nx_u;
	//Ny_u[0] begin area for solving for u
	//Ny_u[1] end area for solving for u
	unsigned int * Ny_u;

	//Nx_v[0] begin area for solving for v
	//Nx_v[1] end area for solving for v
	unsigned int * Nx_v;
	//Ny_v[0] begin area for solving for v
	//Ny_v[1] end area for solving for v
	unsigned int * Ny_v;

	//Nx_others[0] begin area for solving for others
	//Nx_others[1] end area for solving for others
	unsigned int * Nx_others;
	//Ny_others[0] begin area for solving for others
	//Ny_others[1] end area for solving for others
	unsigned int * Ny_others;


	//variables for energy equation
	//Data are stored in storage arrays.
	//Pointer Temper and Temper_pr change the array.
	//At this way is prevented copy of Temper_pr = Temper at the beginning of loop 1 and Temper_old = Temper in loop 3.
	double * storage_Temper_0, * storage_Temper_1,
		* sqrt_T, cT, * Temper,
		* Temper_pr,
		cT1, cT2,
		DT1, DT2, DT3, DT4, FT1, FT2, FT3, FT4,
		* DTx, * DTy,
		aT0 , aT1, aT2, aT3, aT4, bT0,
		* ScT1,
		FiT,
		dudx_T, dvdy_T, dudy_T, dvdx_T,
		//sqrt on boundary of area
		sqrt_T_T_1b, sqrt_T_T_2b, sqrt_T_T_3b, sqrt_T_T_4b,
		MaxError_T;


	//This is Temperature of the fluid on the surface of the body.
	//If There is no temperature jump Temperature of the fluid on the surface of the body
	//is equal to thetemperature of the body.
	//If there is temperature jump the temperature of fluid on the surface of the body will be
	//calculated.
	//Temper_of_fluid_on_the_wall_u[ij] - this is where is u velocity coordinates (x_f[i], y_v[j])
	//Temper_of_fluid_on_the_wall_v[ij] - this is where is v velocity coordinates (x_v[i], y_f[j])
	double * Temper_of_fluid_on_the_wall_u;
	double * Temper_of_fluid_on_the_wall_v;

#ifdef Calculate_p_and_Temper_using_Jakoby_method
	double * storage_Temper_2, * Temper_old;

	//For iteration mathod of Jakoby have to be made minimum two iteration on p and Temper loop
	//to be calculation process convergent. This is accepted from calculated problems.
	const int minimum_iteration_of_p_and_Temper = 2;
#else
	const int minimum_iteration_of_p_and_Temper = 1;
#endif //End of Calculate_p_and_Temper_using_Jakoby_method

	//double max_residual_in_T;

	//double max_residual_in_T_Itre0;

	unsigned int Iter_T, N_I_T;

	//MaxError_to_stop_program - if we have lower valies of
	//MaxError_Velocities, MaxError_p_c, MaxError_T and current Iter == 0
	//the program will stop, becouse that will be equivalent
	//of stationary flow and final rezult.
	//MaxError_to_stop_program == 1e-14 is excellent for 32-bit computers.
	double MaxError_to_stop_program;

	//max_time_for_solve - maximum time to solve hours. After that time program will be stoped.
	//For GRID max time foe solve is 48 hours. Give 45 hours.
	double max_time_for_solve;

	//time_solve_stop - when program must stop
	double time_solve_stop;

	//When solve start. Get the number of seconds elapsed since 00:00 hours, Jan 1, 1970 UTC from the system clock.
	int time_solve_start;


	//solve_is_finished == true  --> solve is finished
	//solve_is_finished == false  --> solve was stopped from some reson. Maybe time_solve_stop is reched.
	bool solve_is_finished;


	//double rho_0_T, rho_1_T, rho_2_T, rho_3_T, rho_4_T;
	//in points
	double p_T_0, p_T_1, p_T_2, p_T_3, p_T_4;
	//on boundaries
	double p_T_1b, p_T_2b, p_T_3b, p_T_4b;



	//dimentionless parametars
	double Ma, Pr, Fr, gamma1, c_mu, Kn, Re;
	double M_local;


	//use diferent type of non domentional volues
	unsigned int TypeOfNonDimensiolization;

	//NonDimentionalGiven_Re_pch_rho_V2_2 - that mean that system depend from
	//Reynolds number:
	//Re = Re_infinity = (rho_ch * V_ch * L / mu_ch) * (rho_non_dimentional_infinity * V_non_dimentional_infinity / mu_non_dimentional_infinity)
	//therefore: cT = (1 / Re) * (rho_non_dimentional_infinity * V_non_dimentional_infinity / mu_non_dimentional_infinity);
	//in program: cT = (1.0 / Re) * (u_given_xb * p_BC_xb / T_BC_xb);
	//and nondimensiolization for pressure is: pch = rho_ch * V_ch * V_ch / 2;
	//therefore: cPch = 0.5;
	static const unsigned int NonDimentionalGiven_Re_pch_rho_V2_2 = 1;

	//NonDimentionalGiven_Re_pch_rho_R_T - that mean that system depend from
	//Reynolds number:
	//Re = Re_infinity = (rho_ch * V_ch * L / mu_ch) * (rho_non_dimentional_infinity * V_non_dimentional_infinity / mu_non_dimentional_infinity)
	//therefore: cT = (1 / Re) * (rho_non_dimentional_infinity * V_non_dimentional_infinity / mu_non_dimentional_infinity);
	//in program: cT = (1.0 / Re) * (u_given_xb * p_BC_xb / T_BC_xb);
	//and nondimensiolization for pressure is: pch = rho_ch * R * T_ch;
	//therefore: cPch = 1 / (gamma1 * Ma * Ma);
	//where:
	//Ma = a / V_ch; Mach number, V_ch - character velocity;
	static const unsigned int NonDimentionalGiven_Re_pch_rho_R_T = 2;

	//NonDimentionalGiven_Kn_pch_rho_R_T - that mean that system depend from
	//Knudsen number:
	//Kn = l_0 / L;
	//where: l_0 - mean-free molecolar path, L - character dimention.
	//therefore: cT = c_mu * Kn / Ma;
	//V_ch = V_th - thermal velocity; V_th = sqrt(2 * R * T);
	//and nondimensiolization for pressure is: pch = rho_ch * R * T_ch;
	//therefore: cPch = 1 / (gamma1 * Ma * Ma);
	//where:
	//Ma = a / V_ch; Mach number, V_ch - character velocity;
	static const unsigned int NonDimentionalGiven_Kn_pch_rho_R_T = 3;

	//gravity coefficient
	//if(SolveGravityField) c_Fr_v = 2.0 / (gamma1 * Fr * Ma * Ma); - force from gravity field
	//else c_Fr_v = 0; - force from gravity field
	double c_Fr_v;
	bool SolveGravityField;

	//gas constant no dimention (nd)
	double u_gas_nd, v_gas_nd, V_gas_nd;

	double cPch, p_gas_nd, rho_gas_nd, T_gas_nd;



	void ApplyBodiesConditionsToVariables(void);



	inline void BoundaryConditions_Velocities(void);

	//Because of the specific of the program velocity slip BC for u are in separate function.
	//This BC can be applied right after calculating rule_vector_ij_bc_u, because
	//in this calculation are included volumes to the wall.
	//This procedure is used only when write data to output files.
	inline void BoundaryConditions_Velocities_SlipBC_u(void);

	//Vector containing the only the indexes if control volumes for field variables (others),
	//which are to the wall and have to be calculated velocity slip boundary condition for
	//horizontal component of velocity u.
	int * rule_vector_ij_BoundaryConditions_Velocities_SlipBC_u,
		* rule_vector_j_BoundaryConditions_Velocities_SlipBC_u,
		Na_BoundaryConditions_Velocities_SlipBC_u;
	//is_rule_vectors_BoundaryConditions_Velocities_SlipBC_u_defined == true  - the rule vectors rule_vector_ij_BoundaryConditions_Velocities_SlipBC_u and rule_vector_j_BoundaryConditions_Velocities_SlipBC_u are defined
	//is_rule_vectors_BoundaryConditions_Velocities_SlipBC_u_defined == false - the rule vectors rule_vector_ij_BoundaryConditions_Velocities_SlipBC_u and rule_vector_j_BoundaryConditions_Velocities_SlipBC_u are NOT defined
	bool is_rule_vectors_BoundaryConditions_Velocities_SlipBC_u_defined;


	//Because of the specific of the program velocity slip BC for v are in separate function.
	//This BC can be applied right after calculating rule_vector_ij_bc_v, because
	//in this calculation are included volumes to the wall.
	//This procedure is used only when write data to output files.
	inline void BoundaryConditions_Velocities_SlipBC_v(void);

	//Vector containing the only the indexes if control volumes for field variables (others),
	//which are to the wall and have to be calculated velocity slip boundary condition for
	//vertical component of velocity v.
	int * rule_vector_ij_BoundaryConditions_Velocities_SlipBC_v,
		* rule_vector_j_BoundaryConditions_Velocities_SlipBC_v,
		Na_BoundaryConditions_Velocities_SlipBC_v;
	//is_rule_vectors_BoundaryConditions_Velocities_SlipBC_v_defined == true  - the rule vectors rule_vector_ij_BoundaryConditions_Velocities_SlipBC_v and rule_vector_j_BoundaryConditions_Velocities_SlipBC_v are defined
	//is_rule_vectors_BoundaryConditions_Velocities_SlipBC_v_defined == false - the rule vectors rule_vector_ij_BoundaryConditions_Velocities_SlipBC_v and rule_vector_j_BoundaryConditions_Velocities_SlipBC_v are NOT defined
	bool is_rule_vectors_BoundaryConditions_Velocities_SlipBC_v_defined;


	//boundary conditions for u dudx
	//dudx_0_BC_xb give: du/dx = 0, on x = x_b
	//dudx_0_BC_xe give: du/dx = 0, on x = x_e
	bool dudx_0_BC_xb, dudx_0_BC_xe;
	//durhodx_0_BC_xb - that give: d(u * rho) / dx = 0, on x = x_b;
	//durhodx_0_BC_xe - that give: d(u * rho) / dx = 0, on x = x_b;
	bool durhodx_0_BC_xb, durhodx_0_BC_xe;
	//drhodt_durhodx_dvrhody_0_BC_xe - that give: drho/dt + d(u * rho)/dx + d(v * rho)/dy = 0; - continuity equation over the x_e boundary
	//drhodt_durhodx_dvrhody_0_BC_xb - that give: drho/dt + d(u * rho)/dx + d(v * rho)/dy = 0; - continuity equation over the x_b boundary
	bool drhodt_durhodx_dvrhody_0_BC_xb, drhodt_durhodx_dvrhody_0_BC_xe;
	//u_uMin_Mout_0_BC_xe give: u[i1j] - u[ij] * Min / Mout = 0, on x = x_e
	//--> u[i1j] - u[ij] * Min / Mout, on x = x_e
	//where:
	//	Min - mass flux coming into domain
	//	Mout - mass flux leaveing domain
	bool u_uMin_Mout_0_BC_xe;
	double Min, Mout;

	//This is calculation u[i_1j] from equation for u for node ij.
	//The target for this kind of calculation is to keep the momentum on OX
	//at the x = x_b.
	bool ToCalculate_ui_1j_from_equation_for_uij_at_x_b;

	//This is calculation u[i1j] from equation for u for node ij.
	//The target for this kind of calculation is to keep the momentum on OX
	//at the x = x_e.
	bool ToCalculate_ui1j_from_equation_for_uij_at_x_e;


	//dudy_0_BC_yb give: du/dy = 0, on y = y_b
	//dudy_0_BC_ye give: du/dy = 0, on y = y_e
	bool dudy_0_BC_yb, dudy_0_BC_ye;

	//boundary conditions for v
	//dvdx_0_BC_xb give: dv/dx = 0, on x = x_b
	//dvdx_0_BC_xe give: dv/dx = 0, on x = x_e
	bool dvdx_0_BC_xb, dvdx_0_BC_xe;
	//dvdy_0_BC_xb give: dv/dy = 0, on y = y_b
	//dvdy_0_BC_xe give: dv/dy = 0, on y = y_e
	bool dvdy_0_BC_yb, dvdy_0_BC_ye;

	//This is calculation v[i_1j] from equation for v for node ij.
	//The target for this kind of calculation is to keep the momentum on OY
	//at the x = x_b.
	bool ToCalculate_vi_1j_from_equation_for_vij_at_x_b;

	//This is calculation v[i1j] from equation for v for node ij.
	//The target for this kind of calculation is to keep the momentum on OY
	//at the x = x_e.
	bool ToCalculate_vi1j_from_equation_for_vij_at_x_e;


	inline void BoundaryConditions_Pressure(void);
	inline void BoundaryConditions_Pressure_begin_iteration_for_it(void);
	//boundary conditions for pressure
	//dpdx_0_BC_xb give dp / dx = 0, on x = x_b
	//dpdx_0_BC_xe give dp / dx = 0, on x = x_e
	//dpdy_0_BC_yb give dp / dy = 0, on y = y_b
	//dpdy_0_BC_ye give dp / dy = 0, on y = y_e
	bool dpdx_0_BC_xb, dpdx_0_BC_xe, dpdy_0_BC_yb, dpdy_0_BC_ye;

	//This is calculation p[i_1j] from equation for p for node ij.
	//The target for this kind of calculation is to keep the mass
	//at the x = x_b.
	bool ToCalculate_pi_1j_from_equation_for_pij_at_x_b;

	//This is calculation p[i1j] from equation for p for node ij.
	//The target for this kind of calculation is to keep the mass
	//at the x = x_e.
	bool ToCalculate_pi1j_from_equation_for_pij_at_x_e;


	inline void BoundaryConditions_Temperature(void);

	//Because of the specific of the program temperature jump BC are in separate function.
	//This BC can be applied right after calculating rule_vector_ij_bc_others, because
	//in this calculation are included volumes to the wall.
	//This procedure is used only when write data to output files.
	inline void BoundaryConditions_Temperature_TemperatureJumpBC(void);

	//Vector containing the only the indexes if control volumes for field variables (others),
	//which are to the wall and have to be calculated temperature jump boundary condition.
	int * rule_vector_ij_BoundaryConditions_Temperature_TemperatureJumpBC,
		* rule_vector_j_BoundaryConditions_Temperature_TemperatureJumpBC,
		Na_BoundaryConditions_Temperature_TemperatureJumpBC;
	//is_rule_vectors_BoundaryConditions_Temperature_TemperatureJumpBC_defined == true  - the rule vectors rule_vector_ij_BoundaryConditions_Temperature_TemperatureJumpBC and rule_vector_j_BoundaryConditions_Temperature_TemperatureJumpBC are defined
	//is_rule_vectors_BoundaryConditions_Temperature_TemperatureJumpBC_defined == false - the rule vectors rule_vector_ij_BoundaryConditions_Temperature_TemperatureJumpBC and rule_vector_j_BoundaryConditions_Temperature_TemperatureJumpBC are NOT defined
	bool is_rule_vectors_BoundaryConditions_Temperature_TemperatureJumpBC_defined;


	//This procedure write temperature to the Temper matrix before write data to file.
	//The target is to have one Temperature of fluid on teh wall in Temper.
	//This procedure can NOT be used to calculation process, because on external vertex
	//is impossible to have two different values in one point.
	inline void BoundaryConditions_Temperature_TemperatureJumpBC_WriteToMassiveToWriteDataToFile(void);

	//boundary conditions for Temperature
	bool dTdx_0_BC_xb, dTdx_0_BC_xe, dTdy_0_BC_yb, dTdy_0_BC_ye;

	//This is calculation Temper[i_1j] from equation for Temper for node ij.
	//The target for this kind of calculation is to keep the energy
	//at the x = x_b.
	bool ToCalculate_Temperi_1j_from_equation_for_Temperij_at_x_b;

	//This is calculation Temper[i1j] from equation for Temper for node ij.
	//The target for this kind of calculation is to keep the energy
	//at the x = x_e.
	bool ToCalculate_Temperi1j_from_equation_for_Temperij_at_x_e;


	//Use Kn local on the wall BC.
	//This are BC for velocity slip and temperature jump.
	//Kn_lical = Kn * rho_local.
	//rho_local is local density, and for tmperature jump is rho[ij]
	//for u is rho_u[ij]
	//for v is rho_v[ij]
	bool To_Use_Kn_local_in_wall_BC;


	bool Periodic_boundary_conditions_about_OX;

	bool Given_velocity_on_xb, Given_velocity_on_xe;

	//When on boundaries x_b and x_e is givven pressure
	bool Pressure_BC;
	//Boundary Conditions for given pressure on bpundary
	double p_BC_xb, p_BC_xe, T_BC_xb, T_BC_xe, rho_BC_xb, rho_BC_xe;

	//p_BC_xb_correction_method == 0 --> no correction for p_xb
	unsigned int p_BC_xb_correction_method;
	//p_BC_xb_correction_method == p_BC_xb_correction_method_Vmax_xb --> V maximum inflow == u_given_xb
	static const unsigned int p_BC_xb_correction_method_Vmax_xb = 1;
	//p_BC_xb_correction_method == p_BC_xb_correction_method_Vmean_xb --> V mean inflow == u_given_xb
	static const unsigned int p_BC_xb_correction_method_Vmean_xb = 2;
	//Pressure_ratio_correct == true - correct pressure ratio
	//Pressure_ratio_correct == false - NOT correct pressure ratio
	bool Pressure_ratio_correct;

	//correct_p_BC_xb == true - correct pressure on x_b
	//correct_p_BC_xb == false - correct pressure on x_e
	bool correct_p_BC_xb;
	//u_given_xb - the velocity, which is given when we have pressure boundary conditions and is given on x = x_b
	double u_given_xb;
	//u_given_xb_in_it - this is the value of velocity in time it on x_b
	double u_given_xb_in_it;
	//u_given_xb_in_it_max - this is the maximum value of velocity in time it on x_b
	double u_given_xb_in_it_max;
	//u_given_xb_in_it_mean - this is the mean value of velocity in time it on x_b
	double u_given_xb_in_it_mean;
	//maximum error from given velocity on x_b
	double u_given_xb_error;
	//p_correction_auto == true - correct p_correction automaticly
	//p_correction_auto == false - p_correction == const
	bool p_correction_auto;
	//minimum and maximum possible values of p_correction
	double p_correction_min, p_correction_max;
	//correct given pressure on x_b or x_e
	double p_correction;

	//ToCheckForFileIsNow_p_u_x_b == true - for first if file exist and the solve is begin will delete exist file
	//ToCheckForFileIsNow_p_u_x_b == false - if solve is begin file will not be deleted
	bool ToCheckForFileIsNow_p_u_x_b;



	//Variables and function for Solve and Write Drug Coefficient
	void WriteDragCoefficient(void);

	//ToCheckForFileIsNow_CD == true - for first if file exist and the solve is begin will delete exist file
	//ToCheckForFileIsNow_CD == false - if solve is begin file will not be deleted
	bool ToCheckForFileIsNow_CD;

	//variable which are needed to solve CD (DragCoefficient)
	unsigned int counter_bottom , counter_front, counter_top, counter_behind,
		Nt_save_DragCoefficient_counter, Nt_save_DragCoefficient;

	double CD_fr_x, CD_p_x,
		CD_fr_bottom, CD_p_bottom,
		CD_fr_front, CD_p_front,
		CD_fr_top, CD_p_top,
		CD_fr_behind, CD_p_behind;
	double S_for_CD_bottom, S_for_CD_front, S_for_CD_top, S_for_CD_behind;

	//test source term
	double Su_abs_max, Sv_abs_max;


	template <typename T_number>
	inline bool positive(T_number& number)
	//Recieve number and return true if number is positive or false if number is 0 or negative.
	{
		if(number > 0) return(true);
		else return(false);
	}


	//Thios procedure solve rho in midlle points - slolve rho_u and rho_v
	inline void Solve_rho_in_middle_points(void);

	//procedure solve Fuxes Fx and Fy
	inline void SolveFluxes(void);


	//solve Continuity Equation
	bool ToSolveContinuityEquation;
	double * ContinuityEquation_stationary, * ContinuityEquation;
	inline void SolveContinuityEquation(void);


	//Number of circle
	int N_circle;
	//Number of polyhedron
	int N_polyhedron;


	//vector contain max residual in p and Temper for last Iter
	double * max_residual_in_p, * max_residual_in_T;

	//vector contain max resudual in p and Temper for first (zero) Iter
	double * max_residual_in_p_Iter_0, * max_residual_in_T_Iter_0;

	//vector contain number of iterations to solve p and Temper for every Iter
	unsigned int * max_Iter_p_c;

	//This variable define to continue iterations ot not
	//to_continue_iterations == true  - iterations contunue
	//to_continue_iterations == false  - iterations stop
	bool to_continue_iterations;


	//weight coefficient for SOR, 1 <= w <= 2
	//Aproximetry:
	//w_u == w_v == 1.11
	//w_p == w_Temper == 1.3
	double w_u, w_v, w_p, w_Temper;

	//residual in equation:
	double r_u, r_v, r_p, r_Temper;

	bool ToStartFromExactSolutionForChannelFlowWithPressureBC;
	//begin index of channel, end index of channel
	unsigned int j_channel_b, j_channel_e;

	double dpdx_PressureBC, L_channel, h_channel;
	unsigned int body_SlipBC_for_exact_solution;


	inline void SolveDiffisionTerms(void);



	//In solving problem we give Reynolds (Re)
	//and Knudsen (Kn) numbers
	//u_given_xb = Re * Kn * sqrt(15.0 * M_PI / 128.0);
	bool GivenReKn;



	//Define MPI variables, functions and procedures

	//Variables, which are defined to simplify the code
	//Number of SubDomains in OX nad OY directions, respectively
	unsigned int N_SubDomains_x, N_SubDomains_y;


	//Variable, which rule the way of write and read data, when is made
	//computation with MPI.
	//
	//ToReadSolvedDataFromSeparateFileForEachSubDomain == true, read data from file for
	//each subdomain separately. For example for u and IJ == 0 and binary file
	//read u from file u.0.bin.
	//ToReadSolvedDataFromSeparateFileForEachSubDomain == false, read data from file for entire
	//computational domain and after that take data for this subdomain.
	bool ToReadSolvedDataFromSeparateFileForEachSubDomain;

	/*
	ToWriteSolvedDataToFileForEntireComputationalDomain == true, write data for u, v, p, Temper
		and vol_inf_fb to one file for u, v, p, Temper and vol_inf_fb, respectively. This file
		is data for entire computational domain.
	ToWriteSolvedDataToFileForEntireComputationalDomain == false, not write data for u, v, p, Temper
		and vol_inf_fb to one file for u, v, p, Temper and vol_inf_fb, respectively.
	*/
	bool ToWriteSolvedDataToFileForEntireComputationalDomain;

	/*
	ToWriteSolvedDataToSeparateFileForEachSubDomain == true, write data to file for
		each subdomain separately and for all computational domain.
		For example for u and IJ == 0 and binary file
		write u to file u.0.bin and to u.bin.
	ToWriteSolvedDataToSeparateFileForEachSubDomain == false, not write data to file
		for each subdomain in files.
	*/
	bool ToWriteSolvedDataToSeparateFileForEachSubDomain;

	//This is option to write only number of subdomains if there is need to make change.
	//The program have to be made to start from begining just for this this change.
	bool ToOnlyWriteDataFor_Nx_Ny_ForAllSubDomains_AndExit;


	//Define BlockLength of swap area in different directions.
	//Block Length of swap area in direction I_1. This is 2,
	//because rho_u[i = 1, j] need rho[i = 0, j] and rho[i = 1, j] to be calculated.
	const unsigned int BlockLengthSwap_I_1 = 2;

	//Block Length of swap area in direction I1. This is 1.
	const unsigned int BlockLengthSwap_I1 = 1;

	//Block Length of swap area in direction J_1. This is 2,
	//because rho_v[i, j = 1] need rho[i, j = 0] and rho[i, j = 1] to be calculated.
	const unsigned int BlockLengthSwap_J_1 = 2;

	//Block Length of swap area in direction J1. This is 1.
	const unsigned int BlockLengthSwap_J1 = 1;
	//End of Define BlockLength of swap area in different directions.


#ifdef COMPILE_MPI


	//Number of domains
	int I, J, IJ, I_1J, I1J, IJ_1, IJ1;
	int N_SubDomains_a;

	//rank - rank of process, it is initialized after MPI_Init.
	//Number_of_processes - Number of processes, it is initialized after MPI_Init.
	int rank, Number_of_processes;

	//N_neighbouhoods - number of neighbouhoods with wich exchange data
	//		- for 1D N_neighbouhoods = 2
	//		- for 2D N_neighbouhoods = 4
	//		- for 3D N_neighbouhoods = 6
	//This is 2D case therefore
	const unsigned int N_neighbouhoods = 4;

	//N_exchanged_variables_with_neighbouhoods - number of axchanged variabled with neighbouhoods
	//		- in this case N_exchanged_variables = 10
	//			(d_p_c_12, apx, bpx, d_p_c_34, apy, bpy, u, v, p, Temper)
	//
	//N_send_variables_to_IJ0 - number of send variables to process IJ == 0
	//		- in this case N_send_variables_to_IJ0 = 4
	//			(max_residual_in_p, max_residual_in_T, max_residual_in_u, max_residual_in_v)
	//
	//N_send_from_IJ0 - number of send variables from process IJ == 0 to all SubDomains.
	//		- in this case N_send_from_IJ0 = 3
	//			(continue_iter_int, countinue_time_step_int)
	const unsigned int N_exchanged_variables_with_neighbouhoods = 10;
	const unsigned int N_send_variables_to_IJ0 = 4;
	const unsigned int N_send_from_IJ0 = 3;

	//Number of nodes in entire computational domain.
	unsigned int Nx_entire_computational_domain, Ny_entire_computational_domain, Na_entire_computational_domain;

	//variables for mesh for entire computational domain
	bool is_defined_array_for_mesh_entire_computational_domain;
	double * x_f_entire_computational_domain, * y_f_entire_computational_domain,
		* hx_entire_computational_domain, * hy_entire_computational_domain,
		* x_v_entire_computational_domain, * y_v_entire_computational_domain;


	//Information for all SubDomans
	DomainDecomposition2D * SubDomain;

	//Maximum request and status variables that are needed are
	//N_max_requests = N_neighbouhoods * N_exchanged_variables * N_of_directions
	//Where:
	//N_of_directions - number of direction for data exchange.
	//In this case data from onde subdomain IJ is send to subdomain I_1J
	//and after that IJ recv data from I_1J.
	//That mean that N_of_directions = 2, for this case.
	const unsigned int N_of_directions = 2;
	unsigned int N_max_requests;


	/*
	For every send variable (for recv variable the tag i sthe same as send, this is
	the meaning) we will have uniqe tag.
	Positives:			reduce the number of possible errors
						give possobility to enlarge the time of echange data

	Negatives:			use more memory. This have to be check, when program is
						ruuning on large naumbers of cpus.

	Send variables are:	d_p_c_12, apx, bpx, d_p_c_34, apy, bpy, u, v, p, Temper,
						max_residual_in_p, max_residual_in_T, max_residual_in_u,
						max_residual_in_v, continue_iter_int, countinue_time_step_int.
						Number of all variables are 16.
						10 to all neighboors, 4 to process IJ == 0 and 2 from
						process IJ == 0 to all SubDomains.
	*/
	//Index of tag according to send variable and neghbour:
	//This variables are 10, therefore N_exchanged_variables_with_neighbouhoods = 10
	//d_p_c_12_send - send d_p_c_12 to neighbour
	const int d_p_c_12_send = 0;
	//apx_send - send apx to neighbour
	const int apx_send = 1;
	//bpx_send - send bpx to neighbour
	const int bpx_send = 2;

	//d_p_c_34_send - send d_p_c_34 to neighbour
	const int d_p_c_34_send = 3;
	//apy_send - send apy to neighbour
	const int apy_send = 4;
	//bpy_send - send bpy to neighbour
	const int bpy_send = 5;

	//p_send - send p to neighbour
	const int p_send = 6;
	//Temper_send - send Temper to neighbour
	const int Temper_send = 7;
	//u_send - send u to neighbour
	const int u_send = 8;
	//v_send - send v to neighbour
	const int v_send = 9;


	//Variable, which send information to process IJ == 0.
	//This variables are 4, therefore N_send_variables_to_IJ0 = 4
	//max_residual_in_p_send_to_IJ0 - send max_residual_in_p to process IJ == 0
	const int max_residual_in_p_send_to_IJ0 = 0;

	//max_residual_in_T_send_to_IJ0 - send max_residual_in_T to process IJ == 0
	const int max_residual_in_T_send_to_IJ0 = 1;

	//max_residual_in_u_send_to_IJ0 - send max_residual_in_u to process IJ == 0
	const int max_residual_in_u_send_to_IJ0 = 2;

	//max_residual_in_v_send_to_IJ0 - send max_residual_in_v to process IJ == 0
	const int max_residual_in_v_send_to_IJ0 = 3;


	//Variable, which send information tfrom IJ == 0 to every one process (SubDomain).
	//This variables are 3, therefore N_send_from_IJ0 = 3
	//continue_iter_int_send_from_IJ0 - send continue_iter_int from process IJ == 0
	//to all SubDomains
	const int continue_iter_int_send_from_IJ0 = 0;

	//countinue_time_step_int_send_from_IJ0 - send countinue_time_step_int from process IJ == 0
	//to all SubDomains
	const int continue_time_step_int_send_from_IJ0 = 1;

	//continue_iter_p_c_int_send_from_IJ0 - send continue_iter_p_c_int from process IJ == 0
	//to all SubDomains
	const int continue_iter_p_c_int_send_from_IJ0 = 2;


	//Send variable to neighbours:
	const int to_I_1J = 0;
	const int to_I1J = 1;
	const int to_IJ_1 = 2;
	const int to_IJ1 = 3;


	//Define index variables rule the requests and status arraies
	//Send or Recv index:
	const int index_send = 0;
	const int index_recv = 1;
	//Number of possible actions. For this case send and recv, this are 2.
	const int N_actions = 2;

	/*
	Target:		Calculate index for request for send/recv variable to/from neighbour.
				This function return value for index, which is calculated using equation
				in depend of variables, neighbour and send_or_recv_data_tmp. The equation
				just arrange the indexes. This index is unique.
	Receive:	variable_tmp - variable which will be send. The variable for, which can
				be used this function are: d_p_c_12_send, apx_send, bpx_send,
				d_p_c_34_send, apy_send, bpy_send, p_send, Temper_send, u_send and
				v_send.

				neighbour_tmp - neighbour to which will be send variable. The variables
				are to_I_1J, to_I1J, to_IJ_1 and to_IJ1.

				send_or_recv_data_tmp - the process is send or recv of data. The variables
				are index_send and index_recv.
	*/
	inline int CalculateIndexRequestForSendOrRecvVariableNeighbour(
		const int& variable_tmp, const int& neighbour_tmp, const int& send_or_recv_data_tmp)
	{
		return(variable_tmp * N_actions * N_neighbouhoods
			+ send_or_recv_data_tmp * N_neighbouhoods
			+ neighbour_tmp);
	};


	/*
	Target:		Calculate index for request for send variable to IJ == 0.
				This function return value for index, which is calculated using equation
				in depend of variables. The equation just arrange the indexes. This index
				is unique.
	Receive:	variable_tmp - variable which will be send. The variable for which can
				be used this function are: max_residual_in_p_send_to_IJ0,
				max_residual_in_T_send_to_IJ0, max_residual_in_u_send_to_IJ0 and
				max_residual_in_v_send_to_IJ0.
	*/
	inline int CalculateIndexRequestForSendVariableToIJ0(
		const int& variable_tmp)
	{
		return(N_exchanged_variables_with_neighbouhoods * N_actions * N_neighbouhoods
			+ variable_tmp);
	};


	/*
	Target:		Calculate index for request for send variable from process IJ == 0
				to all processes (subdomains). This function return value for index
				for request, which is calculated using equation in depend of variable
				The equation just arrange the indexes. This index is unique.
	Receive:	variable_tmp - variable which will be send. The variable for which can
				be used this function are: continue_iter_int_send_from_IJ0 and
				countinue_time_step_int_send_from_IJ0.
	*/
	inline int CalculateIndexRequestForRecvVariableFromIJ0(
		const int& variable_tmp)
	{
		return(N_exchanged_variables_with_neighbouhoods * N_actions * N_neighbouhoods
			+ N_send_variables_to_IJ0
			+ variable_tmp);
	};



	//Variables for MPI request
	MPI_Request * request;

	//Variablee for report status of MPI data exchange
	MPI_Status * status;


	//Variables for MPI request for max_residual_in_u iteration criterion.
	MPI_Request * request_max_residual_in_u;
	//Variablee for report status of MPI data exchange for max_residual_in_u iteration criterion.
	MPI_Status * status_max_residual_in_u;

	//Variables for MPI request for max_residual_in_v iteration criterion.
	MPI_Request * request_max_residual_in_v;
	//Variablee for report status of MPI data exchange for max_residual_in_v iteration criterion.
	MPI_Status * status_max_residual_in_v;

	//Variables for MPI request for max_residual_in_p iteration criterion.
	MPI_Request * request_max_residual_in_p;
	//Variablee for report status of MPI data exchange for max_residual_in_p iteration criterion.
	MPI_Status * status_max_residual_in_p;

	//Variables for MPI request for max_residual_in_T iteration criterion.
	MPI_Request * request_max_residual_in_T;
	//Variablee for report status of MPI data exchange for max_residual_in_T iteration criterion.
	MPI_Status * status_max_residual_in_T;



	//Here have to be created Derived Datatipe to transfer the data.
	//Create derived datatype MPI_TYPE_VECTOR to send data to I_1J
	MPI_Datatype send_to_I_1J;

	//Create derived datatype MPI_TYPE_VECTOR to receive data from I_1J
	MPI_Datatype recv_from_I_1J;

	//Here have to be created Derived Datatipe to transfer the data.
	//Create derived datatype MPI_TYPE_VECTOR to send data to I1J
	MPI_Datatype send_to_I1J;

	//Create derived datatype MPI_TYPE_VECTOR to receive data from I1J
	MPI_Datatype recv_from_I1J;

	//This variable make creation of above MPI_Datatype send_to_I_1J,
	//recv_from_I_1J, send_to_I1J, and recv_from_I1J to be created once, because
	//this MPI_Datatype are the same evry rime, when are needed.
	//This make calclulations much faster
	//is_created_MPI_Datatype_send_to_I_1J == false - this mean that MPI_Datatype send_to_I_1J
	//is not created MUST be created to make data exchange.
	//is_created_MPI_Datatype_send_to_I_1J == true - this mean that MPI_Datatype send_to_I_1J
	//is created and is not necessury to be created again.
	bool is_created_MPI_Datatype_send_to_I_1J, is_created_MPI_Datatype_recv_from_I_1J;
	bool is_created_MPI_Datatype_send_to_I1J, is_created_MPI_Datatype_recv_from_I1J;

	//Number of elements, which will be send to subdomain I_1J
	unsigned int N_send_to_I_1J;

	//Blocklength send to I_1J
	unsigned int blocklength_send_to_I_1J;

	//Stride send to I_1J
	unsigned int stride_send_to_I_1J;


	//Number of elements, which will be receive from subdomain I_1J
	unsigned int N_recv_from_I_1J;

	//Blocklength receive to I_1J
	unsigned int blocklength_recv_from_I_1J;

	//Stride receive to I_1J
	unsigned int stride_recv_from_I_1J;


	//Number of elements, which will be send to subdomain I1J
	unsigned int N_send_to_I1J;

	//Blocklength send to I1J
	unsigned int blocklength_send_to_I1J;

	//Stride send to I1J
	unsigned int stride_send_to_I1J;


	//Number of elements, which will be recv from subdomain I1J
	unsigned int N_recv_from_I1J;

	//Blocklength receive to I1J
	unsigned int blocklength_recv_from_I1J;

	//Stride receive to I1J
	unsigned int stride_recv_from_I1J;


	//Number of elements, which will be send to subdomain IJ_1
	unsigned int N_send_to_IJ_1;

	//Number of elements, which will be recv from subdomain IJ_1
	unsigned int N_recv_from_IJ_1;


	//Number of elements, which will be send to subdomain IJ1
	unsigned int N_send_to_IJ1;

	//Number of elements, which will be recv from subdomain IJ1
	unsigned int N_recv_from_IJ1;


	////All procedures and functions according to the MPI end on "_MPI".
	//Define SubDomains - number of nodes,
	void DefineSubDomains_MPI(void);



	/*
	Target:		Start data exchange with neighbourhoods subdomains.
	Receive:	massive - array to exchange (fi0 or fi1)
				Nx_massive_limits - SubDomain Computation Limits on OX (Nx_u, Nx_v or Nx_others)
				Ny_massive_limits - SubDomain Computation Limits on OY (Ny_u, Ny_v or Ny_others)
	*/
	template <typename T_massive, typename T_Nx_limits_massive, typename T_Ny_limits_massive>
	inline void StartDataExchangeWithNeighbourhoodsSubDomains_MPI(
		T_massive massive[],
		T_Nx_limits_massive Nx_massive_limits[],
		T_Ny_limits_massive Ny_massive_limits[]);


	/*
	Target:		Wait to compete data exchange with neighbourhoods subdomains.
	Receive:	massive - massive to exchange
	*/
	template <typename T_massive>
	inline void WaitToCompleteDataExchangeWithNeighbourhoodsSubDomains_MPI(
		T_massive massive[]);


	//Variables needed to exchange data for convergence criterions.
	//It is made this organization to be used time for data exchange for calculation.
	/*
	Target:		Start exchange data for convergence criterions - max_residual_in_p
				and max_residual_in_T for this subdomain to process IJ == 0.
	Note:		This procedure must be executed, when no other exchange data proceed,
				because in it are used exist tags.
	*/
	inline void StartSendDataForConvergenceCriterionsFromSubDomainsFor_pT_MPI(void);

	/*
	Target:		Start exchange data for convergence criterions - Send max_residual_in_u
				for this subdomain to process IJ == 0.
	Note:		This procedure must be executed, when no other exchange data proceed,
				because in it are used exist tags.
	*/
	inline void StartSendDataForConvergenceCriterionsFromSubDomainsFor_u_MPI(void);

	/*
	Target:		Start exchange data for convergence criterions - Send max_residual_in_v
				for this subdomain to process IJ == 0.
	*/
	inline void StartSendDataForConvergenceCriterionsFromSubDomainsFor_v_MPI(void);



	//Make temporal array where to keep data for convergence criterians
	//from all subdomains.
	double * max_residual_in_u_all_subdomains, * max_residual_in_v_all_subdomains,
		 * max_residual_in_p_all_subdomains, * max_residual_in_T_all_subdomains;

	/*
	Target:		Wait to complete recv data for convergence criterions from subdomains for p and Temper
	*/
	inline void WaitToCompleteRecvDataForConvergenceCriterionsFromSubDomains_for_pT_MPI(void);

	/*
	Target:		Wait to complete recv data for convergence criterions from subdomains for u and v
	*/
	inline void WaitToCompleteRecvDataForConvergenceCriterionsFromSubDomains_for_uv_MPI(void);

	//Find maximum resuduals from all subdomains for p and Temper.
	void FindMaximumResudualsFromAllSubDomains_for_pT_MPI(void);

	//Find maximum resuduals from all subdomains for u and v.
	void FindMaximumResudualsFromAllSubDomains_for_uv_MPI(void);

	/*
	Target:		Start send data for convergence criterions to subdomains.
				When data are received in array continue_iter_all_subdomains,
				the all data are checked and if in all subdomain the convergence
				criterions are satisfyed, is send commant to stop calculating this
				time step. But if in one of all subdomains we do not have satisfaction
				of convergence criterions the calculation is countinue.

	Note:		This procedure must be executed, when no other exchange data proceed,
				because in it are used exist tags.
	*/
	inline void StartSendDataForConvergenceCriterionsToSubDomains_for_pT_MPI(void);


	/*
	Target:		Wait to complete recv data for convergence criterions for p and Temper loop
				from IJ == 0.
	*/
	inline void WaitToCompleteRecvDataForConvergenceCriterionsFromIJ_0_for_pT_MPI(void);

	/*
	Target:		Start send data for convergence criterions to subdomains.
				When data are received in array continue_iter_all_subdomains,
				the all data are checked and if in all subdomain the convergence
				criterions are satisfyed, is send commant to stop calculating this
				time step. But if in one of all subdomains we do not have satisfaction
				of convergence criterions the calculation is countinue.

	Note:		This procedure must be executed, when no other exchange data proceed,
				because in it are used exist tags.
	*/
	inline void StartSendDataForConvergenceCriterionsToSubDomains_for_it_loop_MPI(void);


	/*
	Target:		Wait to complete recv data for convergence criterions from IJ == 0.
	*/
	inline void WaitToCompleteRecvDataForConvergenceCriterionsFromIJ_0_for_it_loop_MPI(void);


	//Write solved data to files for entire computational subdomain:
	/*
	Target:		Write array collected from all subdomains to file.
	Receive:	FileNameWithoutExtension - name of file where data will be written.
					Here is not included extesion, because the extension depend
					from kind of file: binary or text, and extension is putted
					in the procedure.
				massive_entire_computational_domain - massive, which will be written
					to file.
				massive_type - type of massive,
					if massive is from double, then massive_type == double_type
					if massive is from unsigned int, then massive_type == unsigned_int_type
	Note:		If in program is needed to be exchaneged data from other type,
				which is not defined now, the type have to be defined in
				similar way like existing.
	*/
	template <typename T_FileName, typename T_massive>
	void Write_massive_entire_computational_domain_MPI(
		const T_FileName& FileNameWithoutExtension,
		T_massive massive[], const unsigned int& massive_type);

	//Type of data to exchange with MPI.
	//If There is needed of other type of data to be exchanged have to define
	//onother variable.
#define double_type 1
#define unsigned_int_type 2

	//Read solved data from files for entire computational subdomain:
	/*
	Target:		Read array collected from all subdomains to file.
	Receive:	FileNameWithoutExtension - name of file from where data will be readed.
					Here is not included extesion, because the extension depend
					from kind of file: binary or text, and extension is putted
					in the procedure.
				massive_entire_computational_domain - massive, which will be readed from file.
				massive_type - type of massive,
					if massive is from double, then massive_type == double_type
					if massive is from unsigned int, then massive_type == unsigned_int_type
	Note:		If in program is needed to be exchaneged data from other type,
				which is not defined now, the tipe have to be defined in
				similar way like existing.
	*/
	template <typename T_FileName, typename T_massive>
	void Read_massive_entire_computational_domain_MPI(
		const T_FileName& FileNameWithoutExtension,
		T_massive massive[], const unsigned int& massive_type);


	//Write arrays containitng Nx and Ny, respectively, for all subdomains.
	void WriteDataFor_Nx_Ny_ForAllSubDomains(void);

	//Read arrays containitng Nx and Ny, respectively, for all subdomains.
	void ReadDataFor_Nx_Ny_ForAllSubDomains(void);

	//Arrays containitng Nx and Ny, respectively for all subdomains.
	int * Nx_all_SubDomains, * Ny_all_SubDomains;

	bool is_defined_arrays_for_Nx_Ny_for_all_SubDomains;


	//This variables are int, because in MPI is no to send bool type.
	int continue_iter_int, countinue_time_step_int, continue_iter_p_c_int;


	/*
	Target:		Collect data for CD (Drag Coefficient), from all SubDomains (0 < IJ)
				to SubDomain (IJ == 0). The process IJ == 0, will write data for CD to file.
	*/
	void CollectDataForCDfromAllSubDomainsToProcessIJ0(void);




#endif //end of COMPILE_MPI
	//End of Define MPI variables


int main(int argc, char* argv[])
{
	//define times for start of program
	time_solve_start = time(NULL);


	//initiate variables
	is_rule_vectors_defined = false;
	is_defined_array_for_mesh = false;

	is_rule_vectors_BoundaryConditions_Velocities_SlipBC_u_defined = false;
	is_rule_vectors_BoundaryConditions_Velocities_SlipBC_v_defined = false;

	is_rule_vectors_BoundaryConditions_Temperature_TemperatureJumpBC_defined = false;


#ifdef COMPILE_MPI
	//initial variables	according to MPI
	is_defined_arrays_for_Nx_Ny_for_all_SubDomains = false;

	is_defined_array_for_mesh_entire_computational_domain = false;
	is_created_MPI_Datatype_send_to_I_1J = false;
	is_created_MPI_Datatype_recv_from_I_1J = false;
	is_created_MPI_Datatype_send_to_I1J = false;
	is_created_MPI_Datatype_recv_from_I1J = false;


	MPI_Init(&argc, &argv);

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &Number_of_processes);

#endif //End of COMPILE_MPI

	//SIMPLEST_ANL_compr_T_h_var_NoDim_Algorithm  Solver;
	//Solver.Solve();

//void Solve(void)
//{
	bool ToReadEnteredDataFromFile;
	ToReadEnteredDataFromFile = //false;
								true;

	if(ToReadEnteredDataFromFile)
	{
		ReadEnteredDataFromFile();
	}
	else
	{
		//Body for solving - rectangular
		u_gas_nd = 0.25;
		v_gas_nd = 0;
		T_gas_nd = 1;
		p_gas_nd = 1;

		Ma = sqrt(6.0 / 5.0);
		Pr = 2.0 / 3.0;
		gamma1 = 5.0 / 3.0;
		c_mu = (5.0 / 16.0) * sqrt(2.0 * M_PI / gamma1);
		Kn = 0.01;
		Fr = 50;
		SolveGravityField = false;
		Re = 100;
		TypeOfNonDimensiolization = 3;

		Nt = 250000;
		Nx = 200;
		Ny = 150;
		N_I = 1;
		N_I_p_c = 2;
		N_I_T = 2;

		ht = 1e-2;


		x_b = 0;
		x_e = 41;

		y_b = 0;
		y_e = 25;


		kind_of_mesh = 5;

		hxmin = 1.5e-2;
		hxmax = 1.5e-1;
		xfmin = 12.0;
		xfmax = 11.0;
		xbmin = 13.0;
		xbmax = 14.0;

		hymin = 1.5e-2;
		hymax = 1.5e-1;
		ytmin = 13.0;
		ytmax = 14.0;
		ybmin = 12.0;
		ybmax = 11.0;


		Nt_save_solved_data = 1000;
		Nt_save_DragCoefficient = 250;
		Nt_save_temp_data = 250;
		N_I_check = 1;
		WriteDataInShortFormat_kind = 0;

		MaxError_Velocities = 1e-8;
		MaxError_p_c = 1e-9;
		MaxError_T = 1e-9;
		MaxError_to_stop_program = 1e-12;
		max_time_for_solve = 10000;
		solve_is_finished = false;


		ToReadSolvedDataFromFile = false;
		ToReadSolvedDataFromBinaryFile = true;
		ToWriteSolvedDataToBinaryFile = true;
		ToWriteSolvedDataToNewFiles = false;

		ToImport_data_for_circles_from_file_b = false;
		ToImport_data_for_polyhedrons_from_file_b = true;

		ToImport_data_for_circles_from_file_p = false;
		ToImport_data_for_polyhedrons_from_file_p = false;

		ToImport_data_for_circles_from_file_V = false;
		ToImport_data_for_polyhedrons_from_file_V = false;


		Given_velocity_on_xb = false;
		Given_velocity_on_xe = false;
		Periodic_boundary_conditions_about_OX = false;
		Pressure_BC = false;
		ToStartFromExactSolutionForChannelFlowWithPressureBC = true;


		p_BC_xb_correction_method = 0;
		Pressure_ratio_correct = 0;
		GivenReKn = 1;
		correct_p_BC_xb = 0;
		u_given_xb = 0.2;
		u_given_xb_error = 1e-4;
		p_correction_auto = 1;
		p_correction_min = 1e-5;
		p_correction_max = 10;
		p_correction = 1e-4;

		p_BC_xb = 1.5;
		T_BC_xb = 1;

		p_BC_xe = 1;
		T_BC_xe = 1;

		dudx_0_BC_xb = false;
		durhodx_0_BC_xb = false;
		drhodt_durhodx_dvrhody_0_BC_xb = false;
		dudx_0_BC_xe = false;
		durhodx_0_BC_xe = false;
		drhodt_durhodx_dvrhody_0_BC_xe = false;
		u_uMin_Mout_0_BC_xe = false;

		dudy_0_BC_yb = false;
		dudy_0_BC_ye = false;


		dvdx_0_BC_xb = false;
		dvdx_0_BC_xe = false;

		dvdy_0_BC_yb = false;
		dvdy_0_BC_ye = false;


		dpdx_0_BC_xb = false;
		dpdx_0_BC_xe = false;
		dpdy_0_BC_yb = false;
		dpdy_0_BC_ye = false;

		dTdx_0_BC_xb = false;
		dTdx_0_BC_xe = false;
		dTdy_0_BC_yb = false;
		dTdy_0_BC_ye = false;

		ToSolveContinuityEquation = false;

		use_interpolation_for_u_BC_by_time = false;
		t_BC_b = 5;
		t_BC_e = 200;
		u_BC_xb_t_BC_b = 0.25;
		u_BC_xb_t_BC_e = 0.5;


		w_u = 1;
		w_v = 1;
		w_p = 1;
		w_Temper = 1;


		N_SubDomains_x = 1;
		N_SubDomains_y = 1;

		ToReadSolvedDataFromSeparateFileForEachSubDomain = true;
		ToWriteSolvedDataToFileForEntireComputationalDomain = true;
		ToWriteSolvedDataToSeparateFileForEachSubDomain = true;
		ToOnlyWriteDataFor_Nx_Ny_ForAllSubDomains_AndExit = false;
	}

	//This is just temporal to make file Nx_Ny_ForAllSubDomains.txt
	if(ToOnlyWriteDataFor_Nx_Ny_ForAllSubDomains_AndExit)
	{
		ToReadSolvedDataFromFile = false;
	}

	//Solve time for stop of a program
	time_solve_stop = time_solve_start + max_time_for_solve * 3600.0;
	solve_is_finished = false;

	rho_gas_nd = p_gas_nd / T_gas_nd;
	V_gas_nd = sqrt(u_gas_nd * u_gas_nd + v_gas_nd * v_gas_nd);


	//use diferent type of non domentional volues
	if(TypeOfNonDimensiolization == NonDimentionalGiven_Re_pch_rho_V2_2)
	{
		if(GivenReKn)
		{
			u_given_xb = Re * Kn * sqrt(15.0 * M_PI / 128.0);
			u_given_xb_in_it = u_given_xb;

			cT = (1.0 / Re) * (u_given_xb_in_it * p_BC_xb / T_BC_xb);
		}
		else
		{
			u_given_xb_in_it = u_given_xb;

			cT = (1.0 / Re) * (u_given_xb_in_it * p_BC_xb / T_BC_xb);
		}

		cPch = 0.5;
	}
	else if(TypeOfNonDimensiolization == NonDimentionalGiven_Re_pch_rho_R_T)
	{
		if(GivenReKn)
		{
			u_given_xb = Re * Kn * sqrt(15.0 * M_PI / 128.0);
			u_given_xb_in_it = u_given_xb;

			cT = (1.0 / Re) * (u_given_xb_in_it * p_BC_xb / T_BC_xb);
		}
		else
		{
			u_given_xb_in_it = u_given_xb;

			cT = (1.0 / Re) * (u_given_xb_in_it * p_BC_xb / T_BC_xb);
		}

		cPch = u_given_xb_in_it * u_given_xb_in_it / (gamma1 * Ma * Ma);
	}
	else if(TypeOfNonDimensiolization == NonDimentionalGiven_Kn_pch_rho_R_T)
	{
		//Solve parametars
		Ma = sqrt(6.0 / 5.0);
		Pr = 2.0 / 3.0;
		gamma1 = 5.0 / 3.0;
		c_mu = (5.0 / 16.0) * sqrt(2.0 * M_PI / gamma1);

		cT = Kn * 5.0 * sqrt(M_PI) / 16.0;
		cPch = 0.5;
	}

	if(SolveGravityField) c_Fr_v = 2.0 / (gamma1 * Fr * Ma * Ma);
	else c_Fr_v = 0;

	//here parameters for solving equation for Temper are only for case 3 - TypeOfNonDimensiolization == NonDimentionalGiven_Kn_pch_rho_R_T
	cT1 = gamma1 * c_mu * Kn / (Ma * Pr);
	cT2 = gamma1 * (gamma1 - 1.0) * c_mu * Kn * Ma;

	//cT1 = c_mu * Kn / (Ma * Pr);
	//cT2 = (gamma1 - 1.0) * c_mu * Kn * Ma;


	//Solve hx and hy
	if(kind_of_mesh == Par_mesh
	|| kind_of_mesh == Par_mesh_flat)
	{
		xmid = xfmin + 0.5 * (xbmin - xfmin);

		axf2 = 0;
		axf2 = (hxmax - hxmin) / ((xfmax - xfmin) * (xfmax - xfmin));
		axf1 = -2.0 * axf2 * xfmin;
		axf0 = hxmin + axf2 * xfmin * xfmin;

		axb2 = 0;
		axb2 = (hxmax - hxmin) / ((xbmax - xbmin) * (xbmax - xbmin));
		axb1 = -2.0 * axb2 * xbmin;
		axb0 = hxmin + axb2 * xbmin * xbmin;

		hxmid = hx_all_par(xmid);


		ymid = ybmin + 0.5 * (ytmin - ybmin);

		ayt2 = 0;
		ayt2 = (hymax - hymin) / ((ytmax - ytmin) * (ytmax - ytmin));
		ayt1 = -2.0 * ayt2 * ytmin;
		ayt0 = hymin + ayt2 * ytmin * ytmin;

		ayb2 = 0;
		ayb2 = (hymax - hymin) / ((ybmax - ybmin) * (ybmax - ybmin));
		ayb1 = -2.0 * ayb2 * ybmin;
		ayb0 = hymin + ayb2 * ybmin * ybmin;

		hymid = hy_all_par(ymid);
	}
	else if (kind_of_mesh == Polynom3_flat_mesh)
	{
		double axf_tmp;
		axf_tmp = 1.0 / (3.0 * xfmax * pow(xfmin, 2) - pow(xfmin, 3) - 3.0 * pow(xfmax, 2) * xfmin + pow(xfmax, 3));

		axf3 = 2.0 * (hxmin - hxmax) * axf_tmp;
		axf2 = -3.0 * (xfmax + xfmin) * (hxmin - hxmax) * axf_tmp;
		axf1 = 6.0 * xfmin * xfmax * (hxmin - hxmax) * axf_tmp;
		axf0 = (-3.0 * hxmin * pow(xfmax, 2) * xfmin + hxmin * pow(xfmax, 3) + 3.0 * hxmax * pow(xfmin, 2) * xfmax - hxmax * pow(xfmin, 3)) * axf_tmp;


		double axb_tmp;
		axb_tmp = 1.0 / (3.0 * xbmin * pow(xbmax, 2) - pow(xbmax, 3) - 3.0 * pow(xbmin, 2) * xbmax + pow(xbmin, 3));

		axb3 = 2.0 * (hxmax - hxmin) * axb_tmp;
		axb2 = -3.0 * (xbmin + xbmax) * (hxmax - hxmin) * axb_tmp;
		axb1 = 6.0 * xbmax * xbmin * (hxmax - hxmin) * axb_tmp;
		axb0 = (-3.0 * hxmax * pow(xbmin, 2) * xbmax + hxmax * pow(xbmin, 3) + 3.0 * hxmin * pow(xbmax, 2) * xbmin - hxmin * pow(xbmax, 3)) * axb_tmp;



		double ayb_tmp;
		ayb_tmp = 1.0 / (3.0 * ybmax * pow(ybmin, 2) - pow(ybmin, 3) - 3.0 * pow(ybmax, 2) * ybmin + pow(ybmax, 3));

		ayb3 = 2.0 * (hxmin - hxmax) * ayb_tmp;
		ayb2 = -3.0 * (ybmax + ybmin) * (hxmin - hxmax) * ayb_tmp;
		ayb1 = 6.0 * ybmin * ybmax * (hxmin - hxmax) * ayb_tmp;
		ayb0 = (-3.0 * hxmin * pow(ybmax, 2) * ybmin + hxmin * pow(ybmax, 3) + 3.0 * hxmax * pow(ybmin, 2) * ybmax - hxmax * pow(ybmin, 3)) * ayb_tmp;


		double ayt_tmp;
		ayt_tmp = 1.0 / (3.0 * ytmin * pow(ytmax, 2) - pow(ytmax, 3) - 3.0 * pow(ytmin, 2) * ytmax + pow(ytmin, 3));

		ayt3 = 2.0 * (hxmax - hxmin) * ayt_tmp;
		ayt2 = -3.0 * (ytmin + ytmax) * (hxmax - hxmin) * ayt_tmp;
		ayt1 = 6.0 * ytmax * ytmin * (hxmax - hxmin) * ayt_tmp;
		ayt0 = (-3.0 * hxmax * pow(ytmin, 2) * ytmax + hxmax * pow(ytmin, 3) + 3.0 * hxmin * pow(ytmax, 2) * ytmin - hxmin * pow(ytmax, 3)) * ayt_tmp;

	}



	//Define Nx and Ny for all computational domain.
	if(!ToReadSolvedDataFromFile)
	{
		define_Nx_Ny();
		//WriteEnteredDataToFile();
	}

	Na = Nx * Ny;


	//Make mesh for all compuattional dimain.
	//Further if are used subdomains the mesh of subdomain will be
	//reorganized. This is make for simpicity  of define x_b, x_e, y_b and y_e.
	if(!is_defined_array_for_mesh)
	{
		x_f = new double [Nx + 1];
		y_f = new double [Ny + 1];

		hx = new double [Nx];
		hy = new double [Ny];

		x_v = new double [Nx];
		y_v = new double [Ny];

		is_defined_array_for_mesh = true;
	}


	//we want x_v[0] = x_b;
	x_f[0] = x_b - 0.5 * hxmax;
	x_f[1] = x_b + 0.5 * hxmax;
	for(i = 1; i < Nx; i++)
	{
		double x_f_tmp;
		x_f_tmp = x_f[i] + 0.5 * hx_all(x_f[i]);

		x_f[i + 1] = x_f[i] + hx_all(x_f_tmp);
	}

	//we want y_v[0] = y_b;
	y_f[0] = y_b - 0.5 * hymax;
	y_f[1] = y_b + 0.5 * hymax;
	for(j = 1; j < Ny; j++)
	{
		double y_f_tmp;
		y_f_tmp = y_f[j] + 0.5 * hy_all(y_f[j]);

		y_f[j + 1] = y_f[j] + hy_all(y_f_tmp);
	}


	for(i = 0; i < Nx; i++)
	{
		hx[i] = x_f[i + 1] - x_f[i];
		x_v[i] = x_f[i] + 0.5 * hx[i];
	}

	for(j = 0; j < Ny; j++)
	{
		hy[j] = y_f[j + 1] - y_f[j];
		y_v[j] = y_f[j] + 0.5 * hy[j];
	}


	//start iteration to solve grid
	if(N_I_T > 0)
	{
		Iter_T = 0;
		do
		{
			for(i = 1; i < Nx; i++)
			{
				x_f[i + 1] = x_f[i] + hx_all(x_v[i]);
			}

			for(j = 1; j < Ny; j++)
			{
				y_f[j + 1] = y_f[j] + hy_all(y_v[j]);
			}

			for(i = 0; i < Nx; i++)
			{
				hx[i] = x_f[i + 1] - x_f[i];
				x_v[i] = x_f[i] + 0.5 * hx[i];
			}

			for(j = 0; j < Ny; j++)
			{
				hy[j] = y_f[j + 1] - y_f[j];
				y_v[j] = y_f[j] + 0.5 * hy[j];
			}


			Iter_T++;
		}while(Iter_T < N_I_T);

	}


#ifndef COMPILE_MPI
	//write mesh to files
	WriteMassiveToFile_1D("x_f.txt", x_f, Nx+1);
	WriteMassiveToFile_1D("hx.txt", hx, Nx);
	WriteMassiveToFile_1D("x_v.txt", x_v, Nx);

	WriteMassiveToFile_1D("y_f.txt", y_f, Ny+1);
	WriteMassiveToFile_1D("hy.txt", hy, Ny);
	WriteMassiveToFile_1D("y_v.txt", y_v, Ny);


	if(!ToReadSolvedDataFromFile)
	{
		WriteEnteredDataToFile();
	}
#endif

	//Before making mesh and after define Nx and Ny finaly,
	//have to define data of subdomains if the calculations are with MPI.
#ifdef COMPILE_MPI

	//Define SubDomains and mesh of this domain
	DefineSubDomains_MPI();


	//This is option to write only number of subdomains if there is need to make change.
	//The program have to be made to start from begining just for this this change.
	if(ToOnlyWriteDataFor_Nx_Ny_ForAllSubDomains_AndExit)
	{
		//Terminate MPI
		MPI_Finalize();

		return 0;
	}

	//Check condition to continue calculation
	//The condition of continue calculation is to have
	//equality of N_SubDomains_a == N_SubDomains_x * N_SubDomains_y;
	if((N_SubDomains_a == (N_SubDomains_x * N_SubDomains_y))
	&& ((N_SubDomains_x == 1) || (N_SubDomains_y == 1)))
	{
		if(IJ == 0)
		{
			//The condition of continue calculation is to have
			//equality of N_SubDomains_a == N_SubDomains_x * N_SubDomains_y;
			std::cout << "Here is equality of N_SubDomains_a == N_SubDomains_x * N_SubDomains_y" << "\n";
			std::cout << N_SubDomains_a << " == " << N_SubDomains_x * N_SubDomains_y << "\n";
			std::cout << "The program will continue." << "\n";
		}
	}
	else if(N_SubDomains_a != (N_SubDomains_x * N_SubDomains_y))
	{
		if(IJ == 0)
		{
			//Interrupt program:
			std::cout << "The program will be interrupted because:" << "\n";
			std::cout << "N_SubDomains_a != N_SubDomains_x * N_SubDomains_y;" << "\n";

			std::cout << "Run program again and use CPU number equal to: N_SubDomains_x * N_SubDomains_y = "
				<< N_SubDomains_x * N_SubDomains_y << "\n"
				<< "Or change SubDomains." << "\n";
		}


		//Terminate MPI
		MPI_Finalize();

		return 0;
	}
	else
	{
		//Here mean that (1 < N_SubDomains_x) || (1 < N_SubDomains_y)
		//The program is made to separate SubDomains on OX only or OY only!!!
		if(IJ == 0)
		{
			//Interrupt program:
			std::cout << "The program will be interrupted because:" << "\n"
				<< "(N_SubDomains_x == 1) || (N_SubDomains_y == 1);" << "\n"
				<< "The program can separate computational domain on subdomains on only OX or only on OY!" << "\n"
				<< "Can not be made correct domain decomposition if N_SubDomains_x > 1 and N_SubDomains_y > 1 at this time." << "\n"
				<< "This is not done, becouse will increase comunications and code will be more compicated." << "\n"
				<< "This kind of program satisfy good problems, which are calculated and resourses, where can be run." << "\n"
				<< "Run program again and use for number of subdoamins on OX 1 or on OY 1." << "\n";
		}


		//Terminate MPI
		MPI_Finalize();

		return 0;
	}

	//The rank 0 MUST write to file when no other process is read the file.
	//Therefore all processe have to reache this point until rank 0 write data
	//to file EnteredData.txt.
	MPI_Barrier(MPI_COMM_WORLD);


	if((!ToReadSolvedDataFromFile)
	&& (IJ == 0))
	{
		WriteEnteredDataToFile();
		WriteDataFor_Nx_Ny_ForAllSubDomains();

		//Delete arrays Nx_all_SubDomains and Ny_all_SubDomains;
		if(is_defined_arrays_for_Nx_Ny_for_all_SubDomains)
		{
			Nx_all_SubDomains = NULL; delete [] Nx_all_SubDomains;
			Ny_all_SubDomains = NULL; delete [] Ny_all_SubDomains;
		}

		is_defined_arrays_for_Nx_Ny_for_all_SubDomains = false;
	}

	MPI_Barrier(MPI_COMM_WORLD);

#endif //End of COMPILE_MPI


	//Tests for mesh are OK, for non-periodic BC.
	//For PeriodicBC have some problems, which will be solved later.


	//define vol_inf_fb, u, v and p
	define_vol_inf();

	if(ToSolveContinuityEquation)
	{
		ContinuityEquation_stationary = new double[Na];
		ContinuityEquation = new double[Na];
	}

	//sqrt_T_ff_V = new double[Na];
	Gamma_yf = new double[Na];
	Gamma_xf = new double[Na];

	Fx = new double[Na];
	Fy = new double[Na];

	storage_u_0 = new double[Na];
	storage_u_1 = new double[Na];
	u = storage_u_0;
	u_pr = storage_u_1;
	u_pseudo = new double[Na];
	D_ux = new double[Na];
	D_uy = new double[Na];

	storage_v_0 = new double[Na];
	storage_v_1 = new double[Na];
	v = storage_v_0;
	v_pr = storage_v_1;
	v_pseudo = new double[Na];
	D_vx = new double[Na];
	D_vy = new double[Na];

	storage_p_0 = new double[Na];
	storage_p_1 = new double[Na];
	p = storage_p_0;
	p_pr = storage_p_1;

#ifdef Calculate_p_and_Temper_using_Jakoby_method
	storage_p_2 = new double[Na];
	p_old = storage_p_2;
#endif //End of Calculate_p_and_Temper_using_Jakoby_method

	d_p_c_12 = new double[Na];
	d_p_c_34 = new double[Na];

	apx = new double[Na];
	apy = new double[Na];
	bpx = new double[Na];
	bpy = new double[Na];


	storage_rho_0 = new double[Na];
	storage_rho_1 = new double[Na];
	rho = storage_rho_0;
	rho_pr = storage_rho_1;

	rho_u = new double[Na];
	rho_v = new double[Na];


	storage_Temper_0 = new double[Na];
	storage_Temper_1 = new double[Na];
	Temper = storage_Temper_0;
	Temper_pr = storage_Temper_1;

	Temper_of_fluid_on_the_wall_u = new double[Na];
	Temper_of_fluid_on_the_wall_v = new double[Na];

	for(counter = 0; counter < Na; counter++)
	{
		//This is 0, to make easy to find possible errors in program.
		Temper_of_fluid_on_the_wall_u[counter] = 0;
		Temper_of_fluid_on_the_wall_v[counter] = 0;
	}

#ifdef Calculate_p_and_Temper_using_Jakoby_method
	storage_Temper_2 = new double[Na];
	Temper_old = storage_Temper_2;
#endif //End of Calculate_p_and_Temper_using_Jakoby_method

	sqrt_T = new double[Na];
	DTx = new double[Na];
	DTy = new double[Na];
	//FTx = new double[Na];
	//FTy = new double[Na];
	ScT1 = new double[Na];
	//ScT2 = new double[Na];
	//bT0 = new double[Na];
	//sqrt_T_T_bx = new double[Na];
	//sqrt_T_T_by = new double[Na];


	max_residual_in_p = new double[N_I];
	max_residual_in_T = new double[N_I];

	max_residual_in_p_Iter_0 = new double[N_I];
	max_residual_in_T_Iter_0 = new double[N_I];

	max_Iter_p_c = new unsigned int[N_I];


	//Initialize the matrixes
	for(counter = 0; counter < Na; counter++)
	{
		u[counter] = u_gas_nd;
		v[counter] = v_gas_nd;
		p[counter] = p_gas_nd;
		rho[counter] = rho_gas_nd;
		Temper[counter] = T_gas_nd;
	}


	define_area_for_solving();
	define_rule_vectors();

#ifdef DEBUG_PROGRAM
	//Write rule vectors to files
	WriteMassiveToFile_1D("rule_vector_ij_nch_u.txt", rule_vector_ij_nch_u, Na_nch_u, 15);
	WriteMassiveToFile_1D("rule_vector_ij_bc_u.txt", rule_vector_ij_bc_u, Na_bc_u, 15);

	WriteMassiveToFile_1D("rule_vector_ij_nch_v.txt", rule_vector_ij_nch_v, Na_nch_v, 15);
	WriteMassiveToFile_1D("rule_vector_ij_bc_v.txt", rule_vector_ij_bc_v, Na_bc_v, 15);

	WriteMassiveToFile_1D("rule_vector_ij_nch_others.txt", rule_vector_ij_nch_others, Na_nch_others, 15);
	WriteMassiveToFile_1D("rule_vector_ij_bc_others.txt", rule_vector_ij_bc_others, Na_bc_others, 15);
#endif //End of #ifdef DEBUG_PROGRAM


	if(ToReadSolvedDataFromFile)
	{
		ReadSolvedDataFromFile();
	}
	else
	{
		//If the program start from begining calclulation
		//Apply conditions fron areas
		ApplyBodiesConditionsToVariables();
		//ApplyPressureConditionsToPressure();
		//ApplyVelocityConditionsToVelocities();


		//Apply BC for pressure
		rho_BC_xb = p_BC_xb / T_BC_xb;
		rho_BC_xe = p_BC_xe / T_BC_xe;

		if(Pressure_BC)
		{
			//Inflow BC -> x = x_b
			for(j = 0; j < Ny; j++)
			{
				for(i = 0; i < Nx_others[0]; i++)
				{
					node = i + j * Nx;

					if(vol_inf_fb_bool[node]
#ifdef COMPILE_MPI
					//This BC is used only for the inflow of the channel.
					&& (!SubDomain[IJ].is_SubDomain_in_direction_I_1)
#endif //End of COMPILE_MPI
					)
					{
						p[node] = p_BC_xb;
						rho[node] = rho_BC_xb;
						Temper[node] = T_BC_xb;
					}

				}
			}

			//Outflow BC -> x = x_e
			for(j = 0; j < Ny; j++)
			{
				for(i = Nx_others[1]; i < Nx; i++)
				{
					node = i + j * Nx;

					if(vol_inf_fb_bool[node]
#ifdef COMPILE_MPI
					//This BC is used only for the inflow of the channel.
					&& (!SubDomain[IJ].is_SubDomain_in_direction_I1)
#endif //End of COMPILE_MPI
					)
					{
						p[node] = p_BC_xe;
						rho[node] = rho_BC_xe;
						Temper[node] = T_BC_xe;
					}

				}
			}


			if(ToStartFromExactSolutionForChannelFlowWithPressureBC)
			{

				//layer on OX in wich I check begin index of channel
				i = 1;

				for(j = 1; j < Ny; j++)
				{
					node = i + j * Nx;

					if((!vol_inf_fb_bool[node - Nx]) && vol_inf_fb_bool[node])
						j_channel_b = j - 1;

					if(vol_inf_fb_bool[node - Nx] && (!vol_inf_fb_bool[node]))
						j_channel_e = j;

				}

				h_channel = y_v[j_channel_e] - y_v[j_channel_b];

				L_channel = x_v[Nx_others[1]] - x_v[Nx_others[0] - 1];

				dpdx_PressureBC = (p_BC_xe - p_BC_xb) / L_channel;

				body_SlipBC_for_exact_solution = 0;

				double sigma_tmp, SlipBC_tmp, Px_tmp;

				SlipBC_tmp = OforS[gb].M_polyhedrons[body_SlipBC_for_exact_solution].is_VelocitySlipBC_body;

				sigma_tmp = SlipBC_tmp * OforS[gb].M_polyhedrons[body_SlipBC_for_exact_solution].F_VelocitySlip;

				Px_tmp = p_BC_xb / p_BC_xe;

				double * x_v_A, * y_v_A;
				x_v_A = new double [Nx];
				y_v_A = new double [Ny];

				for(i = 0; i < Nx; i++)
				{
					x_v_A[i] = (x_v[i] - x_v[0]) / L_channel;
				}

				for(j = 0; j < Ny; j++)
				{
					if(j <= j_channel_b)
					{
						y_v_A[j] = -0.5 * h_channel;
					}
					else if(j_channel_e <= j)
					{
						y_v_A[j] = 0.5 * h_channel;
					}
					else
					{
						y_v_A[j] = y_v[j + j_channel_b] - (y_v[j_channel_b] + 0.5 *  h_channel);
					}

				}


				for(j = 0; j < Ny; j++)
				{
					for(i = Nx_others[0]; i < Nx_others[1]; i++)
					{
						node = i + j * Nx;

						if(vol_inf_fb_bool[node])
						{
							p[node] = p_BC_xe * (-6.0 * sigma_tmp * Kn +
								sqrt(pow((6.0 * sigma_tmp * Kn), 2) + (1.0 + 12.0 * sigma_tmp * Kn) * x_v_A[i]
								+ (Px_tmp + 12.0 * sigma_tmp * Kn) * Px_tmp * (1.0 - x_v_A[i])));

						}
					}
				}


				double ND_parameter;

				//use diferent type of non domentional volues
				if(TypeOfNonDimensiolization == NonDimentionalGiven_Re_pch_rho_V2_2)
				{
					ND_parameter = Re * h_channel * h_channel / (4.0 * u_given_xb);
				}
				else if(TypeOfNonDimensiolization == NonDimentionalGiven_Re_pch_rho_R_T)
				{
					ND_parameter = Re * h_channel * h_channel / (4.0 * u_given_xb);
				}
				else if(TypeOfNonDimensiolization == NonDimentionalGiven_Kn_pch_rho_R_T)
				{
					ND_parameter = h_channel / (5.0 * sqrt(M_PI * T_BC_xb) * Kn * L_channel);
				}



				for(j = 0; j < Ny; j++)
				{
					for(i = 0; i < Nx; i++)
					{
						node = i + j * Nx;

						if(vol_inf_fb_bool[node])
						{

							u[node] =

								//this is Slip BC (sigma == 0 - NO Slip BC)
								-ND_parameter *
								//dpdx
								((p_BC_xe * (1.0 - Px_tmp * Px_tmp + 12.0 * sigma_tmp * Kn * (1.0 - Px_tmp)))
								/ (2.0 * sqrt(pow((6.0 * sigma_tmp * Kn), 2) + (1.0 + 12.0 * sigma_tmp * Kn) * x_v_A[i]
								+ (Px_tmp + 12.0 * sigma_tmp * Kn) * Px_tmp * (1.0 - x_v_A[i])))) *

								(1.0 - 4.0 * y_v_A[j] * y_v_A[j] / (h_channel * h_channel)
								+ 4.0 * sigma_tmp * Kn / p[node]);

						}

					}
				}

			}


		}


		if(Given_velocity_on_xb && ToStartFromExactSolutionForChannelFlowWithPressureBC)
		{
			//Inflow BC -> x = x_b
			for(j = 0; j < Ny; j++)
			{
				for(i = 0; i < Nx_others[0] + 1; i++)
				{
					node = i + j * Nx;

					if(vol_inf_fb_bool[node]
#ifdef COMPILE_MPI
					//This BC is used only for the inflow of the channel.
					&& (!SubDomain[IJ].is_SubDomain_in_direction_I_1)
#endif //End of COMPILE_MPI
					)
					{
						p[node] = p_BC_xb;
					}
				}
			}


			if(ToStartFromExactSolutionForChannelFlowWithPressureBC)
			{

				//layer on OX in wich I check begin index of channel
				i = 1;

				for(j = 1; j < Ny; j++)
				{
					node = i + j * Nx;

					if((!vol_inf_fb_bool[node - Nx]) && vol_inf_fb_bool[node])
						j_channel_b = j - 1;

					if(vol_inf_fb_bool[node - Nx] && (!vol_inf_fb_bool[node]))
						j_channel_e = j;

				}

				h_channel = y_f[j_channel_e + 1 * (0 < j)] - y_f[j_channel_b];

				L_channel = x_v[Nx_others[1]] - x_v[Nx_others[0] - 1];

				double y_center_of_channel;
				y_center_of_channel = y_v[j_channel_b] + 0.5 * h_channel;

				if(correct_p_BC_xb)
				{
					//This is pressure for u[center_of_channel] = u_given_xb, when we have given profile
					p_BC_xb = p_BC_xe + 4.0 * L_channel * u_given_xb * u_given_xb
						/ (Re * h_channel * h_channel
							* (y_center_of_channel - y_v[j_channel_b] / h_channel
								- (y_center_of_channel - y_v[j_channel_b] / h_channel) * (y_center_of_channel - y_v[j_channel_b] / h_channel)));
				}
				else
				{
					//This is pressure for u[center_of_channel] = u_given_xb, when we have given profile
					p_BC_xe = p_BC_xb - 4.0 * L_channel * u_given_xb * u_given_xb
						/ (Re * h_channel * h_channel
							* (y_center_of_channel - y_v[j_channel_b] / h_channel
								- (y_center_of_channel - y_v[j_channel_b] / h_channel) * (y_center_of_channel - y_v[j_channel_b] / h_channel)));
				}

				dpdx_PressureBC = (p_BC_xe - p_BC_xb) / L_channel;

				body_SlipBC_for_exact_solution = 0;

				double ND_parameter;

				//use diferent type of non domentional volues
				if(TypeOfNonDimensiolization == NonDimentionalGiven_Re_pch_rho_V2_2)
				{
					ND_parameter = Re * h_channel * h_channel / (4.0 * u_given_xb);
				}
				else if(TypeOfNonDimensiolization == NonDimentionalGiven_Re_pch_rho_R_T)
				{
					ND_parameter = Re * h_channel * h_channel / (4.0 * u_given_xb);
				}
				else if(TypeOfNonDimensiolization == NonDimentionalGiven_Kn_pch_rho_R_T)
				{
					ND_parameter = 4.0 * h_channel * h_channel / (5.0 * sqrt(M_PI * T_BC_xb) * Kn);
				}



				for(j = 0; j < Ny; j++)
				{
					//Given_velocity_on_xb - to make only first layers by this method
					for(i = 0; i < Nx_u[0]; i++)
					{
						node = i + j * Nx;

						if(vol_inf_fb_bool[node])
						{
							//this is NO Slip BC
							u[node] = - ND_parameter * dpdx_PressureBC *
								((y_v[j] - y_v[j_channel_b]) /  h_channel
									- ((y_v[j] - y_v[j_channel_b]) /  h_channel) * ((y_v[j] - y_v[j_channel_b]) /  h_channel)

									+ Kn * OforS[gb].M_polyhedrons[body_SlipBC_for_exact_solution].F_VelocitySlip *
									OforS[gb].M_polyhedrons[body_SlipBC_for_exact_solution].is_VelocitySlipBC_body / h_channel//sigmea from polyhedron 0
									);

						}

					}
				}

			}


		}

		if(Given_velocity_on_xe && ToStartFromExactSolutionForChannelFlowWithPressureBC)
		{
			//Outflow BC -> x = x_e
			for(j = 0; j < Ny; j++)
			{
				for(i = Nx_others[1]; i < Nx; i++)
				{
					node = i + j * Nx;

					if(vol_inf_fb_bool[node]
#ifdef COMPILE_MPI
					//This BC is used only for the inflow of the channel.
					&& (!SubDomain[IJ].is_SubDomain_in_direction_I1)
#endif //End of COMPILE_MPI
					)
					{
						p[node] = p_BC_xe;
					}

				}
			}

		}


		for(counter = 0; counter < Na; counter++)
		{
			p_pr[counter] = p[counter];
			u_pr[counter] = u[counter];
			v_pr[counter] = v[counter];

			Temper_pr[counter] = Temper[counter];
			rho_pr[counter] = rho[counter];
		}
		Solve_rho_in_middle_points();


		ApplyBodiesConditionsToVariables();
		BoundaryConditions_Velocities();
		//BoundaryConditions_Velocities_SlipBC_u();
		//BoundaryConditions_Velocities_SlipBC_v();

		for(counter = 0; counter < Na; counter++)
		{
			p_pr[counter] = p[counter];
			u_pr[counter] = u[counter];
			v_pr[counter] = v[counter];
			Temper_pr[counter] = Temper[counter];
			rho_pr[counter] = rho[counter];
		}
		Solve_rho_in_middle_points();

		it_begin = 0;
		it = it_begin;

		WriteSolvedDataToFile();
	}


	//Initialize the matrixes
	for(counter = 0; counter < Na; counter++)
	{
		d_p_c_12[counter] = 0;
		d_p_c_34[counter] = 0;

		apx[counter] = 0;
		apy[counter] = 0;
		bpx[counter] = 0;
		bpy[counter] = 0;


		if(ToSolveContinuityEquation)
		{
			ContinuityEquation_stationary[counter] = 0;
			ContinuityEquation[counter] = 0;
		}

		Gamma_yf[counter] = 0;
		Gamma_xf[counter] = 0;

		Fx[counter] = 0;
		Fy[counter] = 0;


		u_pseudo[counter] = 0;
		v_pseudo[counter] = 0;
	}



	bool Given_velocity_on_xb_bool_tmp, Given_velocity_on_xe_bool_tmp;

	Nt_save_solved_data_counter = 0;
	Nt_save_DragCoefficient_counter = 0;
	Nt_save_temp_data_counter = 0;

	ToCheckForFileIsNow_Tmp = true;
	ToCheckForFileIsNow_Iterations_in_loops = true;
	ToCheckForFileIsNow_CD = true;
	ToCheckForFileIsNow_p_u_x_b = true;



	//inicialize matrixes before start solving
	for(counter = 0; counter < Na; counter++)
	{
		D_ux[counter] = 0;
		D_uy[counter] = 0;

		D_vx[counter] = 0;
		D_vy[counter] = 0;


		DTx[counter] = 0;
		DTy[counter] = 0;
		ScT1[counter] = 0;

		rho_u[counter] = 0;
		rho_v[counter] = 0;
	}


	Iter = 0;
	Iter_pr1 = Iter;
	Iter_pr2 = Iter_pr1;

	max_residual_in_u = 0.0;
	max_residual_in_v = 0.0;

	max_residual_in_u_Iter_0 = 0.0;
	max_residual_in_v_Iter_0 = 0.0;


	//Apply Pressure_BC
	if(Pressure_BC && (!ToReadSolvedDataFromFile))
	{
		for(j = 0; j < Ny; j++)
		{
            //Apply BC on the x_b
			for(i = 0; i < Nx_others[0]; i++)
			{
				node = i + j * Nx;
				if(vol_inf_fb_bool[node]
#ifdef COMPILE_MPI
				//This BC is used only for the inflow of the channel.
				&& (!SubDomain[IJ].is_SubDomain_in_direction_I_1)
#endif //End of COMPILE_MPI
				)
				{
					p[node] = p_BC_xb;
					Temper[node] = T_BC_xb;
				}
			}

            //Apply BC on the x_e
			for(i = Nx_others[1]; i < Nx; i++)
			{
				node = i + j * Nx;
				if(vol_inf_fb_bool[node]
#ifdef COMPILE_MPI
				//This BC is used only for the inflow of the channel.
				&& (!SubDomain[IJ].is_SubDomain_in_direction_I1)
#endif //End of COMPILE_MPI
				)
				{
					p[node] = p_BC_xe;
					Temper[node] = T_BC_xe;
				}
			}

		}
	}


	if(ToContinueFromInterpolation)
	{
		//The calculation of velocity slip and temperature jump boundary
		//conditions are included directly in equations and the values on boundaries.
		ApplyBodiesConditionsToVariables();
		//BoundaryConditions_Velocities();
		BoundaryConditions_Velocities_SlipBC_u();
		BoundaryConditions_Velocities_SlipBC_v();
		BoundaryConditions_Temperature_TemperatureJumpBC();
		BoundaryConditions_Temperature_TemperatureJumpBC_WriteToMassiveToWriteDataToFile();
	}


	for(counter = 0; counter < Na; counter++)
	{
		sqrt_T[counter] = sqrt(Temper[counter]);

		if(vol_inf_fb_bool[counter]) rho[counter] = p[counter] / Temper[counter];
		else rho[counter] = 0;
	}
	Solve_rho_in_middle_points();


#ifdef COMPILE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif //End of COMPILE_MPI

#ifdef Calculate_p_and_Temper_using_Jakoby_method
	memcpy(p_old, p, Na * sizeof (double));
	memcpy(Temper_old, Temper, Na * sizeof (double));
#endif //End of Calculate_p_and_Temper_using_Jakoby_method


	//solve by SIMPLE compressible with h variable
	it = it_begin;
	it_counter = 0;
	countinue_time_step = true;

	for(counter = 0; counter < Na; counter++)
	{
		//Initiate values from previous time step
		p_pr[counter] = p[counter];

		u_pr[counter] = u[counter];
		v_pr[counter] = v[counter];

		Temper_pr[counter] = Temper[counter];
		sqrt_T[counter] = sqrt(Temper[counter]);

		rho_pr[counter] = rho[counter];
	}

#ifdef DEBUG_PROGRAM
	WriteSolvedDataToFile();
#endif


	//Calculate before start the loop 2 - all of the following calculation are placed
	//at the end of loop 2 to hide latency for paralel program using MPI.

	//use different type of non-dimentional volues
	if(TypeOfNonDimensiolization == NonDimentionalGiven_Re_pch_rho_V2_2)
	{
		cT = (1.0 / Re) * (u_given_xb_in_it * p_BC_xb / T_BC_xb);
		cPch = 0.5;
	}
	else if(TypeOfNonDimensiolization == NonDimentionalGiven_Re_pch_rho_R_T)
	{
		cT = (1.0 / Re) * (u_given_xb_in_it * p_BC_xb / T_BC_xb);
		cPch = u_given_xb_in_it * u_given_xb_in_it / (gamma1 * Ma * Ma);
	}
	else if(TypeOfNonDimensiolization == NonDimentionalGiven_Kn_pch_rho_R_T)
	{
		cT = Kn * 5.0 * sqrt(M_PI) / 16.0;
		cPch = 0.5;
	}

	//Calculate rho_u and rho_v after calculation of rho,
	//which is used in calculation of fluxes (Fx and Fy)
	//and procedure BoundaryConditions_Velocities()
	//Solve rho in middle points (rho_u and rho_v) if this is first
	//iteration after starting calculation
	Solve_rho_in_middle_points();

#ifdef DEBUG_PROGRAM
#ifdef COMPILE_MPI
	string rho_u_filename = "rho_u." + toString(IJ) + ".bin";
	string rho_v_filename = "rho_v." + toString(IJ) + ".bin";

	WriteMassiveToBinaryFile(rho_u_filename.c_str(), rho_u, Na);
	WriteMassiveToBinaryFile(rho_v_filename.c_str(), rho_v, Na);
#else
	WriteMassiveToBinaryFile("rho_u.bin", rho_u, Na);
	WriteMassiveToBinaryFile("rho_v.bin", rho_v, Na);
#endif //End of COMPILE_MPI
#endif //End of DEBUG_PROGRAM

	//rho_u have to be calculated
	BoundaryConditions_Velocities_SlipBC_u();
#ifdef DEBUG_PROGRAM
#ifdef COMPILE_MPI
	string u_filename = "u." + toString(IJ) + ".bin";
	WriteMassiveToBinaryFile(u_filename.c_str(), u, Na);

	MPI_Barrier(MPI_COMM_WORLD);
	Write_massive_entire_computational_domain_MPI("u", u, double_type);
	MPI_Barrier(MPI_COMM_WORLD);
#else
	WriteMassiveToBinaryFile("u.bin", u, Na);
#endif //End of COMPILE_MPI
#endif //End of DEBUG_PROGRAM

	//rho_v have to be calculated
	BoundaryConditions_Velocities_SlipBC_v();
#ifdef DEBUG_PROGRAM
#ifdef COMPILE_MPI
	string v_filename = "v." + toString(IJ) + ".bin";
	WriteMassiveToBinaryFile(v_filename.c_str(), v, Na);

	MPI_Barrier(MPI_COMM_WORLD);
	Write_massive_entire_computational_domain_MPI("v", v, double_type);
	MPI_Barrier(MPI_COMM_WORLD);
#else
	WriteMassiveToBinaryFile("v.bin", v, Na);
#endif //End of COMPILE_MPI
#endif //End of DEBUG_PROGRAM

	//Boundary Conditions
	BoundaryConditions_Velocities();
#ifdef DEBUG_PROGRAM
#ifdef COMPILE_MPI
	WriteMassiveToBinaryFile(u_filename.c_str(), u, Na);
	WriteMassiveToBinaryFile(v_filename.c_str(), v, Na);

	MPI_Barrier(MPI_COMM_WORLD);
	Write_massive_entire_computational_domain_MPI("u", u, double_type);
	MPI_Barrier(MPI_COMM_WORLD);
	Write_massive_entire_computational_domain_MPI("v", v, double_type);
	MPI_Barrier(MPI_COMM_WORLD);
#else
	WriteMassiveToBinaryFile("u.bin", u, Na);
	WriteMassiveToBinaryFile("v.bin", v, Na);
#endif //End of COMPILE_MPI
#endif //End of DEBUG_PROGRAM

	//apply BC for correct pressure when we have given
	//velocity in x_b on pressure driven flow
	BoundaryConditions_Pressure_begin_iteration_for_it();

	//Apply temperature jump BC, after calculation of Temper and rho.
	BoundaryConditions_Temperature_TemperatureJumpBC();
#ifdef DEBUG_PROGRAM
#ifdef COMPILE_MPI
	string Temper_of_fluid_on_the_wall_u_filename = "Temper_of_fluid_on_the_wall_u." + toString(IJ) + ".bin";
	string Temper_of_fluid_on_the_wall_v_filename = "Temper_of_fluid_on_the_wall_v." + toString(IJ) + ".bin";

	WriteMassiveToBinaryFile(Temper_of_fluid_on_the_wall_u_filename.c_str(), Temper_of_fluid_on_the_wall_u, Na);
	WriteMassiveToBinaryFile(Temper_of_fluid_on_the_wall_v_filename.c_str(), Temper_of_fluid_on_the_wall_v, Na);
#else
	WriteMassiveToBinaryFile("Temper_of_fluid_on_the_wall_u.bin", Temper_of_fluid_on_the_wall_u, Na);
	WriteMassiveToBinaryFile("Temper_of_fluid_on_the_wall_v.bin", Temper_of_fluid_on_the_wall_v, Na);
#endif //End of COMPILE_MPI
#endif //End of DEBUG_PROGRAM

	BoundaryConditions_Temperature();
	BoundaryConditions_Pressure();

	//Solve fluxes (Fx and Fy) and diffusion terms if this is first iteration after
	//starting calculation
	SolveFluxes();
	SolveDiffisionTerms();

#ifdef DEBUG_PROGRAM
#ifdef COMPILE_MPI
	string Fx_filename = "Fx." + toString(IJ) + ".bin";
	string Fy_filename = "Fy." + toString(IJ) + ".bin";

	//string sqrt_T_ff_V_filename = "sqrt_T_ff_V." + toString(IJ) + ".bin";
	string Gamma_yf_filename = "Gamma_yf." + toString(IJ) + ".bin";
	string Gamma_xf_filename = "Gamma_xf." + toString(IJ) + ".bin";

	string D_ux_filename = "D_ux." + toString(IJ) + ".bin";
	string D_uy_filename = "D_uy." + toString(IJ) + ".bin";

	string D_vx_filename = "D_vx." + toString(IJ) + ".bin";
	string D_vy_filename = "D_vy." + toString(IJ) + ".bin";

	string DTx_filename = "DTx." + toString(IJ) + ".bin";
	string DTy_filename = "DTy." + toString(IJ) + ".bin";


	WriteMassiveToBinaryFile(Fx_filename.c_str(), Fx, Na);
	WriteMassiveToBinaryFile(Fy_filename.c_str(), Fy, Na);

	//WriteMassiveToBinaryFile(sqrt_T_ff_V_filename.c_str(), sqrt_T_ff_V, Na);
	WriteMassiveToBinaryFile(Gamma_yf_filename.c_str(), Gamma_yf, Na);
	WriteMassiveToBinaryFile(Gamma_xf_filename.c_str(), Gamma_xf, Na);

	WriteMassiveToBinaryFile(D_ux_filename.c_str(), D_ux, Na);
	WriteMassiveToBinaryFile(D_uy_filename.c_str(), D_uy, Na);

	WriteMassiveToBinaryFile(D_vx_filename.c_str(), D_vx, Na);
	WriteMassiveToBinaryFile(D_vy_filename.c_str(), D_vy, Na);

	WriteMassiveToBinaryFile(DTx_filename.c_str(), DTx, Na);
	WriteMassiveToBinaryFile(DTy_filename.c_str(), DTy, Na);
#else
	WriteMassiveToBinaryFile("Fx.bin", Fx, Na);
	WriteMassiveToBinaryFile("Fy.bin", Fy, Na);

	//WriteMassiveToBinaryFile("sqrt_T_ff_V.bin", sqrt_T_ff_V, Na);
	WriteMassiveToBinaryFile("Gamma_yf.bin", Gamma_yf, Na);
	WriteMassiveToBinaryFile("Gamma_xf.bin", Gamma_xf, Na);

	WriteMassiveToBinaryFile("D_ux.bin", D_ux, Na);
	WriteMassiveToBinaryFile("D_uy.bin", D_uy, Na);

	WriteMassiveToBinaryFile("D_vx.bin", D_vx, Na);
	WriteMassiveToBinaryFile("D_vy.bin", D_vy, Na);

	WriteMassiveToBinaryFile("DTx.bin", DTx, Na);
	WriteMassiveToBinaryFile("DTy.bin", DTy, Na);
#endif //End of COMPILE_MPI
#endif //End of DEBUG_PROGRAM


#ifdef Calculate_p_and_Temper_using_Jakoby_method
	memcpy(p_old, p, Na * sizeof (double));
	memcpy(Temper_old, Temper, Na * sizeof (double));
#endif //End of Calculate_p_and_Temper_using_Jakoby_method

	for(counter = 0; counter < N_I; counter++)
	{
		max_residual_in_p[counter] = 0;
		max_residual_in_T[counter] = 0;

		max_residual_in_p_Iter_0[counter] = 0;
		max_residual_in_T_Iter_0[counter] = 0;

		max_Iter_p_c[counter] = 0;
	}

#ifdef COMPILE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif //End of COMPILE_MPI

	////define times for start of program
	time_solve_start = time(NULL);
	do
	{
		//This is optimized for speed - here the change pointers is not applicable
		memcpy(p_pr, p, Na * sizeof (double));

		memcpy(u_pr, u, Na * sizeof (double));
		memcpy(v_pr, v, Na * sizeof (double));

		memcpy(Temper_pr, Temper, Na * sizeof (double));
		memcpy(rho_pr, rho, Na * sizeof (double));

		for(counter = 0; counter < Iter; counter++)
		{
			max_residual_in_p[counter] = 0;
			max_residual_in_T[counter] = 0;

			max_residual_in_p_Iter_0[counter] = 0;
			max_residual_in_T_Iter_0[counter] = 0;

			max_Iter_p_c[counter] = 0;
		}

		//it++;
		it += ht;
		it_counter++;

		Iter_pr2 = Iter_pr1;
		Iter_pr1 = Iter;

		Iter = 0;
		continue_iter = true;

		//Iteration process
		do
		{
			Su_abs_max = 0;
			//First are calculations of u_pseudo on BD with subdomains

			//u_pseudo - Boundary Conditions and checks in loop
			for(counter = 0; counter < Na_bc_u; counter++)
			{
				//take node frome rule vectors
				node = rule_vector_ij_bc_u[counter];
				j = rule_vector_j_bc_u[counter];
				i = node - j * Nx;

				if(Periodic_boundary_conditions_about_OX)
				{
					ij_function_PeriodicBC(i, j);
				}
				else
				{
					ij_function(node);
				}


				//solve convective terms
				F_xu1 = 0.5 * (Fx[i_1j] + Fx[ij]);
				F_xu2 = 0.5 * (Fx[ij] + Fx[i1j]);


				////solve source term
				//sqrt_T != const
				Scu = //0;
					cT * (vol_inf_fb_bool[i1j] * u[i1j] * sqrt_T[ij] * hy[j] / (hx[i] * 3.0)
					+ vol_inf_fb_bool[i_2j] * u[i_1j] * sqrt_T[i_1j] * hy[j] / (hx[i_1] * 3.0)
					+ Gamma_yf[ij1] * (v[ij1] * vol_inf_fb_bool[ij1] - v[i_1j1] * vol_inf_fb_bool[i_1j1])
					- Gamma_yf[ij] * (v[ij] * vol_inf_fb_bool[ij_1] - v[i_1j] * vol_inf_fb_bool[i_1j_1])
					- sqrt_T[ij] * (v[ij1] * vol_inf_fb_bool[ij1] - v[ij] * vol_inf_fb_bool[ij_1]) * 2.0 / 3.0
					+ sqrt_T[i_1j] * (v[i_1j1] * vol_inf_fb_bool[i_1j1] - v[i_1j] * vol_inf_fb_bool[i_1j_1]) * 2.0 / 3.0);

				Spu = //0;
					-(D_ux[ij] + D_ux[i1j]) / 3.0;


				////solve full coefficient for equation
				//Upwind
				au1 = maximum(0, F_xu1) + D_ux[ij];
				au2 = maximum(0, -F_xu2) + D_ux[i1j];
				au3 = 0.5 * (maximum(0, Fy[i_1j]) + maximum(0, Fy[ij])) + D_uy[ij];
				au4 = 0.5 * (maximum(0, -Fy[i_1j1]) + maximum(0, -Fy[ij1])) + D_uy[ij1];

				//divide ht
				bu3 = (rho_pr[i_1j] * hx[i_1] + rho_pr[ij] * hx[i]) * u_pr[ij] * hy[j] * 0.5 / ht
					+ Scu;


				//au01 = au0c + D_u1 + D_u2 + D_u3 + D_u4 + aut0 - Spu;
				au0 = au1 + au2 + au3 + au4
					+ 0.5 * (-Fy[i_1j] - Fy[ij] + Fy[i_1j1] + Fy[ij1]) + F_xu2 - F_xu1
					- Spu
					+ (rho[i_1j] * hx[i_1] + rho[ij] * hx[i]) * hy[j] * 0.5 / ht
					;

				u_pseudo[ij] = (au1 * vol_inf_fb_bool[i1j] * u[i_1j] //in case of the external corner in program is not included that fluid can enter the body
				                + au2 * vol_inf_fb_bool[i_2j] * u[i1j] //in case of the external corner in program is not included that fluid can enter the body

				                + au3 * u[ij_1] + au4 * u[ij1] + bu3
								//- cPch * (p_pr[ij] - p_pr[i_1j]) * hy[j]
								)
								/ au0;


				d_p_c_12[ij] = cPch * hy[j] / au0;


				//Solve coefficents for p
				hyj_upwind_rho_1_0 = hy[j] * rho_u[node];

				apx[node] = d_p_c_12[node] * hyj_upwind_rho_1_0;

				bpx[node] = vol_inf_fb_bool[i_1j] * vol_inf_fb_bool[node]
					* hyj_upwind_rho_1_0 * u_pseudo[node];

			}
#ifdef DEBUG_PROGRAM
#ifdef COMPILE_MPI
			string u_pseudo_filename = "u_pseudo." + toString(IJ) + ".bin";
			WriteMassiveToBinaryFile(u_pseudo_filename.c_str(), u_pseudo, Na);

			MPI_Barrier(MPI_COMM_WORLD);
			Write_massive_entire_computational_domain_MPI("u_pseudo", u_pseudo, double_type);
			MPI_Barrier(MPI_COMM_WORLD);
#else
			WriteMassiveToBinaryFile("u_pseudo.bin", u_pseudo, Na);
#endif //End of COMPILE_MPI
#endif //End of DEBUG_PROGRAM

			//The exchange data for d_p_c_12, apx, bpx, d_p_c_34, apy and bpy.
			//is needed to calculate p, u and v.
#ifdef COMPILE_MPI
			//Start data exchange for boundaries of d_p_c_12, apx and bpx.
			//StartDataExchangeWithNeighbourhoodsSubDomains_MPI(d_p_c_12, Nx_u, Ny_u);
			StartDataExchangeWithNeighbourhoodsSubDomains_MPI(apx, Nx_u, Ny_u);
			StartDataExchangeWithNeighbourhoodsSubDomains_MPI(bpx, Nx_u, Ny_u);
#endif //End of COMPILE_MPI



			//First are calculations of v_pseudo on BD with subdomains
			//optimized for speed code
			Sv_abs_max = 0;
			//Boundary Conditions and checks in loop
			for(counter = 0; counter < Na_bc_v; counter++)
			{
				//take node from rule vectors
				node = rule_vector_ij_bc_v[counter];
				j = rule_vector_j_bc_v[counter];
				i = node - j * Nx;


				if(Periodic_boundary_conditions_about_OX)
				{
					ij_function_PeriodicBC(i, j);
				}
				else
				{
					ij_function(node);
				}

				//solve convective terms
				F_yv3 = 0.5 * (Fy[ij_1] + Fy[ij]);
				F_yv4 = 0.5 * (Fy[ij] + Fy[ij1]);


				//solve source term
				//sqrt_T != const
				Scv = //0;
					cT * (vol_inf_fb_bool[ij1] * v[ij1] * sqrt_T[ij] * hx[i] / (hy[j] * 3.0)
					+ vol_inf_fb_bool[ij_1 - Nx * (0 <= (int)(j - 2))] * v[ij_1] * sqrt_T[ij_1] * hx[i] / (hy[j - 1 * (0 <= (int)(j - 1))] * 3.0)
					+ Gamma_xf[i1j] * (u[i1j] * vol_inf_fb_bool[i1j] - u[i1j_1] * vol_inf_fb_bool[i1j_1])
					- Gamma_xf[ij] * (u[ij] * vol_inf_fb_bool[i_1j] - u[ij_1] * vol_inf_fb_bool[i_1j_1])
					- sqrt_T[ij] * (u[i1j] * vol_inf_fb_bool[i1j] - u[ij] * vol_inf_fb_bool[i_1j]) * 2.0 / 3.0
					+ sqrt_T[ij_1] * (u[i1j_1] * vol_inf_fb_bool[i1j_1] - u[ij_1] * vol_inf_fb_bool[i_1j_1]) * 2.0 /3.0)
					- c_Fr_v * (rho[ij_1] * hx[i] * 0.5 * hy[j - 1 * (0 <= (int)(j - 1))] + rho[ij] * hx[i] * 0.5 * hy[j]) //gravity field
					;

				Spv = //0;
					-(D_vy[ij] + D_vy[ij1]) / 3.0;


				////solve full coefficient for equation
				//Upwind
				av1 = 0.5 * (maximum(0, Fx[ij]) + maximum(0, Fx[ij_1])) + D_vx[ij];
				av2 = 0.5 * (maximum(0, -Fx[i1j]) + maximum(0, -Fx[i1j_1])) + D_vx[i1j];
				av3 = maximum(0, F_yv3) + D_vy[ij];
				av4 = maximum(0, -F_yv4) + D_vy[ij1];

				//divide ht
				bv3 = (rho_pr[ij_1] * hy[j - 1] + rho_pr[ij] * hy[j]) * v_pr[ij] * hx[i] * 0.5 / ht
					+ Scv;

				//av01 = av0c + D_v1 + D_v2 + D_v3 + D_v4 + avt0 - Spv;
				av0 = av1 + av2 + av3 + av4
					+ 0.5 * (Fx[i1j] - Fx[ij] + Fx[i1j_1] -Fx[ij_1]) + F_yv4 - F_yv3
					- Spv
					+ (rho[ij_1] * hy[j - 1] + rho[ij] * hy[j]) * hx[i] * 0.5 / ht
					;


				v_pseudo[ij] = (av1 * v[i_1j] + av2 * v[i1j]
				                + av3 * vol_inf_fb_bool[ij_1 - Nx * (0 <= (int)(j - 2))] * v[ij_1] //in case of the external corner in program is not included that fluid can enter the body
				                + av4 * vol_inf_fb_bool[ij1] * v[ij1] //in case of the external corner in program is not included that fluid can enter the body
				                + bv3
								//- cPch * (p_pr[ij] - p_pr[ij_1]) * hx[i]
								)
								/ av0;


				d_p_c_34[ij] = cPch * hx[i] / av0;


				//Solve coefficients for p
				hxi_upwind_rho_3_0 = hx[i] * rho_v[node];

				apy[node] = d_p_c_34[node] * hxi_upwind_rho_3_0;

				bpy[node] = vol_inf_fb_bool[node - Nx] * vol_inf_fb_bool[node]
					* hxi_upwind_rho_3_0 * v_pseudo[node];

			}
#ifdef DEBUG_PROGRAM
#ifdef COMPILE_MPI
			string v_pseudo_filename = "v_pseudo." + toString(IJ) + ".bin";
			WriteMassiveToBinaryFile(v_pseudo_filename.c_str(), v_pseudo, Na);

			MPI_Barrier(MPI_COMM_WORLD);
			Write_massive_entire_computational_domain_MPI("v_pseudo", v_pseudo, double_type);
			MPI_Barrier(MPI_COMM_WORLD);
#else
			WriteMassiveToBinaryFile("v_pseudo.bin", v_pseudo, Na);
#endif //End of COMPILE_MPI
#endif //End of DEBUG_PROGRAM



			//The exchange data for d_p_c_12, apx, bpx, d_p_c_34, apy and bpy.
			//is needed to calculate p, u and v.
#ifdef COMPILE_MPI
			//Start data exchange for boundaries of d_p_c_34, apy and bpy.
			//StartDataExchangeWithNeighbourhoodsSubDomains_MPI(d_p_c_34, Nx_v, Ny_v);
			StartDataExchangeWithNeighbourhoodsSubDomains_MPI(apy, Nx_v, Ny_v);
			StartDataExchangeWithNeighbourhoodsSubDomains_MPI(bpy, Nx_v, Ny_v);
#endif //End of COMPILE_MPI


			//To mask latency make calculations
			//u_pseudo - calculate inside subdomain, with no BC and checks.
			for(counter = 0; counter < Na_nch_u; counter++)
			{
				//take node frome rule vectors
				node = rule_vector_ij_nch_u[counter];
				j = rule_vector_j_nch_u[counter];
				i = node - j * Nx;

				ij_function(node);


				//solve convective terms
				F_xu1 = 0.5 * (Fx[i_1j] + Fx[ij]);
				F_xu2 = 0.5 * (Fx[ij] + Fx[i1j]);


				//solve source term
				//sqrt_T != const
				Scu = //0;
					cT * (u[i1j] * sqrt_T[ij] * hy[j] / (hx[i] * 3.0)
					+ u[i_1j] * sqrt_T[i_1j] * hy[j] / (hx[i_1] * 3.0)
					+ Gamma_yf[ij1] * (v[ij1] - v[i_1j1])
					- Gamma_yf[ij] * (v[ij] - v[i_1j])
					- sqrt_T[ij] * (v[ij1] - v[ij]) * 2.0 / 3.0
					+ sqrt_T[i_1j] * (v[i_1j1] - v[i_1j]) * 2.0 / 3.0);

				Spu = //0;
					-(D_ux[ij] + D_ux[i1j]) / 3.0;


				////solve full coefficient for equation
				//Upwind
				au1 = maximum(0, F_xu1) + D_ux[ij];
				au2 = maximum(0, -F_xu2) + D_ux[i1j];
				au3 = 0.5 * (maximum(0, Fy[i_1j]) + maximum(0, Fy[ij])) + D_uy[ij];
				au4 = 0.5 * (maximum(0, -Fy[i_1j1]) + maximum(0, -Fy[ij1])) + D_uy[ij1];


				//divide ht
				bu3 = (rho_pr[i_1j] * hx[i_1] + rho_pr[ij] * hx[i]) * u_pr[ij] * hy[j] * 0.5 / ht
					+ Scu;


				//au01 = au0c + D_u1 + D_u2 + D_u3 + D_u4 + aut0 - Spu;
				au0 = au1 + au2 + au3 + au4
					+ 0.5 * (-Fy[i_1j] - Fy[ij] + Fy[i_1j1] + Fy[ij1]) + F_xu2 - F_xu1
					- Spu
					+ (rho[i_1j] * hx[i_1] + rho[ij] * hx[i]) * hy[j] * 0.5 / ht
					;


				u_pseudo[ij] = (au1 * u[i_1j] + au2 * u[i1j]
					+ au3 * u[ij_1] + au4 * u[ij1] + bu3
					//- cPch * (p_pr[ij] - p_pr[i_1j]) * hy[j]
					)
					/ au0;


				d_p_c_12[ij] = cPch * hy[j] / au0;


				//Solve coefficents for p
				hyj_upwind_rho_1_0 = hy[j] * rho_u[node];

				apx[node] = d_p_c_12[node] * hyj_upwind_rho_1_0;

				bpx[node] = hyj_upwind_rho_1_0 * u_pseudo[node];

			}
#ifdef DEBUG_PROGRAM
#ifdef COMPILE_MPI
			WriteMassiveToBinaryFile(u_pseudo_filename.c_str(), u_pseudo, Na);

			MPI_Barrier(MPI_COMM_WORLD);
			Write_massive_entire_computational_domain_MPI("u_pseudo", u_pseudo, double_type);
			MPI_Barrier(MPI_COMM_WORLD);
#else
			WriteMassiveToBinaryFile("u_pseudo.bin", u_pseudo, Na);
#endif //End of COMPILE_MPI
#endif //End of DEBUG_PROGRAM



			//v_pseudo - calculate inside subdomain, with no BC and checks.
			for(counter = 0; counter < Na_nch_v; counter++)
			{
				//take node frome rule vectors
				node = rule_vector_ij_nch_v[counter];
				j = rule_vector_j_nch_v[counter];
				i = node - j * Nx;

				ij_function(node);


				//solve convective terms
				F_yv3 = 0.5 * (Fy[ij_1] + Fy[ij]);
				F_yv4 = 0.5 * (Fy[ij] + Fy[ij1]);


				//solve source term
				//sqrt_T != const
				Scv = //0;
					cT * (v[ij1] * sqrt_T[ij] * hx[i] / (hy[j] * 3.0)
					+ v[ij_1] * sqrt_T[ij_1] * hx[i] / (hy[j - 1] * 3.0)
					+ Gamma_xf[i1j] * (u[i1j] - u[i1j_1])
					- Gamma_xf[ij] * (u[ij] - u[ij_1])
					- sqrt_T[ij] * (u[i1j] - u[ij]) * 2.0 / 3.0
					+ sqrt_T[ij_1] * (u[i1j_1] - u[ij_1]) * 2.0 /3.0)
					- c_Fr_v * (rho[ij_1] * hx[i] * 0.5 * hy[j - 1] + rho[ij] * hx[i] * 0.5 * hy[j]) //gravity field
					;

				Spv = //0;
					-(D_vy[ij] + D_vy[ij1]) / 3.0;


				////solve full coefficient for equation
				//Upwind
				av1 = 0.5 * (maximum(0, Fx[ij]) + maximum(0, Fx[ij_1])) + D_vx[ij];
				av2 = 0.5 * (maximum(0, -Fx[i1j]) + maximum(0, -Fx[i1j_1])) + D_vx[i1j];
				av3 = maximum(0, F_yv3) + D_vy[ij];
				av4 = maximum(0, -F_yv4) + D_vy[ij1];


				//divide ht
				bv3 = (rho_pr[ij_1] * hy[j - 1] + rho_pr[ij] * hy[j]) * v_pr[ij] * hx[i] * 0.5 / ht
					+ Scv;

				//av01 = av0c + D_v1 + D_v2 + D_v3 + D_v4 + avt0 - Spv;
				av0 = av1 + av2 + av3 + av4
					+ 0.5 * (Fx[i1j] - Fx[ij] + Fx[i1j_1] - Fx[ij_1]) + F_yv4 - F_yv3
					- Spv
					+ (rho[ij_1] * hy[j - 1] + rho[ij] * hy[j]) * hx[i] * 0.5 / ht
					;


				v_pseudo[ij] = (av1 * v[i_1j] + av2 * v[i1j]
					+ av3 * v[ij_1] + av4 * v[ij1] + bv3
					//- cPch * (p_pr[ij] - p_pr[ij_1]) * hx[i]
					)
					/ av0;


				d_p_c_34[ij] = cPch * hx[i] / av0;


				//Solve coefficients for p
				hxi_upwind_rho_3_0 = hx[i] * rho_v[node];

				apy[node] = d_p_c_34[node] * hxi_upwind_rho_3_0;

				bpy[node] = hxi_upwind_rho_3_0 * v_pseudo[node];


			}
#ifdef DEBUG_PROGRAM
#ifdef COMPILE_MPI
			WriteMassiveToBinaryFile(v_pseudo_filename.c_str(), v_pseudo, Na);

			MPI_Barrier(MPI_COMM_WORLD);
			Write_massive_entire_computational_domain_MPI("v_pseudo", v_pseudo, double_type);
			MPI_Barrier(MPI_COMM_WORLD);
#else
			WriteMassiveToBinaryFile("v_pseudo.bin", v_pseudo, Na);

			WriteMassiveToBinaryFile("d_p_c_12.bin", d_p_c_12, Na);
			WriteMassiveToBinaryFile("apx.bin", apx, Na);
			WriteMassiveToBinaryFile("bpx.bin", bpx, Na);

			WriteMassiveToBinaryFile("d_p_c_34.bin", d_p_c_34, Na);
			WriteMassiveToBinaryFile("apy.bin", apy, Na);
			WriteMassiveToBinaryFile("bpy.bin", bpy, Na);
#endif //End of COMPILE_MPI
#endif //End of DEBUG_PROGRAM


			//Solve Temper coef
			//Solve matrix which are constant when we solve equation for T
			//No checks in loop
			for(counter = 0; counter < Na_nch_others; counter++)
			{
				//take node from rule vectors
				node = rule_vector_ij_nch_others[counter];
				j = rule_vector_j_nch_others[counter];
				i = node - j * Nx;

				ij_function(node);


				dudx_T = (u[i1j] - u[ij]) / hx[i];

				dvdy_T = (v[ij1] - v[ij]) / hy[j];


				//Without check for external corner, when are applied velocity slip boundary conditions
				//Integration of (dudy + dvdx)^2 , without multiplication of hx[i] * hy[j]
				//for subvolume 0
				#define dudy_dvdx_2_subvolume_0 (((u[ij] - u[ij_1]) / (hy[j - 1 * (0 < j)] + hy[j])) + ((v[ij] - v[i_1j]) / (hx[i_1] + hx[i])))
				//for subvolume 1
				#define dudy_dvdx_2_subvolume_1 (((u[i1j] - u[i1j_1]) / (hy[j - 1 * (0 < j)] + hy[j])) + ((v[i1j] - v[ij]) / (hx[i] + hx[i1])))
				//for subvolume 2
				#define dudy_dvdx_2_subvolume_2 (((u[ij1] - u[ij]) / (hy[j] + hy[j + 1 * (j < (Nx - 1))])) + ((v[ij1] - v[i_1j1]) / (hx[i_1] + hx[i])))
				//for subvolume 3
				#define dudy_dvdx_2_subvolume_3 (((u[i1j1] - u[i1j]) / (hy[j] + hy[j + 1 * (j < (Nx - 1))])) + ((v[i1j1] - v[ij1]) / (hx[i] + hx[i1])))



				//This is final - here we multiplicate by ht
				ScT1[ij] =
					ht * sqrt_T[ij] * cT2 * hx[i] * hy[j] *
					(2.0 * (pow(dudx_T, 2) + pow(dvdy_T, 2))

					+ pow(dudy_dvdx_2_subvolume_0, 2)
					+ pow(dudy_dvdx_2_subvolume_1, 2)
					+ pow(dudy_dvdx_2_subvolume_2, 2)
					+ pow(dudy_dvdx_2_subvolume_3, 2)

					- pow((dudx_T + dvdy_T), 2) * 2.0 / 3.0)

					+ p_pr[ij] * hx[i] * hy[j];
			}


			//Temper coefficients on boundaries.
			//Boundary Conditions checks in loop.
			for(counter = 0; counter < Na_bc_others; counter++)
			{
				//take node from rule vectors
				node = rule_vector_ij_bc_others[counter];
				j = rule_vector_j_bc_others[counter];
				i = node - j * Nx;


				if(Periodic_boundary_conditions_about_OX)
				{
					ij_function_PeriodicBC(i, j);
				}
				else
				{
					ij_function(node);
				}



				dudx_T = (vol_inf_fb_bool[i1j] * u[i1j] - vol_inf_fb_bool[i_1j] * u[ij]) / hx[i];

				dvdy_T = (vol_inf_fb_bool[ij1] * v[ij1] - vol_inf_fb_bool[ij_1] * v[ij]) / hy[j];



				//Integration of (dudy + dvdx)^2 , without multiplication of hx[i] * hy[j]
				//for subvolume 0
				#define dudy_dvdx_2_subvolume_0_bc (((u[ij] - u[ij_1]) / ((vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[ij_1]) * hy[j - 1 * (0 < j)] + hy[j])) \
													+ ((v[ij] - v[i_1j]) / ((vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[i_1j]) * hx[i_1] + hx[i])))

				//for subvolume 1
				#define dudy_dvdx_2_subvolume_1_bc (((u[i1j] - u[i1j_1]) / ((vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[i1j_1]) * hy[j - 1 * (0 < j)] + hy[j])) \
													+ ((v[i1j] - v[ij]) / (hx[i] + (vol_inf_fb_bool[i1j_1] && vol_inf_fb_bool[i1j]) * hx[i1])))

				//for subvolume 2
				#define dudy_dvdx_2_subvolume_2_bc (((u[ij1] - u[ij]) / (hy[j] + (vol_inf_fb_bool[i_1j1] && vol_inf_fb_bool[ij1]) * hy[j + 1 * (j < (Nx - 1))])) \
														+ ((v[ij1] - v[i_1j1]) / ((vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[i_1j1]) * hx[i_1] + hx[i])))
				//for subvolume 3
				#define dudy_dvdx_2_subvolume_3_bc (((u[i1j1] - u[i1j]) / (hy[j] + (vol_inf_fb_bool[ij1] && vol_inf_fb_bool[i1j1]) * hy[j + 1 * (j < (Nx - 1))])) \
														+ ((v[i1j1] - v[ij1]) / (hx[i] + (vol_inf_fb_bool[i1j] && vol_inf_fb_bool[i1j1]) * hx[i1])))


				//This is final - here we multiplicate by ht
				ScT1[ij] =
					ht * sqrt_T[ij] * cT2 * hx[i] * hy[j] *
					(2.0 * (pow(dudx_T, 2) + pow(dvdy_T, 2))

					+ pow(dudy_dvdx_2_subvolume_0_bc, 2)
					+ pow(dudy_dvdx_2_subvolume_1_bc, 2)
					+ pow(dudy_dvdx_2_subvolume_2_bc, 2)
					+ pow(dudy_dvdx_2_subvolume_3_bc, 2)

					- pow((dudx_T + dvdy_T), 2) * 2.0 / 3.0)

					+ p_pr[ij] * hx[i] * hy[j];
			}
			//End of Temper coefficients on boundaries.


			//Calculate p and Temper
			//Initiate p_old and Temper_old
#ifdef Calculate_p_and_Temper_using_Jakoby_method
			memcpy(p_old, p, Na * sizeof (double));
			memcpy(Temper_old, Temper, Na * sizeof (double));
#endif //End of #ifdef Calculate_p_and_Temper_using_Jakoby_method


#ifdef COMPILE_MPI
			//Wait to complete data exchange of d_p_c_12, apx and bpx.
			//WaitToCompleteDataExchangeWithNeighbourhoodsSubDomains_MPI(d_p_c_12);
			WaitToCompleteDataExchangeWithNeighbourhoodsSubDomains_MPI(apx);
			WaitToCompleteDataExchangeWithNeighbourhoodsSubDomains_MPI(bpx);

			//Wait to complete data exchange of d_p_c_34, apy and bpy.
			//WaitToCompleteDataExchangeWithNeighbourhoodsSubDomains_MPI(d_p_c_34);
			WaitToCompleteDataExchangeWithNeighbourhoodsSubDomains_MPI(apy);
			WaitToCompleteDataExchangeWithNeighbourhoodsSubDomains_MPI(bpy);
#endif //End of COMPILE_MPI


#ifdef COMPILE_MPI
#ifdef USE_MPI_Barrier_in_iterative_process_to_be_more_stable
			//MPI_Barrier(MPI_COMM_WORLD);
#endif //End of #ifdef USE_MPI_Barrier_in_iterative_process_to_be_more_stable
#endif //End of COMPILE_MPI
			continue_iter_p_c = false;
			Iter_p_c = 0;
			do
			{
#ifdef Calculate_p_and_Temper_using_Jakoby_method
				//Calculate p and Temper using Jakoby method
				Change_two_pointers_double(p, p_old);
//				memcpy(p_old, p, Na * sizeof (double));

				Change_two_pointers_double(Temper, Temper_old);
//				memcpy(Temper_old, Temper, Na * sizeof (double));
#endif //End of Calculate_p_and_Temper_using_Jakoby_method

				max_residual_in_p[Iter] = 0;
				max_residual_in_T[Iter] = 0;



				//Calculate Swap boundary, Boundary Conditions and checks in loop
				for(counter = 0; counter < Na_bc_others; counter++)
				{
					//take node frome rule vectors
					node = rule_vector_ij_bc_others[counter];
					j = rule_vector_j_bc_others[counter];
					i = node - j * Nx;


					if(Periodic_boundary_conditions_about_OX)
					{
						ij_function_PeriodicBC(i, j);
					}
					else
					{
						ij_function(node);
					}


					////Temper
					//upwind
					aT1 = maximum(0, Fx[ij]) + DTx[ij];
					aT2 = maximum(0, -Fx[i1j]) + DTx[i1j];
					aT3 = maximum(0, Fy[ij]) + DTy[ij];
					aT4 = maximum(0, -Fy[ij1]) + DTy[ij1];

					//to save Temper[ij] to calculate r_Temper
#ifndef Calculate_p_and_Temper_using_Jakoby_method
					Temper_tmp = Temper[ij];
#endif //end of Calculate_p_and_Temper_using_Jakoby_method


					//This is final - here we multiplicate by ht
					Temper[ij] =
						(ht * (aT1 * (vol_inf_fb_bool[i_1j] * TEMPER[i_1j] // The control volume vol_inf_fb_bool[i_1j] is fluid
						              + (!vol_inf_fb_bool[i_1j]) * Temper_of_fluid_on_the_wall_u[ij]) // The control volume vol_inf_fb_bool[i_1j] is body

						       + aT2 * (vol_inf_fb_bool[i1j] * TEMPER[i1j] // The control volume vol_inf_fb_bool[i1j] is fluid
						               + (!vol_inf_fb_bool[i1j]) * Temper_of_fluid_on_the_wall_u[i1j]) // The control volume vol_inf_fb_bool[i1j] is body

							   + aT3 * (vol_inf_fb_bool[ij_1] * TEMPER[ij_1] // The control volume vol_inf_fb_bool[ij_1] is fluid
									   + (!vol_inf_fb_bool[ij_1]) * Temper_of_fluid_on_the_wall_v[ij]) // The control volume vol_inf_fb_bool[ij_1] is body

							   + aT4 * (vol_inf_fb_bool[ij1] * TEMPER[ij1] // The control volume vol_inf_fb_bool[ij1] is fluid
									   + (!vol_inf_fb_bool[ij1]) * Temper_of_fluid_on_the_wall_v[ij1]) // The control volume vol_inf_fb_bool[ij1] is body

								//+ sqrt(Temper_old[ij]) * ScT1[ij] + ScT2[ij] //this is bT0
								+ (1.0 - gamma1) * P[ij] * ((u[i1j] - u[ij]) * hy[j] + (v[ij1] - v[ij]) * hx[i]))

								+ ScT1[ij])
							/
							(ht * (aT1 + aT2 + aT3 + aT4
							+ Fx[i1j] - Fx[ij] + Fy[ij1] - Fy[ij])
							+ hx[i] * hy[j] * rho[ij]);


#ifdef Calculate_p_and_Temper_using_Jakoby_method
					r_Temper = Temper[ij] - Temper_old[ij];
#else
					r_Temper = Temper[ij] - Temper_tmp;
#endif //end of Calculate_p_and_Temper_using_Jakoby_method



					////Solve p
					////solve by method Gouse-Zaidel
					//Solve BC if we have Given_velocity_on_xb == true
					Given_velocity_on_xb_bool_tmp = Given_velocity_on_xb * (i == Nx_others[0]);

					//Solve BC if we have Given_velocity_on_xe == true
					Given_velocity_on_xe_bool_tmp = Given_velocity_on_xe * (i == Nx_others[1] - 1);


					//to save p[ij] to calculate r_p
#ifndef Calculate_p_and_Temper_using_Jakoby_method
					p_tmp = p[ij];
#endif //end of Calculate_p_and_Temper_using_Jakoby_method


					//This is final - here we multiplicate by ht
					p[ij] = (ht * (apx[ij] * P[i_1j] + apx[i1j] * P[i1j]
						+ apy[ij] * P[ij_1] + apy[ij1] * P[ij1]
						+ bpx[ij] - bpx[i1j] + bpy[ij] - bpy[ij1]
						+ Given_velocity_on_xb_bool_tmp * Fx[ij]
						- Given_velocity_on_xe_bool_tmp * Fx[i1j]) //this is bp0

						+ p_pr[ij] * hx[i] * hy[j] / (Temper_pr[ij])
						)
						/
						(ht * (apx[ij] + apx[i1j] + apy[ij] + apy[ij1])
						+ hx[i] * hy[j] / (Temper[ij])
						);


#ifdef Calculate_p_and_Temper_using_Jakoby_method
					r_p = p[ij] - p_old[ij];
#else
					r_p = p[ij] - p_tmp;
#endif //end of Calculate_p_and_Temper_using_Jakoby_method

					if(max_residual_in_T[Iter] < fabs(r_Temper))
					{
						max_residual_in_T[Iter] = fabs(r_Temper);
					}

					if(max_residual_in_p[Iter] < fabs(r_p))
					{
						max_residual_in_p[Iter] = fabs(r_p);
					}



				}
				////WriteMassiveToBinaryFile("Temper.bin", Temper, Na);
				////WriteMassiveToBinaryFile("p.bin", p, Na);

#ifdef DEBUG_PROGRAM
#ifdef COMPILE_MPI
				MPI_Barrier(MPI_COMM_WORLD);

				string Temper_filename = "Temper." + toString(IJ) + ".bin";
				string p_filename = "p." + toString(IJ) + ".bin";

				WriteMassiveToBinaryFile(Temper_filename.c_str(), Temper, Na);
				WriteMassiveToBinaryFile(p_filename.c_str(), p, Na);

				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);
#else
				WriteMassiveToBinaryFile("Temper.bin", Temper, Na);
				WriteMassiveToBinaryFile("p.bin", p, Na);
#endif //End of COMPILE_MPI
#endif //End of DEBUG_PROGRAM

#ifdef COMPILE_MPI
				//Start data exchange for boundaries of p and Temper
				StartDataExchangeWithNeighbourhoodsSubDomains_MPI(p, Nx_others, Ny_others);
				StartDataExchangeWithNeighbourhoodsSubDomains_MPI(Temper, Nx_others, Ny_others);
#endif //End of COMPILE_MPI


				//Now make calculations of fi0 (p) and fi1 (Temper) of internal
				//area to mask time for exchanging of data for p and Temper
				//on subdomain boundaries.
				for(counter = 0; counter < Na_nch_others; counter++)
				{
					//take node frome rule vectors
					node = rule_vector_ij_nch_others[counter];
					j = rule_vector_j_nch_others[counter];
					i = node - j * Nx;

					ij_function(node);

					////Temper
					//upwind
					aT1 = maximum(0, Fx[ij]) + DTx[ij];
					aT2 = maximum(0, -Fx[i1j]) + DTx[i1j];
					aT3 = maximum(0, Fy[ij]) + DTy[ij];
					aT4 = maximum(0, -Fy[ij1]) + DTy[ij1];



					//to save Temper[ij] to calculate r_Temper
#ifndef Calculate_p_and_Temper_using_Jakoby_method
					Temper_tmp = Temper[ij];
#endif //end of Calculate_p_and_Temper_using_Jakoby_method


					//This is final - here we multiplicate by ht
					Temper[ij] = //1.0;
							(ht * (aT1 * TEMPER[i_1j] + aT2 * TEMPER[i1j]
								+ aT3 * TEMPER[ij_1] + aT4 * TEMPER[ij1]
								//+ sqrt(Temper_old[ij]) * ScT1[ij] + ScT2[ij] //this is bT0
								+ (1.0 - gamma1) * P[ij] * ((u[i1j] - u[ij]) * hy[j] + (v[ij1] - v[ij]) * hx[i]))

								+ ScT1[ij])
								/
								(ht * (aT1 + aT2 + aT3 + aT4
								+ Fx[i1j] - Fx[ij] + Fy[ij1] - Fy[ij])
								+ hx[i] * hy[j] * rho[ij]);

#ifdef Calculate_p_and_Temper_using_Jakoby_method
					r_Temper = Temper[ij] - Temper_old[ij];
#else
					r_Temper = Temper[ij] - Temper_tmp;
#endif //end of Calculate_p_and_Temper_using_Jakoby_method



					////Solve p
					//to save p[ij] to calculate r_p
#ifndef Calculate_p_and_Temper_using_Jakoby_method
					p_tmp = p[ij];
#endif //end of Calculate_p_and_Temper_using_Jakoby_method


					//This is final - here we multiplicate by ht
					p[ij] = (ht * (apx[ij] * P[i_1j] + apx[i1j] * P[i1j]
						+ apy[ij] * P[ij_1] + apy[ij1] * P[ij1]
						+ bpx[ij] - bpx[i1j] + bpy[ij] - bpy[ij1]) //this is bp0
						+ p_pr[ij] * hx[i] * hy[j] / (Temper_pr[ij])
							)
						/
						(ht * (apx[ij] + apx[i1j] + apy[ij] + apy[ij1])
						+ hx[i] * hy[j] / (Temper[ij])
						);

#ifdef Calculate_p_and_Temper_using_Jakoby_method
					r_p = p[ij] - p_old[ij];
#else
					r_p = p[ij] - p_tmp;
#endif //end of Calculate_p_and_Temper_using_Jakoby_method


					if(max_residual_in_p[Iter] < fabs(r_p))
					{
						max_residual_in_p[Iter] = fabs(r_p);
					}

					if(max_residual_in_T[Iter] < fabs(r_Temper))
					{
						max_residual_in_T[Iter] = fabs(r_Temper);
					}
				}

#ifdef COMPILE_MPI
#ifndef FIXED_NUMBER_OF_ITERATIONS_LOOP_3
				//Start send to IJ == 0 the data for max_residual_in_p and max_residual_in_T
				StartSendDataForConvergenceCriterionsFromSubDomainsFor_pT_MPI();
#endif //End of FIXED_NUMBER_OF_ITERATIONS_LOOP_3

				//Wait to complieate data exchange for p and Temper
				WaitToCompleteDataExchangeWithNeighbourhoodsSubDomains_MPI(p);
				WaitToCompleteDataExchangeWithNeighbourhoodsSubDomains_MPI(Temper);
#endif //End of COMPILE_MPI

				//BoundaryConditions_Pressure();
				//BoundaryConditions_Temperature();


#ifdef DEBUG_PROGRAM
#ifdef COMPILE_MPI
				MPI_Barrier(MPI_COMM_WORLD);

//				string Temper_filename = "Temper." + toString(IJ) + ".bin";
//				string p_filename = "p." + toString(IJ) + ".bin";

				MPI_Barrier(MPI_COMM_WORLD);
				WriteMassiveToBinaryFile(Temper_filename.c_str(), Temper, Na);
				MPI_Barrier(MPI_COMM_WORLD);
				WriteMassiveToBinaryFile(p_filename.c_str(), p, Na);
				MPI_Barrier(MPI_COMM_WORLD);

				Write_massive_entire_computational_domain_MPI("Temper", Temper, double_type);
				MPI_Barrier(MPI_COMM_WORLD);

				Write_massive_entire_computational_domain_MPI("p", p, double_type);
				MPI_Barrier(MPI_COMM_WORLD);
#else
				WriteMassiveToBinaryFile("Temper.bin", Temper, Na);
				WriteMassiveToBinaryFile("p.bin", p, Na);
#endif //End of COMPILE_MPI
#endif //End of DEBUG_PROGRAM



#ifndef FIXED_NUMBER_OF_ITERATIONS_LOOP_3
#ifdef COMPILE_MPI
				WaitToCompleteRecvDataForConvergenceCriterionsFromSubDomains_for_pT_MPI();

				if(IJ == 0)
				{
					FindMaximumResudualsFromAllSubDomains_for_pT_MPI();

					if(Iter_p_c == 0)
					{
						max_residual_in_p_Iter_0[Iter] = max_residual_in_p[Iter];
						max_residual_in_T_Iter_0[Iter] = max_residual_in_T[Iter];
					}
				}
#endif //End of COMPILE_MPI
#endif //End of FIXED_NUMBER_OF_ITERATIONS_LOOP_3

#ifndef COMPILE_MPI
				if(Iter_p_c == 0)
				{
					max_residual_in_p_Iter_0[Iter] = max_residual_in_p[Iter];
					max_residual_in_T_Iter_0[Iter] = max_residual_in_T[Iter];
				}
#endif //End of COMPILE_MPI


				Iter_p_c++;


#ifdef FIXED_NUMBER_OF_ITERATIONS_LOOP_3
				continue_iter_p_c = (Iter_p_c < N_I_p_c);
#else
				continue_iter_p_c =
					(Iter_p_c < MinimalNumberOfIterationInLoopFor_p_and_Temper)
					|| ((Iter_p_c < N_I_p_c)
						&& ((MaxError_p_c < max_residual_in_p[Iter]) || (MaxError_T < max_residual_in_T[Iter]))
						);
#endif //End of FIXED_NUMBER_OF_ITERATIONS_LOOP_3


#ifndef FIXED_NUMBER_OF_ITERATIONS_LOOP_3
#ifdef COMPILE_MPI
				StartSendDataForConvergenceCriterionsToSubDomains_for_pT_MPI();
#endif //End of COMPILE_MPI
#endif //End of FIXED_NUMBER_OF_ITERATIONS_LOOP_3



#ifndef FIXED_NUMBER_OF_ITERATIONS_LOOP_3
#ifdef COMPILE_MPI
				WaitToCompleteRecvDataForConvergenceCriterionsFromIJ_0_for_pT_MPI();
#endif //End of COMPILE_MPI
#endif //End of FIXED_NUMBER_OF_ITERATIONS_LOOP_3

#ifdef COMPILE_MPI
#ifdef USE_MPI_Barrier_in_iterative_process_to_be_more_stable_loop3
				MPI_Barrier(MPI_COMM_WORLD);
#endif //End of #ifdef USE_MPI_Barrier_in_iterative_process_to_be_more_stable_loop3
#endif //End of COMPILE_MPI
				//WriteSolvedDataToFile();
			}
			while(continue_iter_p_c);

			max_Iter_p_c[Iter] = Iter_p_c;


#ifdef DEBUG_PROGRAM
#ifdef COMPILE_MPI
			string Temper_filename = "Temper." + toString(IJ) + ".bin";
			string p_filename = "p." + toString(IJ) + ".bin";
			string rho_filename = "rho." + toString(IJ) + ".bin";

			WriteMassiveToBinaryFile(Temper_filename.c_str(), Temper, Na);
			WriteMassiveToBinaryFile(p_filename.c_str(), p, Na);
			WriteMassiveToBinaryFile(rho_filename.c_str(), rho, Na);


			MPI_Barrier(MPI_COMM_WORLD);
			Write_massive_entire_computational_domain_MPI("Temper", Temper, double_type);
			MPI_Barrier(MPI_COMM_WORLD);

			Write_massive_entire_computational_domain_MPI("p", p, double_type);
			MPI_Barrier(MPI_COMM_WORLD);

			Write_massive_entire_computational_domain_MPI("rho", rho, double_type);
			MPI_Barrier(MPI_COMM_WORLD);
#else
			WriteMassiveToBinaryFile("Temper.bin", Temper, Na);
			WriteMassiveToBinaryFile("p.bin", p, Na);
			WriteMassiveToBinaryFile("rho.bin", rho, Na);
#endif //End of COMPILE_MPI
#endif //End of DEBUG_PROGRAM

#ifdef COMPILE_MPI
#ifdef FIXED_NUMBER_OF_ITERATIONS_LOOP_3
			//Start send to IJ == 0 the data for max_residual_in_p and max_residual_in_T
			StartSendDataForConvergenceCriterionsFromSubDomainsFor_pT_MPI();
#endif //End of FIXED_NUMBER_OF_ITERATIONS_LOOP_3
#endif //End of COMPILE_MPI

			//Start to calculate u and v
			max_residual_in_u = 0;
			//Boundary Conditions checks in loop
			for(counter = 0; counter < Na_bc_u; counter++)
			{
				//take node frome rule vectors
				node = rule_vector_ij_bc_u[counter];
				j = rule_vector_j_bc_u[counter];
				i = node - j * Nx;

				if(Periodic_boundary_conditions_about_OX)
				{
					ij_function_PeriodicBC(i, j);
				}
				else
				{
					ij_function(node);
				}


				u_tmp = u[node];

				u[node] = u_pseudo[node]
					- (p[node] - p[i_1j]) * d_p_c_12[node];

				r_u	= u[node] - u_tmp;

				if(max_residual_in_u < fabs(r_u))
					max_residual_in_u = fabs(r_u);
			}
#ifdef DEBUG_PROGRAM
#ifdef COMPILE_MPI
			WriteMassiveToBinaryFile(u_filename.c_str(), u, Na);
#else
			WriteMassiveToBinaryFile("u.bin", u, Na);
#endif //End of COMPILE_MPI
#endif //End of DEBUG_PROGRAM


#ifdef COMPILE_MPI
			//Start data exchange for boundaries of u
			StartDataExchangeWithNeighbourhoodsSubDomains_MPI(u, Nx_u, Ny_u);
#endif //End of COMPILE_MPI

			max_residual_in_v = 0;
			//Boundary Conditions checks in loop
			for(counter = 0; counter < Na_bc_v; counter++)
			{
				//take node frome rule vectors
				node = rule_vector_ij_bc_v[counter];
				j = rule_vector_j_bc_v[counter];
				i = node - j * Nx;


				if(Periodic_boundary_conditions_about_OX)
				{
					ij_function_PeriodicBC(i, j);
				}
				else
				{
					ij_function(node);
				}


				v_tmp = v[node];

				v[node] = v_pseudo[node]
					- (p[node] - p[node - Nx]) * d_p_c_34[node];

				r_v	= v[node] - v_tmp;

				if(max_residual_in_v < fabs(r_v))
					max_residual_in_v = fabs(r_v);
			}
#ifdef DEBUG_PROGRAM
#ifdef COMPILE_MPI
			WriteMassiveToBinaryFile(v_filename.c_str(), v, Na);
#else
			WriteMassiveToBinaryFile("v.bin", v, Na);
#endif //End of COMPILE_MPI
#endif //End of DEBUG_PROGRAM


#ifdef COMPILE_MPI
			//Start data exchange for boundaries of v
			StartDataExchangeWithNeighbourhoodsSubDomains_MPI(v, Nx_v, Ny_v);
#endif //End of COMPILE_MPI


			//NO checks in loop
			for(counter = 0; counter < Na_nch_u; counter++)
			{
				//take node frome rule vectors
				node = rule_vector_ij_nch_u[counter];

				u_tmp = u[node];

				u[node] = u_pseudo[node]
					- (p[node] - p[node - 1]) * d_p_c_12[node];

				r_u	= u[node] - u_tmp;

				if(max_residual_in_u < fabs(r_u))
					max_residual_in_u = fabs(r_u);
			}
#ifdef DEBUG_PROGRAM
#ifdef COMPILE_MPI
			WriteMassiveToBinaryFile(u_filename.c_str(), u, Na);
#else
			WriteMassiveToBinaryFile("u.bin", u, Na);
#endif //End of COMPILE_MPI
#endif //End of DEBUG_PROGRAM


#ifdef COMPILE_MPI
			//Start send to IJ == 0 the data for max_residual_in_u
			StartSendDataForConvergenceCriterionsFromSubDomainsFor_u_MPI();
#endif //End of COMPILE_MPI


			//No checks in loop
			for(counter = 0; counter < Na_nch_v; counter++)
			{
				//take node frome rule vectors
				node = rule_vector_ij_nch_v[counter];

				v_tmp = v[node];

				v[node] = v_pseudo[node]
					- (p[node] - p[node - Nx]) * d_p_c_34[node];

				r_v	= v[node] - v_tmp;

				if(max_residual_in_v < fabs(r_v))
					max_residual_in_v = fabs(r_v);
			}
#ifdef DEBUG_PROGRAM
#ifdef COMPILE_MPI
			WriteMassiveToBinaryFile(v_filename.c_str(), v, Na);
#else
			WriteMassiveToBinaryFile("v.bin", v, Na);
#endif //End of COMPILE_MPI
#endif //End of DEBUG_PROGRAM


#ifdef COMPILE_MPI
			//Start send to IJ == 0 the data for max_residual_in_v
			StartSendDataForConvergenceCriterionsFromSubDomainsFor_v_MPI();
#endif //End of COMPILE_MPI


			//Make computations to mask time for data exchange.
			//Calculate sqrt_T
			for(counter = 0; counter < Na; counter++)
			{
				sqrt_T[counter] = sqrt(Temper[counter]);

				//solve rho
				if(vol_inf_fb_bool[counter]) rho[counter] = p[counter] / Temper[counter];
				else rho[counter] = 0;
			}
			//End of Make computations to mask time for data exchange.


#ifdef COMPILE_MPI
			//Wait to complete data exchange for u and v
			WaitToCompleteDataExchangeWithNeighbourhoodsSubDomains_MPI(u);
			WaitToCompleteDataExchangeWithNeighbourhoodsSubDomains_MPI(v);

#ifdef USE_MPI_Barrier_in_iterative_process_to_be_more_stable_loop2
			//MPI_Barrier(MPI_COMM_WORLD);
#endif //End of #ifdef USE_MPI_Barrier_in_iterative_process_to_be_more_stable_loop2
#endif //End of COMPILE_MPI


			//Make computations to mask time for data exchange.
			//use different type of non dimentional values
			if(TypeOfNonDimensiolization == NonDimentionalGiven_Re_pch_rho_V2_2)
			{
				cT = (1.0 / Re) * (u_given_xb_in_it * p_BC_xb / T_BC_xb);
				cPch = 0.5;
			}
			else if(TypeOfNonDimensiolization == NonDimentionalGiven_Re_pch_rho_R_T)
			{
				cT = (1.0 / Re) * (u_given_xb_in_it * p_BC_xb / T_BC_xb);
				cPch = u_given_xb_in_it * u_given_xb_in_it / (gamma1 * Ma * Ma);
			}
			else if(TypeOfNonDimensiolization == NonDimentionalGiven_Kn_pch_rho_R_T)
			{
				cT = Kn * 5.0 * sqrt(M_PI) / 16.0;
				cPch = 0.5;
			}

			//Calculate rho_u and rho_v after calculation of rho,
			//which is used in calculation of fluxes (Fx and Fy)
			//and procedure BoundaryConditions_Velocities()
			//Solve rho in middle points (rho_u and rho_v) if this is first
			//iteration after starting calculation
			Solve_rho_in_middle_points();


#ifdef DEBUG_PROGRAM
#ifdef COMPILE_MPI
			string rho_u_filename = "rho_u." + toString(IJ) + ".bin";
			string rho_v_filename = "rho_v." + toString(IJ) + ".bin";

			WriteMassiveToBinaryFile(rho_u_filename.c_str(), rho_u, Na);
			WriteMassiveToBinaryFile(rho_v_filename.c_str(), rho_v, Na);
#else
			WriteMassiveToBinaryFile("rho_u.bin", rho_u, Na);
			WriteMassiveToBinaryFile("rho_v.bin", rho_v, Na);
#endif //End of COMPILE_MPI
#endif //End of DEBUG_PROGRAM
			//End of Make computations to mask time for data exchange.

#ifdef FIXED_NUMBER_OF_ITERATIONS_LOOP_3
#ifdef COMPILE_MPI
				WaitToCompleteRecvDataForConvergenceCriterionsFromSubDomains_for_pT_MPI();

				if(IJ == 0)
				{
					FindMaximumResudualsFromAllSubDomains_for_pT_MPI();

					max_residual_in_p_Iter_0[Iter] = max_residual_in_p[Iter];
					max_residual_in_T_Iter_0[Iter] = max_residual_in_T[Iter];
				}
#endif //End of COMPILE_MPI
#endif //End of FIXED_NUMBER_OF_ITERATIONS_LOOP_3


#ifdef COMPILE_MPI
			WaitToCompleteRecvDataForConvergenceCriterionsFromSubDomains_for_uv_MPI();
#endif //End of COMPILE_MPI



			Iter++;
			//N_I_check_counter++;
			//Check for convergence.
			//This is made only from IJ == 0, when program is executing on MPI.
#ifdef COMPILE_MPI
			//Solve convergent criterions
			//Solve convergent criterions after first iteration

			//Only IJ == 0 check for criterions and send continue_iter to other processes
			if(IJ == 0)
#endif //End of COMPILE_MPI
			{
#ifdef COMPILE_MPI
				FindMaximumResudualsFromAllSubDomains_for_uv_MPI();
#endif //End of COMPILE_MPI

				if(Iter == 1)
				{
					max_residual_in_u_Iter_0 = max_residual_in_u;
					max_residual_in_v_Iter_0 = max_residual_in_v;
				}

				//N_I_check_counter = 0;

				//check for convergence
//				continue_iter = ((Iter < N_I)/* && (it_counter <= 1)*/);

				continue_iter = (Iter < N_I
								&& (MaxError_Velocities < max_residual_in_u
									|| MaxError_Velocities < max_residual_in_v
									|| MaxError_p_c < max_residual_in_p[Iter - 1]
									|| MaxError_T < max_residual_in_T[Iter - 1]));


				//The creterian to continue or not time step is calculated here, because
				//here have calculation to mask data exchange, when using MPI.
				countinue_time_step = (it_counter < Nt
										&& time(NULL) < time_solve_stop

										&& ((!(Iter <= 1 && Iter_p_c <= MinimalNumberOfIterationInLoopFor_p_and_Temper))
											|| MaxError_to_stop_program < max_residual_in_u
											|| MaxError_to_stop_program < max_residual_in_v
											|| MaxError_to_stop_program < max_residual_in_p_Iter_0[Iter - 1]
											|| MaxError_to_stop_program < max_residual_in_T_Iter_0[Iter - 1]));
			}


#ifdef COMPILE_MPI
			//MPI_Barrier(MPI_COMM_WORLD);
			StartSendDataForConvergenceCriterionsToSubDomains_for_it_loop_MPI();
#endif //End of COMPILE_MPI

			//Make computations to mask time for data exchange.
			//rho_u have to be calculated
			BoundaryConditions_Velocities_SlipBC_u();
#ifdef DEBUG_PROGRAM
#ifdef COMPILE_MPI
			string u_filename = "u." + toString(IJ) + ".bin";
			WriteMassiveToBinaryFile(u_filename.c_str(), u, Na);

			MPI_Barrier(MPI_COMM_WORLD);
			Write_massive_entire_computational_domain_MPI("u", u, double_type);
			MPI_Barrier(MPI_COMM_WORLD);
#else
			WriteMassiveToBinaryFile("u.bin", u, Na);
#endif //End of COMPILE_MPI
#endif //End of DEBUG_PROGRAM

			//rho_v have to be calculated
			BoundaryConditions_Velocities_SlipBC_v();
#ifdef DEBUG_PROGRAM
#ifdef COMPILE_MPI
			string v_filename = "v." + toString(IJ) + ".bin";
			WriteMassiveToBinaryFile(v_filename.c_str(), v, Na);

			MPI_Barrier(MPI_COMM_WORLD);
			Write_massive_entire_computational_domain_MPI("v", v, double_type);
			MPI_Barrier(MPI_COMM_WORLD);
#else
			WriteMassiveToBinaryFile("v.bin", v, Na);
#endif //End of COMPILE_MPI
#endif //End of DEBUG_PROGRAM

			//Boundary Conditions
			BoundaryConditions_Velocities();
#ifdef DEBUG_PROGRAM
#ifdef COMPILE_MPI
			WriteMassiveToBinaryFile(u_filename.c_str(), u, Na);
			WriteMassiveToBinaryFile(v_filename.c_str(), v, Na);

			MPI_Barrier(MPI_COMM_WORLD);
			Write_massive_entire_computational_domain_MPI("u", u, double_type);
			MPI_Barrier(MPI_COMM_WORLD);
			Write_massive_entire_computational_domain_MPI("v", v, double_type);
			MPI_Barrier(MPI_COMM_WORLD);
#else
			WriteMassiveToBinaryFile("u.bin", u, Na);
			WriteMassiveToBinaryFile("v.bin", v, Na);
#endif //End of COMPILE_MPI
#endif //End of DEBUG_PROGRAM

			//apply BC for correct pressure when we have given
			//velocity in x_b on pressure driven flow
			BoundaryConditions_Pressure_begin_iteration_for_it();

			//Apply temperature jump BC, after calculation of Temper and rho.
			BoundaryConditions_Temperature_TemperatureJumpBC();
#ifdef DEBUG_PROGRAM
#ifdef COMPILE_MPI
			string Temper_of_fluid_on_the_wall_u_filename = "Temper_of_fluid_on_the_wall_u." + toString(IJ) + ".bin";
			string Temper_of_fluid_on_the_wall_v_filename = "Temper_of_fluid_on_the_wall_v." + toString(IJ) + ".bin";

			WriteMassiveToBinaryFile(Temper_of_fluid_on_the_wall_u_filename.c_str(), Temper_of_fluid_on_the_wall_u, Na);
			WriteMassiveToBinaryFile(Temper_of_fluid_on_the_wall_v_filename.c_str(), Temper_of_fluid_on_the_wall_v, Na);
#else
			WriteMassiveToBinaryFile("Temper_of_fluid_on_the_wall_u.bin", Temper_of_fluid_on_the_wall_u, Na);
			WriteMassiveToBinaryFile("Temper_of_fluid_on_the_wall_v.bin", Temper_of_fluid_on_the_wall_v, Na);
#endif //End of COMPILE_MPI
#endif //End of DEBUG_PROGRAM

			BoundaryConditions_Temperature();
			BoundaryConditions_Pressure();

			//Solve fluxes (Fx and Fy) and diffusion terms if this is first iteration after
			//starting calculation
			SolveFluxes();
			SolveDiffisionTerms();

#ifdef DEBUG_PROGRAM
#ifdef COMPILE_MPI
			string Fx_filename = "Fx." + toString(IJ) + ".bin";
			string Fy_filename = "Fy." + toString(IJ) + ".bin";

			//string sqrt_T_ff_V_filename = "sqrt_T_ff_V." + toString(IJ) + ".bin";
			string Gamma_yf_filename = "Gamma_yf." + toString(IJ) + ".bin";
			string Gamma_xf_filename = "Gamma_xf." + toString(IJ) + ".bin";

			string D_ux_filename = "D_ux." + toString(IJ) + ".bin";
			string D_uy_filename = "D_uy." + toString(IJ) + ".bin";

			string D_vx_filename = "D_vx." + toString(IJ) + ".bin";
			string D_vy_filename = "D_vy." + toString(IJ) + ".bin";

			string DTx_filename = "DTx." + toString(IJ) + ".bin";
			string DTy_filename = "DTy." + toString(IJ) + ".bin";


			WriteMassiveToBinaryFile(Fx_filename.c_str(), Fx, Na);
			WriteMassiveToBinaryFile(Fy_filename.c_str(), Fy, Na);

			//WriteMassiveToBinaryFile(sqrt_T_ff_V_filename.c_str(), sqrt_T_ff_V, Na);
			WriteMassiveToBinaryFile(Gamma_yf_filename.c_str(), Gamma_yf, Na);
			WriteMassiveToBinaryFile(Gamma_xf_filename.c_str(), Gamma_xf, Na);

			WriteMassiveToBinaryFile(D_ux_filename.c_str(), D_ux, Na);
			WriteMassiveToBinaryFile(D_uy_filename.c_str(), D_uy, Na);

			WriteMassiveToBinaryFile(D_vx_filename.c_str(), D_vx, Na);
			WriteMassiveToBinaryFile(D_vy_filename.c_str(), D_vy, Na);

			WriteMassiveToBinaryFile(DTx_filename.c_str(), DTx, Na);
			WriteMassiveToBinaryFile(DTy_filename.c_str(), DTy, Na);
#else
			WriteMassiveToBinaryFile("Fx.bin", Fx, Na);
			WriteMassiveToBinaryFile("Fy.bin", Fy, Na);

			//WriteMassiveToBinaryFile("sqrt_T_ff_V.bin", sqrt_T_ff_V, Na);
			WriteMassiveToBinaryFile("Gamma_yf.bin", Gamma_yf, Na);
			WriteMassiveToBinaryFile("Gamma_xf.bin", Gamma_xf, Na);

			WriteMassiveToBinaryFile("D_ux.bin", D_ux, Na);
			WriteMassiveToBinaryFile("D_uy.bin", D_uy, Na);

			WriteMassiveToBinaryFile("D_vx.bin", D_vx, Na);
			WriteMassiveToBinaryFile("D_vy.bin", D_vy, Na);

			WriteMassiveToBinaryFile("DTx.bin", DTx, Na);
			WriteMassiveToBinaryFile("DTy.bin", DTy, Na);
#endif //End of COMPILE_MPI
#endif //End of DEBUG_PROGRAM
			//End of Make computations to mask time for data exchange.

			if(ToSolveContinuityEquation)
			{
				SolveContinuityEquation();
			}

#ifdef COMPILE_MPI
			WaitToCompleteRecvDataForConvergenceCriterionsFromIJ_0_for_it_loop_MPI();

#ifdef USE_MPI_Barrier_in_iterative_process_to_be_more_stable
			if(continue_iter)
			{
				MPI_Barrier(MPI_COMM_WORLD);
			}
#endif //End of #ifdef USE_MPI_Barrier_in_iterative_process_to_be_more_stable
#endif //End of COMPILE_MPI

#ifdef DEBUG_PROGRAM
			WriteSolvedDataToFile();
#endif

		}while(continue_iter);


		if(use_interpolation_for_u_BC_by_time)
		{
#ifdef COMPILE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif //End of COMPILE_MPI

			BoundaryConditions_function_of_time();

#ifdef COMPILE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif //End of COMPILE_MPI
		}


		Nt_save_solved_data_counter++;
		if(Nt_save_solved_data_counter >= Nt_save_solved_data)
		{
#ifdef COMPILE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif //End of COMPILE_MPI

			WriteSolvedDataToFile();

			//Solve Continuity Equation and write data to files
			if(ToSolveContinuityEquation)
			{
#ifdef COMPILE_MPI
				MPI_Barrier(MPI_COMM_WORLD);
#endif //End of COMPILE_MPI

				SolveContinuityEquation();
			}


			Nt_save_solved_data_counter = 0;

#ifdef COMPILE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif //End of COMPILE_MPI
		}


		Nt_save_DragCoefficient_counter++;
		if(Nt_save_DragCoefficient_counter >= Nt_save_DragCoefficient)
		{
#ifdef COMPILE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif //End of COMPILE_MPI

			WriteDragCoefficient();
			Nt_save_DragCoefficient_counter = 0;

#ifdef COMPILE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif //End of COMPILE_MPI
		}


		Nt_save_temp_data_counter++;
		if(Nt_save_temp_data_counter >= Nt_save_temp_data)
		{
#ifdef COMPILE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif //End of COMPILE_MPI

			WriteTempDataToFile();

#ifdef COMPILE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif //End of COMPILE_MPI

			WriteIterationsInLoopsToFile();
			Nt_save_temp_data_counter = 0;
		}


#ifdef COMPILE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif //End of COMPILE_MPI

	}while(countinue_time_step);


	//check that solve is finished
	solve_is_finished = Iter <= 1 && Iter_p_c <= 1
		&& max_residual_in_u < MaxError_to_stop_program
		&& max_residual_in_v < MaxError_to_stop_program
		&& max_residual_in_p[Iter - 1] < MaxError_to_stop_program
		&& max_residual_in_T[Iter - 1] < MaxError_to_stop_program;

	//Solve work time.
	max_time_for_solve = (time(NULL) - time_solve_start) / 3600.0;

	//write all final data in files
	WriteSolvedDataToFile();

	//Solve Continuity Equation and write data to files
	if(ToSolveContinuityEquation)
	{
		SolveContinuityEquation();
	}

	WriteDragCoefficient();
	WriteTempDataToFile();


#ifdef COMPILE_MPI
	//Wait all proceasses before terminate MPI
	MPI_Barrier(MPI_COMM_WORLD);

	//Terminate MPI
	MPI_Finalize();

#endif //End of COMPILE_MPI

	return 0;
}




//This is begin of file for class SIMPLEST_ANL_compr_T_h_var_NoDim_Algorithm



inline double hxfront(double& x)
{
	return(axf0 + axf1 * x + axf2 * x * x);
}

inline double hxback(double& x)
{
	return(axb0 + axb1 * x + axb2 * x * x + axb3 * x * x * x);
}

inline double hxfront_Polynom3(double& x)
{
	return(axf0 + axf1 * x + axf2 * x * x + axf3 * x * x * x);
}

inline double hxback_Polynom3(double& x)
{
	return(axb0 + axb1 * x + axb2 * x * x + axb3 * x * x * x);
}

inline double hx_all(double& x)
{
	double hx_temp;

	if(kind_of_mesh == Li_mesh) hx_temp = hx_all_linear(x);
	else if(kind_of_mesh == Par_mesh) hx_temp = hx_all_par(x);
	else if(kind_of_mesh == Li_mesh_flat) hx_temp = hx_all_linear_flat(x);
	else if(kind_of_mesh == Par_mesh_flat) hx_temp = hx_all_par_flat(x);
	else hx_temp = hx_all_Polynom3_flat_mesh(x);

	//hy_temp MUST be equal or small of hxmax -> hx_temp <= hxmax
	//and have to be bigger then hxmin
	if(hxmax < hx_temp) hx_temp = hxmax;
	else if (hx_temp < hxmin) hx_temp = hxmin;

	return hx_temp;
}

inline double hx_all_par(double& x)
{
	if(x < xfmax)
		return hxmax;
	else if((xfmax <= x) && (x <= xmid))
		return hxfront(x);
	else if((xmid < x) && (x <= xbmax))
		return hxback(x);
	else
		return hxmax;
}

inline double hx_all_par_flat(double& x)
{
	if(x < xfmax)
		return hxmax;
	else if((xfmax <= x) && (x <= xfmin))
		return hxfront(x);
	else if((xfmin < x) && (x <= xbmin))
		return (hxmin);
	else if((xbmin < x) && (x <= xbmax))
		return hxback(x);
	else
		return hxmax;
}

inline double hx_all_linear(double& x)
{
	if(x < xfmax)
		return (hxmax);
	else if((xfmax <= x) && (x <= xfmin))
		return (((hxmin - hxmax) / (xfmin - xfmax)) * (x - xfmax) + hxmax);
	else if((xfmin < x) && (x <= xmid))
		return (((hxmid - hxmin) / (xmid - xfmin)) * (x - xfmin) + hxmin);
	else if((xmid < x) && (x <= xbmin))
		return (((hxmin - hxmid) / (xbmin - xmid)) * (x - xmid) + hxmid);
	else if((xbmin < x) && (x <= xbmax))
		return (((hxmax - hxmin) / (xbmax - xbmin)) * (x - xbmin) + hxmin);
	else
		return (hxmax);
}

inline double hx_all_linear_flat(double& x)
{
	if(x < xfmax)
		return (hxmax);
	else if((xfmax <= x) && (x <= xfmin))
		return (((hxmin - hxmax) / (xfmin - xfmax)) * (x - xfmax) + hxmax);
	else if((xfmin < x) && (x <= xbmin))
		return (hxmin);
	else if((xbmin < x) && (x <= xbmax))
		return (((hxmax - hxmin) / (xbmax - xbmin)) * (x - xbmin) + hxmin);
	else
		return (hxmax);
}

inline double hx_all_Polynom3_flat_mesh(double& x)
{
	if(x < xfmax)
		return hxmax;
	else if((xfmax <= x) && (x <= xfmin))
		return hxfront_Polynom3(x);
	else if((xfmin < x) && (x <= xbmin))
		return (hxmin);
	else if((xbmin < x) && (x <= xbmax))
		return hxback_Polynom3(x);
	else
		return hxmax;
}


inline double hytop(double& y)
{
	return(ayt0 + ayt1 * y + ayt2 * y * y);
}

inline double hybottom(double& y)
{
	return(ayb0 + ayb1 * y + ayb2 * y * y);
}

inline double hytop_Polynom3(double& y)
{
	return(ayt0 + ayt1 * y + ayt2 * y * y + ayt3 * y * y * y);
}

inline double hybottom_Polynom3(double& y)
{
	return(ayb0 + ayb1 * y + ayb2 * y * y + ayb3 * y * y * y);
}

inline double hy_all(double& y)
{
	double hy_temp;

	if(kind_of_mesh == Li_mesh) hy_temp = hy_all_linear(y);
	else if(kind_of_mesh == Par_mesh) hy_temp = hy_all_par(y);
	else if(kind_of_mesh == Li_mesh_flat) hy_temp = hy_all_linear_flat(y);
	else if(kind_of_mesh == Par_mesh_flat) hy_temp = hy_all_par_flat(y);
	else hy_temp = hy_all_Polynom3_flat_mesh(y);

	//hy_temp MUST be equal or small of hymax -> hy_temp <= hymax
	//and have to be bigger then hymin
	if(hymax < hy_temp) hy_temp = hymax;
	else if (hy_temp < hymin) hy_temp = hymin;

	return hy_temp;
}

inline double hy_all_par(double& y)
{
	if(y < ybmax)
		return hymax;
	else if((ybmax <= y) && (y <= ymid))
		return hybottom(y);
	else if((ymid < y) && (y <= ytmax))
		return hytop(y);
	else
		return hymax;
}

inline double hy_all_par_flat(double& y)
{
	if(y < ybmax)
		return hymax;
	else if((ybmax <= y) && (y <= ybmin))
		return hybottom(y);
	else if((ybmin < y) && (y <= ytmin))
		return (hymin);
	else if((ytmin < y) && (y <= ytmax))
		return hytop(y);
	else
		return hymax;
}

inline double hy_all_linear(double& y)
{
	if(y < ybmax)
		return (hymax);
	else if((ybmax <= y) && (y <= ybmin))
		return (((hymin - hymax) / (ybmin - ybmax)) * (y - ybmax) + hymax);
	else if((ybmin < y) && (y <= ymid))
		return (((hymid - hymin) / (ymid - ybmin)) * (y - ybmin) + hymin);
	else if((ymid < y) && (y <= ytmin))
		return (((hymin - hymid) / (ytmin - ymid)) * (y - ymid) + hymid);
	else if((ytmin < y) && (y <= ytmax))
		return (((hymax - hymin) / (ytmax - ytmin)) * (y - ytmin) + hymin);
	else
		return (hymax);
}

inline double hy_all_linear_flat(double& y)
{
	if(y < ybmax)
		return (hymax);
	else if((ybmax <= y) && (y <= ybmin))
		return (((hymin - hymax) / (ybmin - ybmax)) * (y - ybmax) + hymax);
	else if((ybmin < y) && (y <= ytmin))
		return (hymin);
	else if((ytmin < y) && (y <= ytmax))
		return (((hymax - hymin) / (ytmax - ytmin)) * (y - ytmin) + hymin);
	else
		return (hymax);
}


inline double hy_all_Polynom3_flat_mesh(double& y)
{
	if(y < ybmax)
		return hymax;
	else if((ybmax <= y) && (y <= ybmin))
		return hybottom_Polynom3(y);
	else if((ybmin < y) && (y <= ytmin))
		return (hymin);
	else if((ytmin < y) && (y <= ytmax))
		return hytop_Polynom3(y);
	else
		return hymax;
}


void define_Nx_Ny(void)
{
	double x_e_temp;
	double x_e_TOL;//tolerance of the end of area on x: x_e_TOL = x_e - x_e_temp;
	x_e_TOL = 5 * hxmax;

	double x_e_TOL_temp;
	bool cont_solve_Nx;

	cont_solve_Nx = true;
	do
	{
		x_e_temp = x_b;
		for(counter = 0; counter < Nx; counter++)
		{
			x_e_temp = x_e_temp + hx_all(x_e_temp);
		}

		x_e_TOL_temp = x_e_temp - x_e;
		if(fabs(x_e_TOL_temp) > x_e_TOL)
		{
			//Nx must be corrected
			if(x_e_TOL_temp < 0)
			{
				//x_e_temp is inside area. Nx must be increased.
				Nx = Nx + int(fabs(x_e_TOL_temp / hxmax));
			}

			if(x_e_TOL_temp > 0)
			{
				//x_e_temp is inside area. Nx must be decreased.
				if((Nx - int(fabs(x_e_TOL_temp / hxmax))) > 2)
					Nx = Nx - int(fabs(x_e_TOL_temp / hxmax));
				else
					Nx = int(0.75 * Nx);
			}
		}
		else
		{
			cont_solve_Nx = false;
		}

	}while(cont_solve_Nx);


	double y_e_temp;
	double y_e_TOL;//tolerance of the end of area on y: y_e_TOL = y_e - y_e_temp;
	y_e_TOL = 5 * hxmax;

	double y_e_TOL_temp;
	bool cont_solve_Ny;

	cont_solve_Ny = true;
	do
	{
		y_e_temp = y_b;
		for(counter = 0; counter < Ny; counter++)
		{
			y_e_temp = y_e_temp + hy_all(y_e_temp);
		}

		y_e_TOL_temp = y_e_temp - y_e;
		if(fabs(y_e_TOL_temp) > y_e_TOL)
		{
			//Ny must be corrected
			if(y_e_TOL_temp < 0)
			{
				//y_e_temp is inside area. Ny must be increased.
				Ny = Ny + int(fabs(y_e_TOL_temp / hymax));
			}

			if(y_e_TOL_temp > 0)
			{
				//y_e_temp is inside area. Ny must be decreased.
				if((Ny - int(fabs(y_e_TOL_temp / hymax))) > 2)
					Ny = Ny - int(fabs(y_e_TOL_temp / hymax));
				else
					Ny = int(0.75 * Ny);
			}
		}
		else
		{
			cont_solve_Ny = false;
		}

	}while(cont_solve_Ny);
}


void WriteEnteredDataToFile(void)
{
	using namespace std;

	ofstream output_data;
	output_data.open("EnteredData.txt");

	if(output_data.is_open())
	{
		output_data
			<< std::setprecision(15)
			<< u_gas_nd << "	//u_gas_nd" << endl
			<< v_gas_nd << "	//v_gas_nd" << endl
			<< T_gas_nd << "	//T_gas_nd" << endl
			<< p_gas_nd << "	//p_gas_nd" << endl

			<< Ma << "	//Ma - solve in program   Ma = sqrt(6.0 / 5.0);" << endl
			<< Pr << "	//Pr - solve in program   Pr = 2.0 / 3.0;" << endl
			<< gamma1 << "	//gamma1 - solve in program   gamma1 = 5.0 / 3.0;" << endl
			<< c_mu << "	//c_mu - solve in program   c_mu = (5.0 / 16.0) * sqrt(2.0 * M_PI / gamma1);" << endl
			<< Kn << "	//Kn" << endl
			<< Fr << "	//Fr" << endl
			<< SolveGravityField << "	//SolveGravityField" << endl
			<< Re << "	//Re" << endl
			<< TypeOfNonDimensiolization << "	//TypeOfNonDimensiolization = 1 == NonDimentionalGiven_Re_pch_rho_V2_2, TypeOfNonDimensiolization = 2 == NonDimentionalGiven_Re_pch_rho_R_T, TypeOfNonDimensiolization = 3 == NonDimentionalGiven_Kn_pch_rho_R_T" << endl

			<< Nt << "	//Nt" << endl
#ifndef COMPILE_MPI
			<< Nx << "	//Nx" << endl
			<< Ny << "	//Ny" << endl
#endif
#ifdef COMPILE_MPI
			//When the computational domain is separated to subdomains the number of nodes of enire computational domain are contaning in
			//
			<< Nx_entire_computational_domain << "	//Nx" << endl
			<< Ny_entire_computational_domain << "	//Ny" << endl
#endif

			<< N_I << "	//N_I" << endl
			<< N_I_p_c << "	//N_I_p_c - if(N_I == 1) MUST N_I_p_c >= 2, if(N_I > 1) N_I_p_c can be 1" << endl
			<< N_I_T << "	//N_I_T - if(N_I == 1) MUST N_I_T >= 2, if(N_I > 1) N_I_T can be 1" << endl

			<< ht << "	//ht" << endl


			<< x_b << "	//x_b" << endl
			<< x_e << "	//x_e" << endl

			<< y_b << "	//y_b" << endl
			<< y_e << "	//y_e" << endl

			<< kind_of_mesh << "	//kind_of_mesh == Li_mesh = 1, Par_mesh = 2, Li_mesh_flat = 3 li, Par_mesh_flat = 4, Polynom3_flat_mesh = 5" << endl

			<< hxmin << "	//hxmin" << endl
			<< hxmax << "	//hxmax" << endl
			<< xfmin << "	//xfmin" << endl
			<< xfmax << "	//xfmax" << endl
			<< xbmin << "	//xbmin" << endl
			<< xbmax << "	//xbmax" << endl
			<< xmid << "	//xmid " << endl

			<< hymin << "	//hymin" << endl
			<< hymax << "	//hymax" << endl
			<< ytmin << "	//ytmin" << endl
			<< ytmax << "	//ytmax" << endl
			<< ybmin << "	//ybmin" << endl
			<< ybmax << "	//ybmax" << endl
			<< ymid << "	//ymid " << endl

			<< MaxError_Velocities << "	//MaxError_Velocities" << endl
			<< MaxError_p_c << "	//MaxError_p_c" << endl
			<< MaxError_T << "	//MaxError_T" << endl
			<< MaxError_to_stop_program << "	//MaxError_to_stop_program" << endl
			<< max_time_for_solve << "	//max_time_for_solve" << endl
			<< solve_is_finished << "	//solve_is_finished" << endl

			<< Nt_save_solved_data << "	//Nt_save_solved_data" << endl
			<< Nt_save_DragCoefficient << "	//Nt_save_DragCoefficient" << endl
			<< Nt_save_temp_data << "	//Nt_save_temp_data" << endl
			<< N_I_check << "	//N_I_check" << endl
			<< WriteDataInShortFormat_kind << "	//WriteDataInShortFormat_kind == 0; ==> to not write data in short format; //if WriteDataInShortFormat_kind == 1; ==> write: Re Kn CD CD_p u_gas_nd q_xb q_xe p_BC_xb p_BC_xe;" << endl

			<< ToReadSolvedDataFromFile << "	//ToReadSolvedDataFromFile" << endl
			<< ToReadSolvedDataFromBinaryFile << "	//ToReadSolvedDataFromBinaryFile" << endl
			<< ToWriteSolvedDataToBinaryFile << "	//ToWriteSolvedDataToBinaryFile" << endl
			<< ToWriteSolvedDataToNewFiles << "	//ToWriteSolvedDataToNewFiles" << endl
			<< ToContinueFromInterpolation << "	//ToContinueFromInterpolation" << endl

			<< ToImport_data_for_circles_from_file_b << "	//ToImport_data_for_circles_from_file_b" << endl
			<< ToImport_data_for_polyhedrons_from_file_b << "	//ToImport_data_for_polyhedrons_from_file_b" << endl

			<< ToImport_data_for_circles_from_file_p << "	//ToImport_data_for_circles_from_file_p" << endl
			<< ToImport_data_for_polyhedrons_from_file_p << "	//ToImport_data_for_polyhedrons_from_file_p" << endl

			<< ToImport_data_for_circles_from_file_V << "	//ToImport_data_for_circles_from_file_V" << endl
			<< ToImport_data_for_polyhedrons_from_file_V << "	//ToImport_data_for_polyhedrons_from_file_V" << endl


			<< Given_velocity_on_xb << "	//Given_velocity_on_xb" << endl
			<< Given_velocity_on_xe << "	//Given_velocity_on_xe" << endl
			<< Periodic_boundary_conditions_about_OX << "	//Periodic_boundary_conditions_about_OX" << endl
			<< Pressure_BC << "	//Pressure_BC" << endl
			<< ToStartFromExactSolutionForChannelFlowWithPressureBC << "	//ToStartFromExactSolutionForChannelFlowWithPressureBC" << endl

			<< p_BC_xb_correction_method << "	//p_BC_xb_correction_method == 0 --> no correction for p_BC_xb; p_BC_xb_correction_method == 1 == p_BC_xb_correction_method_Vmax_xb --> V maximum inflow == 1, p_BC_xb_correction_method == 2 == p_BC_xb_correction_method_Vmean_xb --> V mean inflow == 1" << endl
			<< Pressure_ratio_correct << "	//Pressure_ratio_correct == 0 --> not correct pressure ratio; Pressure_ratio_correct == 1 --> correct pressure ratio;" << endl
			<< GivenReKn << "	//GivenReKn -> 	u_given_xb = Re * Kn * sqrt(15.0 * M_PI / 128.0);" << endl
			<< correct_p_BC_xb << "	//correct_p_BC_xb" << endl
			<< u_given_xb << "	//u_given_xb" << endl
			<< u_given_xb_error << "	//u_given_xb_error" << endl
			<< p_correction_auto << "	//p_correction_auto" << endl
			<< p_correction_min << "	//p_correction_min" << endl
			<< p_correction_max << "	//p_correction_max" << endl
			<< p_correction << "	//p_correction" << endl

			<< p_BC_xb << "	//p_BC_xb" << endl
			<< T_BC_xb << "	//T_BC_xb" << endl

			<< p_BC_xe << "	//p_BC_xe" << endl
			<< T_BC_xe << "	//T_BC_xe" << endl

			<< dudx_0_BC_xb << "	//dudx_0_BC_xb" << endl
			<< durhodx_0_BC_xb << "	//durhodx_0_BC_xb" << endl
			<< drhodt_durhodx_dvrhody_0_BC_xb << "	//drhodt_durhodx_dvrhody_0_BC_xb" << endl
			<< ToCalculate_ui_1j_from_equation_for_uij_at_x_b << "	//ToCalculate_ui_1j_from_equation_for_uij_at_x_b" << endl

			<< dudx_0_BC_xe << "	//dudx_0_BC_xe" << endl
			<< durhodx_0_BC_xe << "	//durhodx_0_BC_xe" << endl
			<< drhodt_durhodx_dvrhody_0_BC_xe << "	//drhodt_durhodx_dvrhody_0_BC_xe" << endl
			<< ToCalculate_ui1j_from_equation_for_uij_at_x_e << "	//ToCalculate_ui1j_from_equation_for_uij_at_x_e" << endl
			<< u_uMin_Mout_0_BC_xe << "	//u_uMin_Mout_0_BC_xe" << endl

			<< dudy_0_BC_yb << "	//dudy_0_BC_yb" << endl
			<< dudy_0_BC_ye << "	//dudy_0_BC_ye" << endl

			<< dvdx_0_BC_xb << "	//dvdx_0_BC_xb" << endl
			<< ToCalculate_vi_1j_from_equation_for_vij_at_x_b << "	//ToCalculate_vi_1j_from_equation_for_vij_at_x_b" << endl

			<< dvdx_0_BC_xe << "	//dvdx_0_BC_xe" << endl
			<< ToCalculate_vi1j_from_equation_for_vij_at_x_e << "	//ToCalculate_vi1j_from_equation_for_vij_at_x_e" << endl

			<< dvdy_0_BC_yb << "	//dvdy_0_BC_yb" << endl
			<< dvdy_0_BC_ye << "	//dvdy_0_BC_ye" << endl

			<< dpdx_0_BC_xb << "	//dpdx_0_BC_xb" << endl
			<< ToCalculate_pi_1j_from_equation_for_pij_at_x_b << "	//ToCalculate_pi_1j_from_equation_for_pij_at_x_b" << endl

			<< dpdx_0_BC_xe << "	//dpdx_0_BC_xe" << endl
			<< ToCalculate_pi1j_from_equation_for_pij_at_x_e << "	//ToCalculate_pi1j_from_equation_for_pij_at_x_e" << endl

			<< dpdy_0_BC_yb << "	//dpdy_0_BC_yb" << endl
			<< dpdy_0_BC_ye << "	//dpdy_0_BC_ye" << endl

			<< dTdx_0_BC_xb << "	//dTdx_0_BC_xb" << endl
			<< ToCalculate_Temperi_1j_from_equation_for_Temperij_at_x_b << "	//ToCalculate_Temperi_1j_from_equation_for_Temperij_at_x_b" << endl

			<< dTdx_0_BC_xe << "	//dTdx_0_BC_xe" << endl
			<< ToCalculate_Temperi1j_from_equation_for_Temperij_at_x_e << "	//ToCalculate_Temperi1j_from_equation_for_Temperij_at_x_e" << endl

			<< dTdy_0_BC_yb << "	//dTdy_0_BC_yb" << endl
			<< dTdy_0_BC_ye << "	//dTdy_0_BC_ye" << endl

			<< To_Use_Kn_local_in_wall_BC << "	//To_Use_Kn_local_in_wall_BC" << endl

			<< ToSolveContinuityEquation << "	//ToSolveContinuityEquation" << endl

			<< use_interpolation_for_u_BC_by_time << "	//use_interpolation_for_u_BC_by_time" << endl
			<< t_BC_b << "	//t_BC_b" << endl
			<< t_BC_e << "	//t_BC_e" << endl
			<< u_BC_xb_t_BC_b << "	//u_BC_xb_t_BC_b" << endl
			<< u_BC_xb_t_BC_e << "	//u_BC_xb_t_BC_e" << endl

			<< w_u << "	//w_u" << endl
			<< w_v << "	//w_v" << endl
			<< w_p << "	//w_p" << endl
			<< w_Temper << "	//w_Temper" << endl

			<< N_SubDomains_x << "	//N_SubDomains_x" << endl
			<< N_SubDomains_y << "	//N_SubDomains_y" << endl

			<< ToReadSolvedDataFromSeparateFileForEachSubDomain << "	//ToReadSolvedDataFromSeparateFileForEachSubDomain" << endl
			<< ToWriteSolvedDataToFileForEntireComputationalDomain << "	//ToWriteSolvedDataToFileForEntireComputationalDomain" << endl
			<< ToWriteSolvedDataToSeparateFileForEachSubDomain << "	//ToWriteSolvedDataToSeparateFileForEachSubDomain" << endl
			<< ToOnlyWriteDataFor_Nx_Ny_ForAllSubDomains_AndExit << "	//ToOnlyWriteDataFor_Nx_Ny_ForAllSubDomains_AndExit" << endl
			;

		output_data.close();
	}
	else
	{
		cout << "The program can not open file EnteredData.txt to write data." << endl;
	}


}


void ReadEnteredDataFromFile(void)
{
	using namespace std;

	const int MaxCharRead = 5000;
	char buffer[MaxCharRead];

	ifstream input_data;
	input_data.open("EnteredData.txt");

	if(input_data.is_open())
	{
		input_data >> u_gas_nd;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> v_gas_nd;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> T_gas_nd;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> p_gas_nd;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> Ma;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> Pr;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> gamma1;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> c_mu;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> Kn;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> Fr;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> SolveGravityField;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> Re;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> TypeOfNonDimensiolization;
		input_data.getline(&buffer[0], MaxCharRead, '\n');


		input_data >> Nt;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> Nx;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> Ny;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> N_I;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> N_I_p_c;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> N_I_T;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> ht;
		input_data.getline(&buffer[0], MaxCharRead, '\n');


		input_data >> x_b;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> x_e;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> y_b;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> y_e;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> kind_of_mesh;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> hxmin;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> hxmax;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> xfmin;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> xfmax;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> xbmin;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> xbmax;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> xmid;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> hymin;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> hymax;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ytmin;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ytmax;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ybmin;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ybmax;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ymid;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> MaxError_Velocities;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> MaxError_p_c;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> MaxError_T;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> MaxError_to_stop_program;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> max_time_for_solve;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> solve_is_finished;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> Nt_save_solved_data;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> Nt_save_DragCoefficient;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> Nt_save_temp_data;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> N_I_check;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> WriteDataInShortFormat_kind;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> ToReadSolvedDataFromFile;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ToReadSolvedDataFromBinaryFile;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ToWriteSolvedDataToBinaryFile;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ToWriteSolvedDataToNewFiles;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ToContinueFromInterpolation;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> ToImport_data_for_circles_from_file_b;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ToImport_data_for_polyhedrons_from_file_b;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> ToImport_data_for_circles_from_file_p;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ToImport_data_for_polyhedrons_from_file_p;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> ToImport_data_for_circles_from_file_V;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ToImport_data_for_polyhedrons_from_file_V;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> Given_velocity_on_xb;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> Given_velocity_on_xe;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> Periodic_boundary_conditions_about_OX;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> Pressure_BC;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ToStartFromExactSolutionForChannelFlowWithPressureBC;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> p_BC_xb_correction_method;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> Pressure_ratio_correct;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> GivenReKn;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> correct_p_BC_xb;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> u_given_xb;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> u_given_xb_error;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> p_correction_auto;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> p_correction_min;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> p_correction_max;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> p_correction;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> p_BC_xb;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> T_BC_xb;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> p_BC_xe;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> T_BC_xe;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> dudx_0_BC_xb;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> durhodx_0_BC_xb;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> drhodt_durhodx_dvrhody_0_BC_xb;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ToCalculate_ui_1j_from_equation_for_uij_at_x_b;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> dudx_0_BC_xe;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> durhodx_0_BC_xe;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> drhodt_durhodx_dvrhody_0_BC_xe;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ToCalculate_ui1j_from_equation_for_uij_at_x_e;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> u_uMin_Mout_0_BC_xe;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> dudy_0_BC_yb;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> dudy_0_BC_ye;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> dvdx_0_BC_xb;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ToCalculate_vi_1j_from_equation_for_vij_at_x_b;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> dvdx_0_BC_xe;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ToCalculate_vi1j_from_equation_for_vij_at_x_e;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> dvdy_0_BC_yb;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> dvdy_0_BC_ye;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> dpdx_0_BC_xb;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ToCalculate_pi_1j_from_equation_for_pij_at_x_b;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> dpdx_0_BC_xe;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ToCalculate_pi1j_from_equation_for_pij_at_x_e;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> dpdy_0_BC_yb;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> dpdy_0_BC_ye;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> dTdx_0_BC_xb;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ToCalculate_Temperi_1j_from_equation_for_Temperij_at_x_b;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> dTdx_0_BC_xe;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ToCalculate_Temperi1j_from_equation_for_Temperij_at_x_e;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> dTdy_0_BC_yb;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> dTdy_0_BC_ye;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> To_Use_Kn_local_in_wall_BC;
		input_data.getline(&buffer[0], MaxCharRead, '\n');


		input_data >> ToSolveContinuityEquation;
		input_data.getline(&buffer[0], MaxCharRead, '\n');


		input_data >> use_interpolation_for_u_BC_by_time;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> t_BC_b;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> t_BC_e;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> u_BC_xb_t_BC_b;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> u_BC_xb_t_BC_e;
		input_data.getline(&buffer[0], MaxCharRead, '\n');


		input_data >> w_u;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> w_v;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> w_p;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> w_Temper;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> N_SubDomains_x;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> N_SubDomains_y;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

		input_data >> ToReadSolvedDataFromSeparateFileForEachSubDomain;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ToWriteSolvedDataToFileForEntireComputationalDomain;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ToWriteSolvedDataToSeparateFileForEachSubDomain;
		input_data.getline(&buffer[0], MaxCharRead, '\n');
		input_data >> ToOnlyWriteDataFor_Nx_Ny_ForAllSubDomains_AndExit;
		input_data.getline(&buffer[0], MaxCharRead, '\n');


		input_data.close();
	}
	else
	{
		cout << "The program can not open file File EnteredData.txt to read data"

#ifdef COMPILE_MPI
			<< " process rank = " << rank
			<< " Number_of_processes = " << Number_of_processes
#endif //End of COMPILE_MPI

			<< endl;
	}

}


void WriteSolvedDataToFile(void)
{
	BoundaryConditions_Temperature_TemperatureJumpBC();
	BoundaryConditions_Temperature_TemperatureJumpBC_WriteToMassiveToWriteDataToFile();

#ifdef COMPILE_MPI
	if(ToWriteSolvedDataToSeparateFileForEachSubDomain)
#endif //End of COMPILE_MPI
	{
		if(ToWriteSolvedDataToBinaryFile)
		{
#ifndef COMPILE_MPI
			string u_filename = "u.bin";
			string v_filename = "v.bin";
			string p_filename = "p.bin";
			string Temper_filename = "Temper.bin";
			string vol_inf_fb_filename = "vol_inf_fb.bin";
#endif //End of COMPILE_MPI

#ifdef COMPILE_MPI
			string u_filename = "u." + toString(IJ) + ".bin";
			string v_filename = "v." + toString(IJ) + ".bin";
			string p_filename = "p." + toString(IJ) + ".bin";
			string Temper_filename = "Temper." + toString(IJ) + ".bin";
			string vol_inf_fb_filename = "vol_inf_fb." + toString(IJ) + ".bin";
#endif //End of COMPILE_MPI

			////Write massive to file
			//string v_filename2Dtxt = "v.2D." + toString(IJ) + ".txt";
			//WriteMassiveToFile_2D(v_filename2Dtxt.c_str(), v, Nx, Ny);

			//string vol_inf_fb_filename2Dtxt = "vol_inf_fb.2D." + toString(IJ) + ".txt";
			//WriteMassiveToFile_2D(vol_inf_fb_filename2Dtxt.c_str(), vol_inf_fb, Nx, Ny);


			WriteMassiveToBinaryFile(u_filename.c_str(), u, Na);
			WriteMassiveToBinaryFile(v_filename.c_str(), v, Na);
			WriteMassiveToBinaryFile(p_filename.c_str(), p, Na);
			WriteMassiveToBinaryFile(Temper_filename.c_str(), Temper, Na);
			WriteMassiveToBinaryFile(vol_inf_fb_filename.c_str(), vol_inf_fb, Na);
			//WriteMassiveToBinaryFile("vol_inf_p.bin", vol_inf_p, Na);
			//WriteMassiveToBinaryFile("vol_inf_V.bin", vol_inf_V, Na);

			if(ToWriteSolvedDataToNewFiles)
			{
#ifndef COMPILE_MPI
				string u_filename_new = "u." + toString(it) + ".bin";
				string v_filename_new = "v." + toString(it) + ".bin";
				string p_filename_new = "p." + toString(it) + ".bin";
				string Temper_filename_new = "Temper." + toString(it) + ".bin";
				string vol_inf_fb_filename_new = "vol_inf_fb." + toString(it) + ".bin";
#endif //End of COMPILE_MPI

#ifdef COMPILE_MPI
				string u_filename_new = "u." + toString(IJ) + "." + toString(it) + ".bin";
				string v_filename_new = "v." + toString(IJ) + "." + toString(it) + ".bin";
				string p_filename_new = "p." + toString(IJ) + "." + toString(it) + ".bin";
				string Temper_filename_new = "Temper." + toString(IJ) + "." + toString(it) + ".bin";
				string vol_inf_fb_filename_new = "vol_inf_fb." + toString(IJ) + "." + toString(it) + ".bin";
#endif //End of COMPILE_MPI


				WriteMassiveToBinaryFile(u_filename_new.c_str(), u, Na);
				WriteMassiveToBinaryFile(v_filename_new.c_str(), v, Na);
				WriteMassiveToBinaryFile(p_filename_new.c_str(), p, Na);
				WriteMassiveToBinaryFile(Temper_filename_new.c_str(), Temper, Na);
				WriteMassiveToBinaryFile(vol_inf_fb_filename_new.c_str(), vol_inf_fb, Na);
			}
		}
		else
		{
			//Filenames:
#ifndef COMPILE_MPI
			string u_filename = "u.txt";
			string v_filename = "v.txt";
			string p_filename = "p.txt";
			string Temper_filename = "Temper.txt";
			string vol_inf_fb_filename = "vol_inf_fb.txt";
#endif //End of COMPILE_MPI

#ifdef COMPILE_MPI
			string u_filename = "u." + toString(IJ) + ".txt";
			string v_filename = "v." + toString(IJ) + ".txt";
			string p_filename = "p." + toString(IJ) + ".txt";
			string Temper_filename = "Temper." + toString(IJ) + ".txt";
			string vol_inf_fb_filename = "vol_inf_fb." + toString(IJ) + ".txt";
#endif //End of COMPILE_MPI


			WriteMassiveToFile_1D(u_filename.c_str(), u, Na);
			WriteMassiveToFile_1D(v_filename.c_str(), v, Na);
			WriteMassiveToFile_1D(p_filename.c_str(), p, Na);
			WriteMassiveToFile_1D(Temper_filename.c_str(), Temper, Na);
			WriteMassiveToFile_1D(vol_inf_fb_filename.c_str(), vol_inf_fb, Na);
			//WriteMassiveToFile_1D("vol_inf_p.txt", vol_inf_p, Na);
			//WriteMassiveToFile_1D("vol_inf_V.txt", vol_inf_V, Na);

			if(ToWriteSolvedDataToNewFiles)
			{
#ifndef COMPILE_MPI
				string u_filename_new = "u." + toString(it) + ".txt";
				string v_filename_new = "v." + toString(it) + ".txt";
				string p_filename_new = "p." + toString(it) + ".txt";
				string Temper_filename_new = "Temper." + toString(it) + ".txt";
				string vol_inf_fb_filename_new = "vol_inf_fb.txt." + toString(it) + ".txt";
#endif //End of COMPILE_MPI

#ifdef COMPILE_MPI
				string u_filename_new = "u." + toString(IJ) + "." + toString(it) + ".txt";
				string v_filename_new = "v." + toString(IJ) + "." + toString(it) + ".txt";
				string p_filename_new = "p." + toString(IJ) + "." + toString(it) + ".txt";
				string Temper_filename_new = "Temper." + toString(IJ) + "." + toString(it) + ".txt";
				string vol_inf_fb_filename_new = "vol_inf_fb." + toString(IJ) + "." + toString(it) + ".txt";
#endif //End of COMPILE_MPI


				WriteMassiveToFile_1D(u_filename_new.c_str(), u, Na);
				WriteMassiveToFile_1D(v_filename_new.c_str(), v, Na);
				WriteMassiveToFile_1D(p_filename_new.c_str(), p, Na);
				WriteMassiveToFile_1D(Temper_filename_new.c_str(), Temper, Na);
				WriteMassiveToFile_1D(vol_inf_fb_filename_new.c_str(), vol_inf_fb, Na);
			}
		}
	}

#ifndef COMPILE_MPI
	WriteEnteredDataToFile();
	Write_it_ToFile();
#endif //End of COMPILE_MPI

#ifdef COMPILE_MPI
	if(ToWriteSolvedDataToFileForEntireComputationalDomain)
	{
		MPI_Barrier(MPI_COMM_WORLD);

		////Write data for entire compuattional domain to files
		Write_massive_entire_computational_domain_MPI("u", u, double_type);
		MPI_Barrier(MPI_COMM_WORLD);

		Write_massive_entire_computational_domain_MPI("v", v, double_type);
		MPI_Barrier(MPI_COMM_WORLD);

		Write_massive_entire_computational_domain_MPI("p", p, double_type);
		MPI_Barrier(MPI_COMM_WORLD);

		Write_massive_entire_computational_domain_MPI("Temper", Temper, double_type);
		MPI_Barrier(MPI_COMM_WORLD);

		Write_massive_entire_computational_domain_MPI("vol_inf_fb", vol_inf_fb, unsigned_int_type);
	}

	MPI_Barrier(MPI_COMM_WORLD);


	if(IJ == 0)
	{
		//Write data only from rank 0, otherwise file recovary is possible!!!
		WriteEnteredDataToFile();
		Write_it_ToFile();
	}
	MPI_Barrier(MPI_COMM_WORLD);


#endif //End of COMPILE_MPI

	//Fix the fields before continiue calculation
	//ApplyBodiesConditionsToVariables();
}


void ReadSolvedDataFromFile(void)
{
#ifdef COMPILE_MPI
	if(ToReadSolvedDataFromSeparateFileForEachSubDomain)
#endif //End of COMPILE_MPI
	{
		if(ToReadSolvedDataFromBinaryFile)
		{
			//Filenames:
#ifndef COMPILE_MPI
			string u_filename = "u.bin";
			string v_filename = "v.bin";
			string p_filename = "p.bin";
			string Temper_filename = "Temper.bin";
			string vol_inf_fb_filename = "vol_inf_fb.bin";
#endif //End of COMPILE_MPI

#ifdef COMPILE_MPI
			string u_filename = "u." + toString(IJ) + ".bin";
			string v_filename = "v." + toString(IJ) + ".bin";
			string p_filename = "p." + toString(IJ) + ".bin";
			string Temper_filename = "Temper." + toString(IJ) + ".bin";
			string vol_inf_fb_filename = "vol_inf_fb." + toString(IJ) + ".bin";
#endif //End of COMPILE_MPI


			ReadMassiveFromBinaryFile(u_filename.c_str(), u, Na);
			ReadMassiveFromBinaryFile(v_filename.c_str(), v, Na);
			ReadMassiveFromBinaryFile(p_filename.c_str(), p, Na);
			ReadMassiveFromBinaryFile(Temper_filename.c_str(), Temper, Na);
			ReadMassiveFromBinaryFile(vol_inf_fb_filename.c_str(), vol_inf_fb, Na);
			//ReadMassiveFromBinaryFile("vol_inf_p.bin", vol_inf_p, Na);
			//ReadMassiveFromBinaryFile("vol_inf_V.bin", vol_inf_V, Na);
			//ReadMassiveFromBinaryFile("rho.bin", rho, Na);
		}
		else
		{
			//Filenames:
#ifndef COMPILE_MPI
			string u_filename = "u.txt";
			string v_filename = "v.txt";
			string p_filename = "p.txt";
			string Temper_filename = "Temper.txt";
			string vol_inf_fb_filename = "vol_inf_fb.txt";
#endif //End of COMPILE_MPI

#ifdef COMPILE_MPI
			string u_filename = "u." + toString(IJ) + ".txt";
			string v_filename = "v." + toString(IJ) + ".txt";
			string p_filename = "p." + toString(IJ) + ".txt";
			string Temper_filename = "Temper." + toString(IJ) + ".txt";
			string vol_inf_fb_filename = "vol_inf_fb." + toString(IJ) + ".txt";
#endif //End of COMPILE_MPI


			ReadMassiveFromFile_1D(u_filename.c_str(), u, Na);
			ReadMassiveFromFile_1D(v_filename.c_str(), v, Na);
			ReadMassiveFromFile_1D(p_filename.c_str(), p, Na);
			ReadMassiveFromFile_1D(Temper_filename.c_str(), Temper, Na);
			ReadMassiveFromFile_1D(vol_inf_fb_filename.c_str(), vol_inf_fb, Na);
			//ReadMassiveFromFile_1D("vol_inf_p.txt", vol_inf_p, Na);
			//ReadMassiveFromFile_1D("vol_inf_V.txt", vol_inf_V, Na);
			//ReadMassiveFromFile_1D("rho.txt", rho, Na);

		}

	}
#ifdef COMPILE_MPI
	else
	{
		MPI_Barrier(MPI_COMM_WORLD);

		//Read solved data from file for entire computational domain
		Read_massive_entire_computational_domain_MPI("u", u, double_type);
		MPI_Barrier(MPI_COMM_WORLD);

		Read_massive_entire_computational_domain_MPI("v", v, double_type);
		MPI_Barrier(MPI_COMM_WORLD);

		Read_massive_entire_computational_domain_MPI("p", p, double_type);
		MPI_Barrier(MPI_COMM_WORLD);

		Read_massive_entire_computational_domain_MPI("Temper", Temper, double_type);
		MPI_Barrier(MPI_COMM_WORLD);

		Read_massive_entire_computational_domain_MPI("vol_inf_fb", vol_inf_fb, unsigned_int_type);
		MPI_Barrier(MPI_COMM_WORLD);
	}
#endif //End of COMPILE_MPI


	//Rho will be calculated from p and Temper
	for(counter = 0; counter < Na; counter++)
	{
		if(vol_inf_fb_bool[counter]) rho[counter] = p[counter] / Temper[counter];
		else rho[counter] = 0;
	}
	//Update rho in middle points
	Solve_rho_in_middle_points();

	Read_it_FromFile();
}


void WriteTempDataToFile(void)
{
/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!! IMPORTANT !!!   !!! IMPORTANT !!!   !!! IMPORTANT !!!   !!! IMPORTANT !!!

	DO NOT NOT EXCHANGE DATA WITH MPI_Gather IN THE SAME PROCEDURE WHERE WRITE DATA

	THAT CAN MAY LEAD TO DIFFICULTIES TO OPEN FILE FOR WRITTING!!!

	THIS IS THE CASE WITH Tmp.txt


	CAN NOT OPEN THE SECOND FILE TO WRITE IN THE SAME PROCEDURE - I DO NOT

	KNOW WAY!!!


	!!! IMPORTANT !!!   !!! IMPORTANT !!!   !!! IMPORTANT !!!   !!! IMPORTANT !!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/



#ifdef COMPILE_MPI
	if(IJ == 0) //Only rank 0 write data to file.
#endif //End of COMPILE_MPI
	{
		using namespace std;

		ofstream output_data;


		//If we continue the old solve the file will not be deleted
		if(ToCheckForFileIsNow_Tmp)
		{
			if(ToReadSolvedDataFromFile)
				output_data.open("Tmp.txt", ios::out | ios::app);
			else
				output_data.open("Tmp.txt");
				//output_data.open(FileName_Tmp_tmp.c_str());
		}
		else
			output_data.open("Tmp.txt", ios::out | ios::app);


		if(output_data.is_open())
		{
			output_data
				<< "it = " << it
				<< "	it_counter = " << it_counter
				<< "	Iter = " << Iter
				<< "	Iter_p_c = " << Iter_p_c
				<< "	Iter_T = " << Iter_T
				<< endl;

			unsigned int sum_max_Iter_p_c = 0;

			for(counter = 0; counter < Iter; counter++)
			{
				output_data
					<< "	" << counter
					<< "	" << max_Iter_p_c[counter]

					<< std::setiosflags(std::ios::scientific)
					<< std::setprecision(15)

					<< "	" << max_residual_in_p_Iter_0[counter] << "	" << max_residual_in_p[counter]
					<< "	" << max_residual_in_T_Iter_0[counter] << "	" << max_residual_in_T[counter]
					<< endl;

					sum_max_Iter_p_c += max_Iter_p_c[counter];
			}

			output_data
					<< "		" << sum_max_Iter_p_c
					<< endl;



			output_data
				<< std::setiosflags(std::ios::scientific)
				<< std::setprecision(15)

				<< "	" << max_residual_in_u_Iter_0
				<< "	" << max_residual_in_u << endl

				<< "	" << max_residual_in_v_Iter_0
				<< "	" << max_residual_in_v << endl;


			output_data.close();
		}
		else
		{
			cout << "The program can not open file Tmp.txt to write data." << endl;
		}

	}

	ToCheckForFileIsNow_Tmp = false;

}


void WriteIterationsInLoopsToFile(void)
{
#ifdef COMPILE_MPI
	if(IJ == 0) //Only rank 0 write data to file.
#endif //End of COMPILE_MPI
	{
		using namespace std;


		//Write number of iterations Iter_p_c, Iter_T and Iter in separate file
		ofstream output_data1;

		//If we continue the old solve the file will not be erased
		if(ToCheckForFileIsNow_Iterations_in_loops)
		{
			if(ToReadSolvedDataFromFile)
				output_data1.open("Iterations_in_loops.txt", ios::out | ios::app);
			else
			{
				output_data1.open("Iterations_in_loops.txt");

				if(output_data1.is_open())
				{
					output_data1 << "it"
						<< "	it_counter"
						<< "	Iter"
						<< "	sum_max_Iter_p_c"
						<< endl;
				}

			}
		}
		else
			output_data1.open("Iterations_in_loops.txt", ios::out | ios::app);

		if(output_data1.is_open())
		{
			unsigned int sum_max_Iter_p_c = 0;

			for(counter = 0; counter < Iter; counter++)
			{
				sum_max_Iter_p_c += max_Iter_p_c[counter];
			}

			output_data1
				<< it
				<< "	" << it_counter
				<< "		" << Iter
				<< "	" << sum_max_Iter_p_c
				<< endl;

			output_data1.close();
		}
		else
		{
			cout << "The program can not open file Iterations_in_loops.txt to write data." << endl;
		}

	}

	ToCheckForFileIsNow_Iterations_in_loops = false;
}


void define_area_for_solving(void)
//Define area for solving
{
	Nx_u = new unsigned int [ND];
	Ny_u = new unsigned int [ND];

	Nx_v = new unsigned int [ND];
	Ny_v = new unsigned int [ND];

	Nx_others = new unsigned int [ND];
	Ny_others = new unsigned int [ND];


	//Here the Subdomains are already defined, if they exist.
	//If in solving problem are subdomains all Subdomains boundaries have to be
	//Nx_u[0] = 1; Nx_u[0] = Nx - 1; Ny_u[0] = 1; Ny_u[0] = Ny - 1;
	//Nx_v[0] = 1; Nx_v[1] = Nx - 1; Ny_v[0] = 1; Ny_v[0] = Ny - 1;
	//Nx_others[0] = 1; Nx_others[1] = Nx - 1; Ny_others[0] = 1; Ny_others[0] = Ny - 1;
	//
	//If the solving problem is with Periodic_boundary_conditions_about_OX and for I == 0,
	//because of the swap area, which is needed.
	//Nx_u[0] = 1; Ny_u[0] = 1;
	//Nx_v[0] = 1; Ny_v[0] = 1;
	//Nx_others[0] = 1; Ny_others[0] = 1;
	//
	//If the solving problem is with Periodic_boundary_conditions_about_OX and for I == N_SubDomain_x,
	//because of the swap area, which is needed.
	//Nx_u[0] = Nx - 1;Ny_u[0] = Ny - 1;
	//Nx_v[1] = Nx - 1;Ny_v[0] = Ny - 1;
	//Nx_others[1] = Nx - 1;Ny_others[0] = Ny - 1;
	//
	//


	//The boundary on x_b
#ifdef COMPILE_MPI
	if(SubDomain[IJ].is_SubDomain_in_direction_I_1)
	{
		Nx_u[0] = BlockLengthSwap_I_1;
		Nx_v[0] = BlockLengthSwap_I_1;
		Nx_others[0] = BlockLengthSwap_I_1;
	}
	else
#endif //End of #ifdef COMPILE_MPI
	{
		//There is no other subdomanin in x - 1 direction => this is first column
		if(Periodic_boundary_conditions_about_OX)
		{
			Nx_u[0] = 0;
			Nx_v[0] = 0;
			Nx_others[0] = 0;
		}
		else if(Pressure_BC)
		{
			Nx_u[0] = 1;
			Nx_v[0] = 1;
			Nx_others[0] = 1;
		}
		else
		{
			if(Given_velocity_on_xb)
			{
				Nx_u[0] = 2;
				Nx_v[0] = 1;
				Nx_others[0] = 1;
			}
			else
			{
				Nx_u[0] = 1;
				Nx_v[0] = 1;
				Nx_others[0] = 1;
			}
		}

	}


	//The boundary on x_e
#ifdef COMPILE_MPI
	if(SubDomain[IJ].is_SubDomain_in_direction_I1)
	{
		Nx_u[1] = Nx - BlockLengthSwap_I1;
		Nx_v[1] = Nx - BlockLengthSwap_I1;
		Nx_others[1] = Nx - BlockLengthSwap_I1;
	}
	else
#endif //End of #ifdef COMPILE_MPI
	{
		//There is no other subdomanin in x + 1 direction => this is last column
		if(Periodic_boundary_conditions_about_OX)
		{
			Nx_u[1] = Nx;// - 1;
			Nx_v[1] = Nx;// - 1;
			Nx_others[1] = Nx;// - 1;
		}
		else if(Pressure_BC)
		{
			Nx_u[1] = Nx - 1;
			Nx_v[1] = Nx - 2;
			Nx_others[1] = Nx - 2;
		}
		else
		{
			if(Given_velocity_on_xe)
			{
				Nx_u[1] = Nx - 1;
				Nx_v[1] = Nx - 2;
				Nx_others[1] = Nx - 1;
			}
			else
			{
				Nx_u[1] = Nx - 1;
				Nx_v[1] = Nx - 2;
				Nx_others[1] = Nx - 2;
			}
		}
	}


	//The boundary on y_b
#ifdef COMPILE_MPI
	if(SubDomain[IJ].is_SubDomain_in_direction_J_1)
	{
		Ny_u[0] = BlockLengthSwap_J_1;
		Ny_v[0] = BlockLengthSwap_J_1;
		Ny_others[0] = BlockLengthSwap_J_1;
	}
	else
#endif //End of #ifdef COMPILE_MPI
	{
		//There is no other subdomanin in y - 1 direction => this is first raw
		if(Periodic_boundary_conditions_about_OX)
		{
			Ny_u[0] = 1;
			Ny_v[0] = 1;
			Ny_others[0] = 1;
		}
		else if(Pressure_BC)
		{
			Ny_u[0] = 1;
			Ny_v[0] = 1;
			Ny_others[0] = 1;
		}
		else
		{
			Ny_u[0] = 1;
			Ny_v[1] = Ny - 1;
			Ny_others[0] = 1;
		}
	}


	//The boundary on y_e
#ifdef COMPILE_MPI
	if(SubDomain[IJ].is_SubDomain_in_direction_J1)
	{
		Ny_u[1] = Ny - BlockLengthSwap_J1;
		Ny_v[1] = Ny - BlockLengthSwap_J1;
		Ny_others[1] = Ny - BlockLengthSwap_J1;
	}
	else
#endif //End of #ifdef COMPILE_MPI
	{
		//There is no other subdomanin in y - 1 direction => this is last raw
		if(Periodic_boundary_conditions_about_OX)
		{
			Ny_u[1] = Ny - 2;
			Ny_v[1] = Ny - 1;
			Ny_others[1] = Ny - 2;
		}
		else if(Pressure_BC)
		{
			Ny_u[1] = Ny - 2;
			Ny_v[1] = Ny - 1;
			Ny_others[1] = Ny - 2;
		}
		else
		{
			Ny_u[1] = Ny - 2;
			Ny_v[0] = 1;
			Ny_others[1] = Ny - 2;
		}
	}


//OLD
//	if(Periodic_boundary_conditions_about_OX)
//	{
//		Nx_u[0] = 0 + i_BLswap_I_1;
//		Nx_u[1] = Nx - i_swap_I1;// - 1;
//		Ny_u[0] = 1 + j_swap_J_1;
//		Ny_u[1] = Ny - 2 + j_swap_J1;
//
//		Nx_v[0] = 0 + i_BLswap_I_1;
//		Nx_v[1] = Nx - i_swap_I1;// - 1;
//		Ny_v[0] = 1 + j_swap_J_1;
//		Ny_v[1] = Ny - 1;
//
//		Nx_others[0] = 0 + i_BLswap_I_1;
//		Nx_others[1] = Nx - i_swap_I1;// - 1;
//		Ny_others[0] = 1 + j_swap_J_1;
//		Ny_others[1] = Ny - 2 + j_swap_J1;
//	}
//	else if(Pressure_BC)
//	{
//		Nx_u[0] = 1 + i_swap_I_1;
//		Nx_u[1] = Nx - 1;
//		Ny_u[0] = 1 + j_swap_J_1;
//		Ny_u[1] = Ny - 2 + j_swap_J1;
//
//		Nx_v[0] = 1 + i_swap_I_1;
//		Nx_v[1] = Nx - 2 + i_swap_I1;
//		Ny_v[0] = 1 + j_swap_J_1;
//		Ny_v[1] = Ny - 1;
//
//		Nx_others[0] = 1 + i_swap_I_1;
//		Nx_others[1] = Nx - 2 + i_swap_I1;
//		Ny_others[0] = 1 + j_swap_J_1;
//		Ny_others[1] = Ny - 2 + j_swap_J1;
//	}
//	else
//	{
//		if(Given_velocity_on_xb)
//		{
//			Nx_u[0] = 2 - i_swap_I_1;
//			Nx_v[0] = 1 + i_swap_I_1;
//			Nx_others[0] = 1 + i_swap_I_1;
//		}
//		else
//		{
//			Nx_u[0] = 1 + i_swap_I_1;
//			Nx_v[0] = 1 + i_swap_I_1;
//			Nx_others[0] = 1 + i_swap_I_1;
//		}
//
//		if(Given_velocity_on_xe)
//		{
//			Nx_u[1] = Nx - 1;
//			Nx_v[1] = Nx - 2 + i_swap_I1;
//			Nx_others[1] = Nx - 1;
//		}
//		else
//		{
//			Nx_u[1] = Nx - 1;
//			Nx_v[1] = Nx - 2 + i_swap_I1;
//			Nx_others[1] = Nx - 2 + i_swap_I1;
//		}
//
//		Ny_u[0] = 1 + j_swap_J_1;
//		Ny_u[1] = Ny - 2 + j_swap_J1;
//
//		Ny_v[0] = 1 + j_swap_J_1;
//		Ny_v[1] = Ny - 1;
//
//		Ny_others[0] = 1 + j_swap_J_1;
//		Ny_others[1] = Ny - 2 + j_swap_J1;
//	}

}

void define_OforS(void)
{
	N_vol_inf = 3;

	OforS = new Objects_for_Solving[N_vol_inf];

    char file_name_solve_polyhedrons_b [] = "solve_polyhedrons_b.txt";
    char file_name_solve_polyhedrons_p [] = "solve_polyhedrons_p.txt";
    char file_name_solve_polyhedrons_V [] = "solve_polyhedrons_V.txt";


	//Import data for given body
	if(ToImport_data_for_circles_from_file_b)
		OforS[gb].Import_data_for_circles_from_file("solve_circles_b.txt");

	if(ToImport_data_for_polyhedrons_from_file_b)
		OforS[gb].Import_data_for_polyhedrons_from_file(file_name_solve_polyhedrons_b, gb);


	//Import data for given pressure
	if(ToImport_data_for_circles_from_file_p)
		OforS[gp].Import_data_for_circles_from_file("solve_circles_p.txt");

	if(ToImport_data_for_polyhedrons_from_file_p)
		OforS[gp].Import_data_for_polyhedrons_from_file(file_name_solve_polyhedrons_p, gp);


	//Import data for given Velocity
	if(ToImport_data_for_circles_from_file_V)
		OforS[gV].Import_data_for_circles_from_file("solve_circles_V.txt");

	if(ToImport_data_for_polyhedrons_from_file_V)
		OforS[gV].Import_data_for_polyhedrons_from_file(file_name_solve_polyhedrons_V, gV);

}


void define_vol_inf(void)
{
	//Import data for solving bodies and areas
	define_OforS();


	//DEFINE vol_inf_fb
	if(ToImport_data_for_circles_from_file_b || ToImport_data_for_polyhedrons_from_file_b)
	{
		vol_inf_fb = new unsigned int[Na];

		//To solve rule matrixes vol_inf
		for(counter = 0; counter < Na; counter++)
		{
			vol_inf_fb[counter] = fluid;
		}

		//Number of all objects for solving body
		unsigned int N_all_objects_b = 0;
		N_all_objects_b = OforS[gb].N_circles + OforS[gb].N_polyhedrons;

		//Define body in solved area, if we want to solve it
		if(N_all_objects_b > 0)
		{
			//Verifing circls
			for(counter = 0; counter < OforS[gb].N_circles; counter++)
			{
				for (j = 0; j < Ny; j++)
				{
					for (i = 0; i < Nx; i++)
					{
						if(OforS[gb].M_circles[counter].is_point_inside_circle(x_v[i], y_v[j]))
							vol_inf_fb[i + j * Nx] = Ncp_beg_b + counter;
					}
				}
			}


			//Verifing polyhedrond
			for(counter = 0; counter < OforS[gb].N_polyhedrons; counter++)
			{
				for (j = 0; j < Ny; j++)
				{
					for (i = 0; i < Nx; i++)
					{
						if(OforS[gb].M_polyhedrons[counter].is_point_inside_polyhedron(x_v[i], y_v[j]))
							vol_inf_fb[i + j * Nx] = Ncp_beg_b + OforS[gb].N_circles + counter;
					}
				}
			}


		}

	}


	//DEFINE vol_inf_p
	if(ToImport_data_for_circles_from_file_p || ToImport_data_for_polyhedrons_from_file_p)
	{
		vol_inf_p = new unsigned int[Na];

		//To solve rule matrixes vol_inf
		for(counter = 0; counter < Na; counter++)
		{
			vol_inf_p[counter] = no_gp;
		}

		//Number of all objects for solving pressure
		unsigned int N_all_objects_p = 0;
		N_all_objects_p = OforS[gp].N_circles + OforS[gp].N_polyhedrons;

		//Define body in solved area, if we want to solve it
		if(N_all_objects_p > 0)
		{
			//Verifing circls
			for(counter = 0; counter < OforS[gp].N_circles; counter++)
			{
				for (j = 0; j < Ny; j++)
				{
					for (i = 0; i < Nx; i++)
					{
						if(OforS[gp].M_circles[counter].is_point_inside_circle(x_v[i], y_v[j]))
							vol_inf_p[i + j * Nx] = Ncp_beg_p + counter;
					}
				}
			}


			//Verifing polyhedrond
			for(counter = 0; counter < OforS[gp].N_polyhedrons; counter++)
			{
				for (j = 0; j < Ny; j++)
				{
					for (i = 0; i < Nx; i++)
					{
						if(OforS[gp].M_polyhedrons[counter].is_point_inside_polyhedron(x_v[i], y_v[j]))
							vol_inf_p[i + j * Nx] = Ncp_beg_p + OforS[gp].N_circles + counter;
					}
				}
			}


		}

	}


	//DEFINE vol_inf_V
	if(ToImport_data_for_circles_from_file_V || ToImport_data_for_polyhedrons_from_file_V)
	{
		vol_inf_V = new unsigned int[Na];

		//To solve rule matrixes vol_inf
		for(counter = 0; counter < Na; counter++)
		{
			vol_inf_V[counter] = no_gV;
		}

		//Number of all objects for solving Velocity
		unsigned int N_all_objects_V = 0;
		N_all_objects_V = OforS[gV].N_circles + OforS[gV].N_polyhedrons;

		//Define body in solved area, if we want to solve it
		if(N_all_objects_V > 0)
		{
			//Verifing circls
			for(counter = 0; counter < OforS[gV].N_circles; counter++)
			{
				for (j = 0; j < Ny; j++)
				{
					for (i = 0; i < Nx; i++)
					{
						if(OforS[gV].M_circles[counter].is_point_inside_circle(x_v[i], y_v[j]))
							vol_inf_V[i + j * Nx] = Ncp_beg_V + counter;
					}
				}
			}


			//Verifing polyhedrond
			for(counter = 0; counter < OforS[gV].N_polyhedrons; counter++)
			{
				for (j = 0; j < Ny; j++)
				{
					for (i = 0; i < Nx; i++)
					{
						if(OforS[gV].M_polyhedrons[counter].is_point_inside_polyhedron(x_v[i], y_v[j]))
							vol_inf_V[i + j * Nx] = Ncp_beg_V + OforS[gV].N_circles + counter;
					}
				}
			}


		}

	}



}


//Apply bodies velocities to massives for u, v, T, rho and p
void ApplyBodiesConditionsToVariables(void)
{
	//conditions for u
	for (i = 0; i < Nx; i++)
	{
		for (j = 0; j < Ny; j++)
		{

			node = i + j * Nx;
			ij_function_i_and_j_to_the_boundaries(node);

			if (Ncp_beg_b <= vol_inf_fb[ij])
			{
				//Verifing is point in circle or in polyhedron
				if(vol_inf_fb[ij] < (Ncp_beg_b + OforS[gb].N_circles))
				{
					//solving in which circle is point
					counter = vol_inf_fb[ij] - Ncp_beg_b;

					u[ij] = OforS[gb].M_circles[counter].u_body;
					if(i < Nx) u[i1j] = OforS[gb].M_circles[counter].u_body;
				}
				else
				{
					//solving in which polyhedron is point
					counter = vol_inf_fb[ij] - (Ncp_beg_b + OforS[gb].N_circles);

					u[ij] = OforS[gb].M_polyhedrons[counter].u_body;
					if(i < Nx) u[i1j] = OforS[gb].M_polyhedrons[counter].u_body;
				}

			}

		}

	}



	//conditions for v
	for (i = 0; i < Nx; i++)
	{
		for (j = 0; j < Ny; j++)
		{

			node = i + j * Nx;
			ij_function_i_and_j_to_the_boundaries(node);

			if (Ncp_beg_b <= vol_inf_fb[ij])
			{
				//Verifing is point in circle or in polyhedron
				if(vol_inf_fb[ij] < (Ncp_beg_b + OforS[gb].N_circles))
				{
					//solving in which circle is point
					counter = vol_inf_fb[ij] - Ncp_beg_b;

					v[ij] = OforS[gb].M_circles[counter].v_body;
					if(j < Ny) v[ij1] = OforS[gb].M_circles[counter].v_body;
				}
				else
				{
					//solving in which polyhedron is point
					counter = vol_inf_fb[ij] - (Ncp_beg_b + OforS[gb].N_circles);

					v[ij] = OforS[gb].M_polyhedrons[counter].v_body;
					if(j < Ny) v[ij1] = OforS[gb].M_polyhedrons[counter].v_body;
				}

			}

		}

	}


	//conditions for p_pr, rho and Temper in body
	for (i = 0; i < Nx; i++)
	{
		for (j = 0; j < Ny; j++)
		{
			node = i + j * Nx;

			if (Ncp_beg_b <= vol_inf_fb[node])
			{
				//Verifing is point in circle or in polyhedron
				if(vol_inf_fb[node] < (Ncp_beg_b + OforS[gb].N_circles))
				{
					//solving in which circle is point
					counter = vol_inf_fb[node] - Ncp_beg_b;

					p[node] = OforS[gb].M_circles[counter].p_body;
					Temper[node] = OforS[gb].M_circles[counter].T_body;

					if(vol_inf_fb_bool[node]) rho[node] = p[node] / Temper[node];
					else rho[node] = 0;
				}
				else
				{
					//solving in which polyhedron is point
					counter = vol_inf_fb[node] - (Ncp_beg_b + OforS[gb].N_circles);

					p[node] = OforS[gb].M_polyhedrons[counter].p_body;
					Temper[node] = OforS[gb].M_polyhedrons[counter].T_body;

					if(vol_inf_fb_bool[node]) rho[node] = p[node] / Temper[node];
					else rho[node] = 0;
				}

			}

		}

	}

}

void ApplyPressureConditionsToPressure(void)
{
	//conditions for p_pr and rho in area for given pressure
	for (i = 0; i < Nx; i++)
	{
		for (j = 0; j < Ny; j++)
		{

			node = i + j * Nx;
			ij_function(node);

			if (Ncp_beg_p <= vol_inf_p[ij])
			{
				//Verifing is point in circle or in polyhedron
				if(vol_inf_p[ij] < (Ncp_beg_p + OforS[gp].N_circles))
				{
					//solving in which circle is point
					counter = vol_inf_p[ij] - Ncp_beg_p;

					p[node] = OforS[gp].M_circles[counter].p_body;

					if(vol_inf_fb_bool[node]) rho[node] = OforS[gp].M_circles[counter].p_body / OforS[gp].M_circles[counter].T_body;
					else rho[node] = 0;

				}
				else
				{
					//solving in which polyhedron is point
					counter = vol_inf_p[ij] - (Ncp_beg_p + OforS[gp].N_circles);

					p[node] = OforS[gp].M_polyhedrons[counter].p_body;

					if(vol_inf_fb_bool[node]) rho[node] = OforS[gp].M_polyhedrons[counter].p_body / OforS[gp].M_polyhedrons[counter].T_body;
					else rho[node] = 0;
				}

			}

		}

	}

}

//Apply bodies velocities to massives for u, v, T, rho and p
void ApplyVelocityConditionsToVelocities(void)
{
	//conditions for u
	for (i = 0; i < Nx - 1; i++)
	{
		for (j = 0; j < Ny; j++)
		{

			node = i + j * Nx;
			ij_function(node);

			if (Ncp_beg_V <= vol_inf_V[ij])
			{
				//Verifing is point in circle or in polyhedron
				if(vol_inf_V[ij] < (Ncp_beg_b + OforS[gV].N_circles))
				{
					//solving in which circle is point
					counter = vol_inf_fb[ij] - Ncp_beg_V;

					u[ij] = OforS[gV].M_circles[counter].u_body;
					u[i1j] = OforS[gV].M_circles[counter].u_body;
				}
				else
				{
					//solving in which polyhedron is point
					counter = vol_inf_V[ij] - (Ncp_beg_b + OforS[gV].N_circles);

					u[ij] = OforS[gV].M_polyhedrons[counter].u_body;
					u[i1j] = OforS[gV].M_polyhedrons[counter].u_body;
				}

			}

		}

	}



	//conditions for v
	for (i = 0; i < Nx; i++)
	{
		for (j = 0; j < Ny - 1; j++)
		{

			node = i + j * Nx;
			ij_function(node);

			if (Ncp_beg_V <= vol_inf_V[ij])
			{
				//Verifing is point in circle or in polyhedron
				if(vol_inf_V[ij] < (Ncp_beg_b + OforS[gV].N_circles))
				{
					//solving in which circle is point
					counter = vol_inf_fb[ij] - Ncp_beg_V;

					v[ij] = OforS[gV].M_circles[counter].v_body;
					v[ij1] = OforS[gV].M_circles[counter].v_body;
				}
				else
				{
					//solving in which polyhedron is point
					counter = vol_inf_V[ij] - (Ncp_beg_V + OforS[gV].N_circles);

					v[ij] = OforS[gV].M_polyhedrons[counter].v_body;
					v[ij1] = OforS[gV].M_polyhedrons[counter].v_body;
				}
			}
		}

	}

	//conditions for p, rho and T
	for (i = 0; i < Nx; i++)
	{
		for (j = 0; j < Ny - 1; j++)
		{

			node = i + j * Nx;
			ij_function(node);

			if (Ncp_beg_V <= vol_inf_V[ij])
			{
				//Verifing is point in circle or in polyhedron
				if(vol_inf_V[ij] < (Ncp_beg_b + OforS[gV].N_circles))
				{
					//solving in which circle is point
					counter = vol_inf_fb[ij] - Ncp_beg_V;

					p[ij] = OforS[gV].M_circles[counter].p_body;
					Temper[ij] = OforS[gV].M_circles[counter].T_body;

				}
				else
				{
					//solving in which polyhedron is point
					counter = vol_inf_V[ij] - (Ncp_beg_V + OforS[gV].N_circles);

					p[ij] = OforS[gV].M_polyhedrons[counter].p_body;
					Temper[ij] = OforS[gV].M_polyhedrons[counter].T_body;

					if(vol_inf_fb_bool[node]) rho[node] = p[node] / Temper[node];
					else rho[node] = 0;
				}
			}
		}

	}

}

void Write_it_ToFile(void)
{
	using namespace std;

	ofstream output_data;
	output_data.open("it_begin.txt");

	if(output_data.is_open())
	{
		output_data << it << endl;

		output_data.close();
	}
	else
	{
		cout << "The program can not open file it_begin.txt to write data." << endl;
	}

}

void Read_it_FromFile(void)
{
	using namespace std;

	ifstream input_data;
	input_data.open("it_begin.txt");

	if(input_data.is_open())
	{
		input_data >> it_begin;
		it = it_begin;

		input_data.close();
	}
	else
	{
		cout << "File it_begin.txt, can not be opened to read data."

#ifdef COMPILE_MPI
			<< " process rank = " << rank
			<< " Number_of_processes = " << Number_of_processes
#endif //End of COMPILE_MPI

			<< endl;
	}
}


inline void BoundaryConditions_Velocities(void)
/*
Target:	For this element are applyed Bounary Conditions For Velocities In Solved Area.
*/
{
	if(!Periodic_boundary_conditions_about_OX)
	{
		//Before solve BoundaryConditions_Velocities for u we
		//MUST be sure that rho_u and rho_v are calculated - Solve_rho_in_middle_points()

		//First must solve v, becouse is used when we solve u,
		//and rho, becouse is used when we solve u in some cases
		//dvdx_0_BC_xb give dv / dx = 0, on x = x_b
		if(dvdx_0_BC_xb
#ifdef COMPILE_MPI
				&& (!SubDomain[IJ].is_SubDomain_in_direction_I_1)
#endif //End of COMPILE_MPI
			)
		{
			for (j = 0; j < Ny; j++)
			{
				node = Nx_v[0] + j * Nx;
				for (i = 0; i < Nx_v[0]; i++)
				{
					if(vol_inf_fb_bool[node])
						v[i + j * Nx] = v[node];
				}
			}
		}

		//dvdx_0_BC_xe give dv / dx = 0, on x = x_e
		if(dvdx_0_BC_xe
#ifdef COMPILE_MPI
				&& (!SubDomain[IJ].is_SubDomain_in_direction_I1)
#endif //End of COMPILE_MPI
			)
		{
			for (j = 0; j < Ny; j++)
			{
				node = Nx_v[1] - 1 + j * Nx;
				for (i = Nx_v[1]; i < Nx; i++)
				{
					if(vol_inf_fb_bool[node])
					{
						//dvdx=0;
						v[i + j * Nx] = v[node];
					}
				}
			}
		}


		//dvdy_0_BC_xb give: dv/dy = 0, on y = y_b
		if(dvdy_0_BC_yb
#ifdef COMPILE_MPI
				&& (!SubDomain[IJ].is_SubDomain_in_direction_J_1)
#endif //End of COMPILE_MPI
			)
		{
			for(i = 0; i < Nx; i++)
			{
				node = i + Ny_v[0] * Nx;
				for(j = 0; j < Ny_v[0]; j++)
				{
					if(vol_inf_fb_bool[node])
						v[i + j * Nx] = v[node];
				}
			}
		}

		//dvdy_0_BC_xe give: dv/dy = 0, on y = y_e
		if(dvdy_0_BC_ye
#ifdef COMPILE_MPI
			&& (!SubDomain[IJ].is_SubDomain_in_direction_J1)
#endif //End of COMPILE_MPI
			)
		{
			for(i = 0; i < Nx; i++)
			{
				node = i + (Ny_v[1] - 1) * Nx;
				for(j = Ny_v[1]; j < Ny; j++)
				{
					if(vol_inf_fb_bool[node])
						v[i + j * Nx] = v[node];
				}
			}

		}


		//boundary conditions for u
		//dudx_0_BC_xb give du / dx = 0, on x = x_b
		if(dudx_0_BC_xb
#ifdef COMPILE_MPI
			&& (!SubDomain[IJ].is_SubDomain_in_direction_I_1)
#endif //End of COMPILE_MPI
			)
		{
			for (j = 0; j < Ny; j++)
			{
				//boundary conditions for u in xb
				node = Nx_u[0] + j * Nx;
				for (i = 0; i < Nx_u[0]; i++)
				{
					if(vol_inf_fb_bool[node])
						u[i + j * Nx] = u[node];
						//u[1 + j * Nx] = u[node];
				}
			}

		}

		//durhodx_0_BC_xb - that give: d(u * rho) / dx = 0, on x = x_b;
		if(durhodx_0_BC_xb
#ifdef COMPILE_MPI
			&& (!SubDomain[IJ].is_SubDomain_in_direction_I_1)
#endif //End of COMPILE_MPI
			)
		{
			for (j = 0; j < Ny; j++)
			{
				//boundary conditions for u in xe
				node = Nx_u[0] - 1 + j * Nx;

				if(vol_inf_fb_bool[node])
				{
					u[node] = (rho_u[node + 1] * u[node + 1]) / rho_u[node];
				}

			}
		}


		//drhodt_durhodx_dvrhody_0_BC_xb - that give: drho/dt + d(u * rho)/dx + d(v * rho)/dy = 0; - continuity equation over the x_b boundary
		if(drhodt_durhodx_dvrhody_0_BC_xb
#ifdef COMPILE_MPI
			&& (!SubDomain[IJ].is_SubDomain_in_direction_I_1)
#endif //End of COMPILE_MPI
			)
		{
			i = Nx_u[0] - 1;
			for (j = 0; j < Ny; j++)
			{
				//boundary conditions for u in xe
				node = i + j * Nx;

				if(vol_inf_fb_bool[node])
				{
					u[node] =
						(
						(rho[node] - rho_pr[node]) * hx[i] * hy[j] / ht
						+ (rho_v[node + Nx] * v[node + Nx] - rho_v[node] * v[node]) * hx[i]
						+ rho_u[node + 1] * u[node + 1] * hy[j]
						) / (rho_u[node] * hy[j]);
				}

			}

		}

		//dudx_0_BC_xe give du / dx = 0, on x = x_e
		if(dudx_0_BC_xe
#ifdef COMPILE_MPI
				&& (!SubDomain[IJ].is_SubDomain_in_direction_I1)
#endif //End of COMPILE_MPI
			)
		{
			for (j = 0; j < Ny; j++)
			{
				//boundary conditions for u in xe
				node = Nx_u[1] - 1 + j * Nx;
				for (i = Nx_u[1]; i < Nx; i++)
				{
					if(vol_inf_fb_bool[node])
					{
						//dudx=0;
						u[i + j * Nx] = u[node];
					}
				}
			}
		}

		//durhodx_0_BC_xe - that give: d(u * rho) / dx = 0, on x = x_b;
		if(durhodx_0_BC_xe
#ifdef COMPILE_MPI
			&& (!SubDomain[IJ].is_SubDomain_in_direction_I1)
#endif //End of COMPILE_MPI
			)
		{
			for (j = 0; j < Ny; j++)
			{
				//boundary conditions for u in xe
				node = Nx_u[1] - 1 + j * Nx;

				if(vol_inf_fb_bool[node])
				{
					u[node + 1] = (rho_u[node] * u[node]) / rho_u[node + 1];
				}

			}
		}

		//drhodt_durhodx_dvrhody_0_BC_xe - that give: drho/dt + d(u * rho)/dx + d(v * rho)/dy = 0; - continuity equation over the x_e boundary
		if(drhodt_durhodx_dvrhody_0_BC_xe
#ifdef COMPILE_MPI
			&& (!SubDomain[IJ].is_SubDomain_in_direction_I1)
#endif //End of COMPILE_MPI
			)
		{
			i = Nx_u[1] - 1;
			for (j = 0; j < Ny; j++)
			{
				//boundary conditions for u in xe
				node = i + j * Nx;

				if(vol_inf_fb_bool[node])
				{
					u[node + 1] =
						(
						-(rho[node] - rho_pr[node]) * hx[i] * hy[j] / ht
						- (rho_v[node + Nx] * v[node + Nx] - rho_v[node] * v[node]) * hx[i]
						+ rho_u[node] * u[node] * hy[j]
						) / (rho_u[node + 1] * hy[j]);
				}

			}
		}

#ifndef COMPILE_MPI
		//u_uMin_Mout_0_BC_xe give: u[i1j] - u[ij] * Min / Mout = 0, on x = x_e
		//--> u[i1j] - u[ij] * Min / Mout, on x = x_e
		//where:
		//	Min - mass flux coming into domain
		//	Mout - mass flux leaveing domain
		if(u_uMin_Mout_0_BC_xe)
		{
			//solve Min
			Min = 0;
			i = Nx_u[0];
			for(j = 0; j < Ny; j++)
			{
				//Min = Fx[i + j * Nx];
				Min = rho_u[i + j * Nx];
			}

			//solve Mout
			Mout = 0;
			i = Nx_u[1];
			for(j = 0; j < Ny; j++)
			{
				//Mout = Fx[i + j * Nx];
				Mout = rho_u[i + j * Nx];
			}


			//solve u out
			i = Nx_u[1] - 1;
			for (j = 0; j < Ny; j++)
			{
				//boundary conditions for u in xe
				node = Nx_u[1] - 1 + j * Nx;

				u[node + 1] = u[node] * Min / Mout;
			}
		}
#endif //End of COMPILE_MPI


		//dudy_0_BC_yb give: du/dy = 0, on y = y_b
		if(dudy_0_BC_yb
#ifdef COMPILE_MPI
			&& (!SubDomain[IJ].is_SubDomain_in_direction_J_1)
#endif //End of COMPILE_MPI
			)
		{
			for(i = 0; i < Nx; i++)
			{
				node = i + Ny_u[0] * Nx;
				for(j = 0; j < Ny_u[0]; j++)
				{
					if(vol_inf_fb_bool[node])
						u[i + j * Nx] = u[node];
				}
			}
		}

		//dudy_0_BC_ye give: du/dy = 0, on y = y_e
		if(dudy_0_BC_ye
#ifdef COMPILE_MPI
			&& (!SubDomain[IJ].is_SubDomain_in_direction_J1)
#endif //End of COMPILE_MPI
			)
		{
			for(i = 0; i < Nx; i++)
			{
				node = i + (Ny_u[1] - 1) * Nx;
				for(j = Ny_u[1]; j < Ny; j++)
				{
					if(vol_inf_fb_bool[node])
						u[i + j * Nx] = u[node];
				}
			}
		}

	}

}

inline void BoundaryConditions_Velocities_SlipBC_u(void)
{
	//Velocity Slip Boundary Conditions for u
	//Here have to be checked all calculated nodes in this SubDomain.
	//Because of the case, when we have channel wall only in i == 0,
	//we have to check this node. For the case, when we have SubDomains on
	//surface boundary this is also good, because the internal poits of SubDomain
	//are calculated using point from other SubDomain.
	for(counter = 0; counter < Na_BoundaryConditions_Velocities_SlipBC_u; counter++)
	{
		node = rule_vector_ij_BoundaryConditions_Velocities_SlipBC_u[counter];
		j = rule_vector_j_BoundaryConditions_Velocities_SlipBC_u[counter];
		i = node - j * Nx;

		if(Periodic_boundary_conditions_about_OX)
		{
			ij_function_PeriodicBC_i_and_j_to_the_boundaries(i, j);
		}
		else
		{
			ij_function_i_and_j_to_the_boundaries(node);
		}

		is_VelocitySlipBC_body_tmp = false;

		//Now check the kind of body: circle or polyhedron
		if((vol_inf_fb[ij] >= Ncp_beg_b)
			&& (vol_inf_fb[ij] < (Ncp_beg_b + OforS[gb].N_circles)))
		{
			//The body is circle
			//Velocity Slip Boundary Conditions for circles.
			//The number of circle in this node is:
			N_circle = vol_inf_fb[ij] - Ncp_beg_b;

			u_body_tmp = OforS[gb].M_circles[N_circle].u_body;

			is_VelocitySlipBC_body_tmp = OforS[gb].M_circles[N_circle].is_VelocitySlipBC_body;
			w_VelocitySlipBC_tmp = OforS[gb].M_circles[N_circle].w_VelocitySlipBC;
			F_VelocitySlip_tmp = OforS[gb].M_circles[N_circle].F_VelocitySlip;
		}
		else if(vol_inf_fb[ij] >= (Ncp_beg_b + OforS[gb].N_circles))
		{
			//Velocity Slip Boundary Conditions for polyhedrons.
			//The number of polyhedron in this node is:
			N_polyhedron = vol_inf_fb[ij] - (Ncp_beg_b + OforS[gb].N_circles);

			u_body_tmp = OforS[gb].M_polyhedrons[N_polyhedron].u_body;

			is_VelocitySlipBC_body_tmp = OforS[gb].M_polyhedrons[N_polyhedron].is_VelocitySlipBC_body;
			w_VelocitySlipBC_tmp = OforS[gb].M_polyhedrons[N_polyhedron].w_VelocitySlipBC;
			F_VelocitySlip_tmp = OforS[gb].M_polyhedrons[N_polyhedron].F_VelocitySlip;
		}


		//Apply velocity slip BC
		//This condition here is enough, because if there is no body to this volume
		//is_VelocitySlipBC_body_tmp, will left false.
		if(is_VelocitySlipBC_body_tmp)
		{
			if(fluid == vol_inf_fb[ij_1])
			{
				//In this case fluid is bottom from the wall
				if(To_Use_Kn_local_in_wall_BC)
				{
					//Use Kn local
					a_tmp_ij = (F_VelocitySlip_tmp * Kn) / (rho_u[ij_1] * 0.5 * hy[j - 1 * (0 < j)]);
					a_tmp_i1j = (F_VelocitySlip_tmp * Kn) / (rho_u[i1j_1] * 0.5 * hy[j - 1 * (0 < j)]);
				}
				else
				{
					//Use Kn
					a_tmp_ij = (F_VelocitySlip_tmp * Kn) / (0.5 * hy[j - 1 * (0 < j)]);
					a_tmp_i1j = (F_VelocitySlip_tmp * Kn) / (0.5 * hy[j - 1 * (0 < j)]);
				}

				//Calculate velocities:
				//Here are made two calculations, but this is simplest case.
				//If the check are made to make one calculation, the checks are much more.
				if(vol_inf_fb_bool[i_1j_1])
					u[ij] = (u[ij_1] * a_tmp_ij + u_body_tmp) / (1.0 + a_tmp_ij);

				if(vol_inf_fb_bool[i1j_1])
					u[i1j] = (u[i1j_1] * a_tmp_i1j + u_body_tmp) / (1.0 + a_tmp_i1j);
			}
			else
			{
				//In this case fluid is over the wall
				if(To_Use_Kn_local_in_wall_BC)
				{
					//Use Kn local
					a_tmp_ij = (F_VelocitySlip_tmp * Kn) / (rho_u[ij1] * 0.5 * hy[j + 1 * (j < (Ny - 1))]);
					a_tmp_i1j = (F_VelocitySlip_tmp * Kn) / (rho_u[i1j1] * 0.5 * hy[j + 1 * (j < (Ny - 1))]);
				}
				else
				{
					//Use Kn
					a_tmp_ij = (F_VelocitySlip_tmp * Kn) / (0.5 * hy[j + 1 * (j < (Ny - 1))]);
					a_tmp_i1j = (F_VelocitySlip_tmp * Kn) / (0.5 * hy[j + 1 * (j < (Ny - 1))]);
				}

				//Calculate velocities:
				//Here are made two calculations, but this is simplest case.
				//If the check are made to make one calculation, the checks are much more.
				if(vol_inf_fb_bool[i_1j1])
					u[ij] = (u[ij1] * a_tmp_ij + u_body_tmp) / (1.0 + a_tmp_ij);

				if(vol_inf_fb_bool[i1j1])
					u[i1j] = (u[i1j1] * a_tmp_i1j + u_body_tmp) / (1.0 + a_tmp_i1j);
			}

		}


	}

}


inline void BoundaryConditions_Velocities_SlipBC_v(void)
{
	//Velocity Slip Boundary Conditions for v
	//Here have to be checked all calculated nodes in this SubDomain.
	//Because of the case, when we have channel wall only in i == 0,
	//we have to check this node. For the case, when we have SubDomains on
	//surface boundary this is also good, because the internal points of SubDomain
	//are calculated using point from other SubDomain.
	for(counter = 0; counter < Na_BoundaryConditions_Velocities_SlipBC_v; counter++)
	{
		node = rule_vector_ij_BoundaryConditions_Velocities_SlipBC_v[counter];
		j = rule_vector_j_BoundaryConditions_Velocities_SlipBC_v[counter];
		i = node - j * Nx;

		if(Periodic_boundary_conditions_about_OX)
		{
			ij_function_PeriodicBC_i_and_j_to_the_boundaries(i, j);
		}
		else
		{
			ij_function_i_and_j_to_the_boundaries(node);
		}

		is_VelocitySlipBC_body_tmp = false;

		//Now check the kind of body: circle or polyhedron
		if((vol_inf_fb[ij] >= Ncp_beg_b)
			&& (vol_inf_fb[ij] < (Ncp_beg_b + OforS[gb].N_circles)))
		{
			//The body is circle
			//Velocity Slip Boundary Conditions for circles.
			//The number of circle in this node is:
			N_circle = vol_inf_fb[ij] - Ncp_beg_b;

			v_body_tmp = OforS[gb].M_circles[N_circle].v_body;

			is_VelocitySlipBC_body_tmp = OforS[gb].M_circles[N_circle].is_VelocitySlipBC_body;
			w_VelocitySlipBC_tmp = OforS[gb].M_circles[N_circle].w_VelocitySlipBC;
			F_VelocitySlip_tmp = OforS[gb].M_circles[N_circle].F_VelocitySlip;
		}
		else if(vol_inf_fb[ij] >= (Ncp_beg_b + OforS[gb].N_circles))
		{
			//Velocity Slip Boundary Conditions for polyhedrons.
			//The number of polihedron in this node is:
			N_polyhedron = vol_inf_fb[ij] - (Ncp_beg_b + OforS[gb].N_circles);

			v_body_tmp = OforS[gb].M_polyhedrons[N_polyhedron].v_body;

			is_VelocitySlipBC_body_tmp = OforS[gb].M_polyhedrons[N_polyhedron].is_VelocitySlipBC_body;
			w_VelocitySlipBC_tmp = OforS[gb].M_polyhedrons[N_polyhedron].w_VelocitySlipBC;
			F_VelocitySlip_tmp = OforS[gb].M_polyhedrons[N_polyhedron].F_VelocitySlip;
		}


		//Apply velocity slip BC
		//This condition here is enough, because if there is no body to this volume
		//is_VelocitySlipBC_body_tmp, will left false.
		if(is_VelocitySlipBC_body_tmp)
		{
			if(fluid == vol_inf_fb[i_1j])
			{
				//In this case fluid is back the wall
				if(To_Use_Kn_local_in_wall_BC)
				{
					//Use Kn local
					a_tmp_ij = (F_VelocitySlip_tmp * Kn) / (rho_v[i_1j] * 0.5 * hx[i_1]);
					a_tmp_ij1 = (F_VelocitySlip_tmp * Kn) / (rho_v[i_1j1] * 0.5 * hx[i_1]);
				}
				else
				{
					//Use Kn
					a_tmp_ij = (F_VelocitySlip_tmp * Kn) / (0.5 * hx[i_1]);
					a_tmp_ij1 = (F_VelocitySlip_tmp * Kn) / (0.5 * hx[i_1]);
				}

				//Calculate velocities:
				//Here are made two calculations, but this is simplest case.
				//If the check are made to make one calculation, the checks are much more.
				if(vol_inf_fb_bool[i_1j_1])
					v[ij] = (v[i_1j] * a_tmp_ij + v_body_tmp) / (1.0 + a_tmp_ij);

				if(vol_inf_fb_bool[i_1j1])
					v[ij1] = (v[i_1j1] * a_tmp_ij1 + v_body_tmp) / (1.0 + a_tmp_ij1);
			}
			else
			{
				//In this case fluid is in front of the wall
				if(To_Use_Kn_local_in_wall_BC)
				{
					//Use Kn local
					a_tmp_ij = (F_VelocitySlip_tmp * Kn) / (rho_v[i1j] * 0.5 * hx[i1]);
					a_tmp_ij1 = (F_VelocitySlip_tmp * Kn) / (rho_v[i1j1] * 0.5 * hx[i1]);
				}
				else
				{
					//Use Kn
					a_tmp_ij = (F_VelocitySlip_tmp * Kn) / (0.5 * hx[i1]);
					a_tmp_ij1 = (F_VelocitySlip_tmp * Kn) / (0.5 * hx[i1]);
				}

				//Calculate velocities:
				//Here are made two calculations, but this is simplest case.
				//If the check are made to make one calculation, the checks are much more.
				if(vol_inf_fb_bool[i1j_1])
					v[ij] = (v[i1j] * a_tmp_ij + v_body_tmp) / (1.0 + a_tmp_ij);

				if(vol_inf_fb_bool[i1j1])
					v[ij1] = (v[i1j1] * a_tmp_ij1 + v_body_tmp) / (1.0 + a_tmp_ij1);
			}
		}


	}


}



inline void BoundaryConditions_Pressure(void)
/*
Target:	For this element are applyed Bounary Conditions For Pressure In Solved Area.
*/
{
	if(!Periodic_boundary_conditions_about_OX)
	{
		//dpdx_0_BC_xb give dp / dx = 0, on x = x_b
		if(dpdx_0_BC_xb
#ifdef COMPILE_MPI
				&& (!SubDomain[IJ].is_SubDomain_in_direction_I_1)
#endif //End of COMPILE_MPI
			)
		{
			for (j = 0; j < Ny; j++)
			{
				//boundary conditions for others in x_b
				node = Nx_others[0] + j * Nx;
				for (i = 0; i < Nx_others[0]; i++)
				{
					if(vol_inf_fb_bool[node])
						p[i + j * Nx] = p[node];
				}
			}
		}

		//dpdx_0_BC_xe give dp / dx = 0, on x = x_e
		if(dpdx_0_BC_xe
#ifdef COMPILE_MPI
				&& (!SubDomain[IJ].is_SubDomain_in_direction_I1)
#endif //End of COMPILE_MPI
			)
		{
			for (j = 0; j < Ny; j++)
			{
				//boundary conditions for others in x_e
				node = Nx_others[1] - 1 + j * Nx;
				for (i = Nx_others[1]; i < Nx; i++)
				{
					if(vol_inf_fb_bool[node])
					{
						//dpdx=0
						p[i + j * Nx] = p[node];
					}

				}

			}
		}

		//dpdy_0_BC_yb give dp / dy = 0, on y = y_b
		if(dpdy_0_BC_yb
#ifdef COMPILE_MPI
				&& (!SubDomain[IJ].is_SubDomain_in_direction_J_1)
#endif //End of COMPILE_MPI
			)
		{
			for (i = 0; i < Nx; i++)
			{
				//boundary conditions for others in y_b
				node = i + Ny_others[0] * Nx;
				for (j = 0; j < Ny_others[0]; j++)
				{
					if(vol_inf_fb_bool[node])
						p[i + j * Nx] = p[node];
				}
			}
		}

		//dpdy_0_BC_ye give dp / dy = 0, on y = y_e
		if(dpdy_0_BC_ye
#ifdef COMPILE_MPI
				&& (!SubDomain[IJ].is_SubDomain_in_direction_J1)
#endif //End of COMPILE_MPI
			)
		{
			for (i = 0; i < Nx; i++)
			{
				//boundary conditions for others in y_b
				node = i + (Ny_others[1] - 1) * Nx;
				for (j = Ny_others[1]; j < Ny; j++)
				{
					if(vol_inf_fb_bool[node])
						p[i + j * Nx] = p[node];
				}

			}
		}

	////in BC with given velocity must have one node given pressure inside area
	//p[Nx_others[0] + 5 * Nx] = p_gas_nd;
	}

}

inline void BoundaryConditions_Pressure_begin_iteration_for_it(void)
/*
This boundary condition is applyed when we start iteration for it.
In begin of every iteration
*/
{
	if(!Periodic_boundary_conditions_about_OX)
	{
#ifdef COMPILE_MPI
		//This is prevent, because u_given_xb_in_it is needed
		//to calculate cT and Drag Coefficient (CD).
		u_given_xb_in_it = u_given_xb;
#endif //End of COMPILE_MPI



#ifdef COMPILE_MPI
			//This BC is used only for the inflow of the channel.
		if(!SubDomain[IJ].is_SubDomain_in_direction_I_1)
	#endif //End of COMPILE_MPI
		{
			//layer on OX where we solve mass inflow
			i = Nx_u[0];

			//find u_max
			u_given_xb_in_it_max = u[0 + 0 * Nx];
			for(j = 0; j < Ny; j++)
			{
				node = i + j * Nx;

				if(vol_inf_fb_bool[node])
				{
					u_given_xb_in_it_max = maximum(u_given_xb_in_it_max, u[node]);
				}
			}

			double h_channel_inflow;

			//layer on OX where we solve mass inflow
			i = Nx_u[0];

			u_given_xb_in_it_mean = 0;
			h_channel_inflow = 0;
			counter = 0;
			for(j = 0; j < Ny; j++)
			{
				node = i + j * Nx;

				if(vol_inf_fb_bool[node])
				{
					u_given_xb_in_it_mean += u[node] * hy[j];
					h_channel_inflow += hy[j];
					counter++;
				}
			}

			u_given_xb_in_it_mean = u_given_xb_in_it_mean / h_channel_inflow;

			//Solve given velocity if we give Re and Kn numbers
			if(GivenReKn)
			{
				u_given_xb = Re * Kn * sqrt(15.0 * M_PI / 128.0);
			}


			//choose fixed velocity
			if(p_BC_xb_correction_method == 0)
			{
				u_given_xb_in_it = u_given_xb;

				//We need pressure in the begining of solved area to define density in infinity.
				//We need that for compressible fluid only!
				//We need that only if we have Given_velocity_on_xb.
				if(Given_velocity_on_xb)
				{
					//The pressure which is solved.
					i = Nx_others[0];

					//The pressure near to the center of the solve area on OY axis.
					j = int((Ny_others[1] - 1.0) - Ny_others[0]);

					p_BC_xb = p[i + j * Nx];
				}
			}
			else if(p_BC_xb_correction_method == p_BC_xb_correction_method_Vmax_xb)
			{
				u_given_xb_in_it = u_given_xb_in_it_max;
			}
			else if(p_BC_xb_correction_method == p_BC_xb_correction_method_Vmean_xb)
			{
				u_given_xb_in_it = u_given_xb_in_it_mean;
			}


			if(Pressure_BC
				//&& (p_BC_xb_correction_method > 0)
				//&& (Iter_p_c == 0)
				)
			{

				if(Pressure_ratio_correct)
				{
					//correct p_correction
					if(p_correction_auto)
					{
						double p_correction_tmp;

						if(Iter_pr1 <= Iter_pr2 && Iter_pr1 < (N_I - 1))
							p_correction_tmp = p_correction * (1.0 + sqrt(max_residual_in_u * max_residual_in_u + max_residual_in_v * max_residual_in_v));
						else
							p_correction_tmp = p_correction * (1.0 - sqrt(max_residual_in_u * max_residual_in_u + max_residual_in_v * max_residual_in_v));

						if(p_correction_tmp < p_correction_min)
							p_correction_tmp = p_correction_min;
						else if(p_correction_max < p_correction_tmp)
							p_correction_tmp = p_correction_max;

						p_correction = p_correction_tmp;

					}


					//maximum variation of p_BC_xb per iteration is +- p_x_b_correction
					if(u_given_xb_in_it < (u_given_xb * (1.0 - u_given_xb_error)))
					{
						if(correct_p_BC_xb)
						{
							p_BC_xb = p_BC_xb * (1.0 + p_correction * (u_given_xb - u_given_xb_in_it)* (u_given_xb - u_given_xb_in_it));
						}
						else
						{
							double p_BC_tmp;
							p_BC_tmp = p_BC_xe * (1.0 - p_correction * (u_given_xb - u_given_xb_in_it)* (u_given_xb - u_given_xb_in_it));

							if(p_BC_tmp < 0.0)
								p_BC_xe = 0.0;
							else
								p_BC_xe = p_BC_tmp;
						}
					}
					else if((u_given_xb * (1.0 + u_given_xb_error)) < u_given_xb_in_it)
					{
						if(correct_p_BC_xb)
						{
							double p_BC_tmp;
							p_BC_tmp = p_BC_xb * (1.0 - p_correction * (u_given_xb - u_given_xb_in_it)* (u_given_xb - u_given_xb_in_it));

							//p_XB_xb > p_BC_xe - if we want to have pressure driven flow
							if(p_BC_tmp > p_BC_xe)
								p_BC_xb = p_BC_tmp;
							else
								p_BC_xb = p_BC_xe;
						}
						else
						{
							double p_BC_tmp;
							p_BC_tmp = p_BC_xb * (1.0 + p_correction * (u_given_xb - u_given_xb_in_it)* (u_given_xb - u_given_xb_in_it));

							if(p_BC_tmp > p_BC_xb)
								p_BC_xe = p_BC_xb;
							else
								p_BC_xe = p_BC_tmp;

						}
					}


					//Apply pressure boundary conditions
					if(correct_p_BC_xb)
					{
						//Inflow BC -> x = x_b
						for(j = 0; j < Ny; j++)
						{
							for(i = 0; i < Nx_others[0]; i++)
							{
								node = i + j * Nx;
								if(vol_inf_fb_bool[node])
								{
									p[node] = p_BC_xb;
									rho[node] = p_BC_xb / T_BC_xb;
								}

							}
						}

					}
					else
					{
						//Outflow BC -> x = x_e
						for(j = 0; j < Ny; j++)
						{
							for(i = Nx_others[1]; i < Nx; i++)
							{
								node = i + j * Nx;
								if(vol_inf_fb_bool[node])
								{
									p[node] = p_BC_xe;
									rho[node] = p_BC_xe / T_BC_xe;
								}

							}
						}

					}
					////WriteMassiveToBinaryFile("p.bin", p, Na);

				}


			}


		}

	}


}


inline void BoundaryConditions_Temperature(void)
/*
Target:	For this element are applyed Bounary Conditions For Temperature In Solved Area.
*/
{
	if(!Periodic_boundary_conditions_about_OX)
	{
		//dTdx_0_BC_xb give dTemper / dx = 0, on x = x_b
		if(dTdx_0_BC_xb
#ifdef COMPILE_MPI
				&& (!SubDomain[IJ].is_SubDomain_in_direction_I_1)
#endif //End of COMPILE_MPI
			)
		{
			for (j = 0; j < Ny; j++)
			{
				//boundary conditions for others in x_b
				node = Nx_others[0] + j * Nx;
				for (i = 0; i < Nx_others[0]; i++)
				{
					if(vol_inf_fb_bool[node])
						Temper[i + j * Nx] = Temper[node];
				}
			}
		}

		//dTdx_0_BC_xe give dTemper / dx = 0, on x = x_e
		if(dTdx_0_BC_xe
#ifdef COMPILE_MPI
				&& (!SubDomain[IJ].is_SubDomain_in_direction_I1)
#endif //End of COMPILE_MPI
			)
		{
			for (j = 0; j < Ny; j++)
			{
				//boundary conditions for others in x_e
				node = Nx_others[1] - 1 + j * Nx;
				for (i = Nx_others[1]; i < Nx; i++)
				{
					if(vol_inf_fb_bool[node])
					{
						//dTdx=0;
						Temper[i + j * Nx] = Temper[node];
					}

				}

			}
		}

		//dTdy_0_BC_yb give dTemper / dy = 0, on y = y_b
		if(dTdy_0_BC_yb
#ifdef COMPILE_MPI
				&& (!SubDomain[IJ].is_SubDomain_in_direction_J_1)
#endif //End of COMPILE_MPI
			)
		{
			for (i = 0; i < Nx; i++)
			{
				//boundary conditions for others in y_b
				node = i + Ny_others[0] * Nx;
				for (j = 0; j < Ny_others[0]; j++)
				{
					if(vol_inf_fb_bool[node])
						Temper[i + j * Nx] = Temper[node];
				}
			}
		}

		//dTdy_0_BC_ye give dTemper / dy = 0, on y = y_e
		if(dTdy_0_BC_ye
#ifdef COMPILE_MPI
				&& (!SubDomain[IJ].is_SubDomain_in_direction_J1)
#endif //End of COMPILE_MPI
			)
		{
			for (i = 0; i < Nx; i++)
			{
				//boundary conditions for others in y_b
				node = i + (Ny_others[1] - 1) * Nx;
				for (j = Ny_others[1]; j < Ny; j++)
				{
					if(vol_inf_fb_bool[node])
						Temper[i + j * Nx] = Temper[node];
				}

			}
		}

	}

}


inline void BoundaryConditions_Temperature_TemperatureJumpBC(void)
{
	for(counter = 0; counter < Na_BoundaryConditions_Temperature_TemperatureJumpBC; counter++)
	{
		node = rule_vector_ij_BoundaryConditions_Temperature_TemperatureJumpBC[counter];
		j = rule_vector_j_BoundaryConditions_Temperature_TemperatureJumpBC[counter];
		i = node - j * Nx;

		if(Periodic_boundary_conditions_about_OX)
		{
			ij_function_PeriodicBC_i_and_j_to_the_boundaries(i, j);
		}
		else
		{
			ij_function_i_and_j_to_the_boundaries(node);
		}


		//Check for body in node vol_inf_fb_bool[i_1j]
		if(!vol_inf_fb_bool[i_1j])
		{
			//There is body in volume i_1j
			//Temperature jump BC for external vertex:
			is_TemperatureJumpBC_body_tmp = false;

			//Now check the kind of body: circle or polyhedron
			if((vol_inf_fb[i_1j] >= Ncp_beg_b)
				&& (vol_inf_fb[i_1j] < (Ncp_beg_b + OforS[gb].N_circles)))
			{
				//The body is circle
				//The number of circle in this node is:
				N_circle = vol_inf_fb[i_1j] - Ncp_beg_b;

				T_body_tmp = OforS[gb].M_circles[N_circle].T_body;

				is_TemperatureJumpBC_body_tmp = OforS[gb].M_circles[N_circle].is_TemperatureJumpBC_body;
				dTdn_on_wall_0_tmp = OforS[gb].M_circles[N_circle].dTdn_on_wall_0;
				F_TemperatureJump_tmp = OforS[gb].M_circles[N_circle].F_TemperatureJump;
			}
			else if(vol_inf_fb[i_1j] >= (Ncp_beg_b + OforS[gb].N_circles))
			{
				//Temperature Jump Boundary Conditions for polyhedrons.
				//The number of polihedron in this node is:
				N_polyhedron = vol_inf_fb[i_1j] - (Ncp_beg_b + OforS[gb].N_circles);

				T_body_tmp = OforS[gb].M_polyhedrons[N_polyhedron].T_body;

				is_TemperatureJumpBC_body_tmp = OforS[gb].M_polyhedrons[N_polyhedron].is_TemperatureJumpBC_body;
				dTdn_on_wall_0_tmp = OforS[gb].M_polyhedrons[N_polyhedron].dTdn_on_wall_0;
				F_TemperatureJump_tmp = OforS[gb].M_polyhedrons[N_polyhedron].F_TemperatureJump;
			}


			if(is_TemperatureJumpBC_body_tmp)
			{
				//If the Kn is not local this is true for all variants,
				//but F_TemperatureJump_tmp is different for all variants
				//and have to calculate this on every variant.
				if(To_Use_Kn_local_in_wall_BC)
				{
					F_TemperatureJump_Kn_for_this_BC_tmp = (F_TemperatureJump_tmp * Kn) / rho[ij];
				}
				else
				{
					F_TemperatureJump_Kn_for_this_BC_tmp = F_TemperatureJump_tmp * Kn;
				}

				//Apply BC for temperature jump according:
				a_tmp = F_TemperatureJump_Kn_for_this_BC_tmp / (0.5 * hx[i]);

				//From dimentional equation for Slip BC case NonDimentionalGiven_Kn_pch_rho_R_T
				Temper_of_fluid_on_the_wall_u[ij] = (a_tmp * Temper[ij] + T_body_tmp) / (1.0 + a_tmp);
			}
			else if(dTdn_on_wall_0_tmp)
			{
				//If normal derivative of the temperature on the wall is zero
				Temper_of_fluid_on_the_wall_u[ij] = Temper[ij];
			}
			else
			{
				//There is no temperature jump on this surface and tempetarure of the fluid on the
				//surface is equal to the temperature of the body.
				Temper_of_fluid_on_the_wall_u[ij] = T_body_tmp;
			}

		}


		//Check for body in node vol_inf_fb_bool[i1j]
		if(!vol_inf_fb_bool[i1j])
		{
			//There is body in volume i1j
			//Temperature jump BC for external vertex:
			is_TemperatureJumpBC_body_tmp = false;

			//Now check the kind of body: circle or polyhedron
			if((vol_inf_fb[i1j] >= Ncp_beg_b)
				&& (vol_inf_fb[i1j] < (Ncp_beg_b + OforS[gb].N_circles)))
			{
				//The body is circle
				//The number of circle in this node is:
				N_circle = vol_inf_fb[i1j] - Ncp_beg_b;

				T_body_tmp = OforS[gb].M_circles[N_circle].T_body;

				is_TemperatureJumpBC_body_tmp = OforS[gb].M_circles[N_circle].is_TemperatureJumpBC_body;
				dTdn_on_wall_0_tmp = OforS[gb].M_circles[N_circle].dTdn_on_wall_0;
				F_TemperatureJump_tmp = OforS[gb].M_circles[N_circle].F_TemperatureJump;
			}
			else if(vol_inf_fb[i1j] >= (Ncp_beg_b + OforS[gb].N_circles))
			{
				//Temperature Jump Boundary Conditions for polyhedrons.
				//The number of polihedron in this node is:
				N_polyhedron = vol_inf_fb[i1j] - (Ncp_beg_b + OforS[gb].N_circles);

				T_body_tmp = OforS[gb].M_polyhedrons[N_polyhedron].T_body;

				is_TemperatureJumpBC_body_tmp = OforS[gb].M_polyhedrons[N_polyhedron].is_TemperatureJumpBC_body;
				dTdn_on_wall_0_tmp = OforS[gb].M_polyhedrons[N_polyhedron].dTdn_on_wall_0;
				F_TemperatureJump_tmp = OforS[gb].M_polyhedrons[N_polyhedron].F_TemperatureJump;
			}


			if(is_TemperatureJumpBC_body_tmp)
			{
				//If the Kn is not local this is true for all variants,
				//but F_TemperatureJump_tmp is different for all variants
				//and have to calculate this on every variant.
				if(To_Use_Kn_local_in_wall_BC)
				{
					F_TemperatureJump_Kn_for_this_BC_tmp = (F_TemperatureJump_tmp * Kn) / rho[ij];
				}
				else
				{
					F_TemperatureJump_Kn_for_this_BC_tmp = F_TemperatureJump_tmp * Kn;
				}

				//Apply BC for temperature jump according:
				a_tmp = F_TemperatureJump_Kn_for_this_BC_tmp / (0.5 * hx[i]);

				//From dimentional equation for Slip BC case NonDimentionalGiven_Kn_pch_rho_R_T
				Temper_of_fluid_on_the_wall_u[i1j] = (a_tmp * Temper[ij] + T_body_tmp) / (1.0 + a_tmp);
			}
			else if(dTdn_on_wall_0_tmp)
			{
				//If normal derivative of the temperature on the wall is zero
				Temper_of_fluid_on_the_wall_u[i1j] = Temper[ij];
			}
			else
			{
				//There is no temperature jump on this surface and tempetarure of the fluid on the
				//surface is equal to the temperature of the body.
				Temper_of_fluid_on_the_wall_u[i1j] = T_body_tmp;
			}

		}


		//Check for body in node vol_inf_fb_bool[ij_1]
		if(!vol_inf_fb_bool[ij_1])
		{
			//There is body in volume ij_1
			//Temperature jump BC for external vertex:
			is_TemperatureJumpBC_body_tmp = false;

			//Now check the kind of body: circle or polyhedron
			if((vol_inf_fb[ij_1] >= Ncp_beg_b)
				&& (vol_inf_fb[ij_1] < (Ncp_beg_b + OforS[gb].N_circles)))
			{
				//The body is circle
				//The number of circle in this node is:
				N_circle = vol_inf_fb[ij_1] - Ncp_beg_b;

				T_body_tmp = OforS[gb].M_circles[N_circle].T_body;

				is_TemperatureJumpBC_body_tmp = OforS[gb].M_circles[N_circle].is_TemperatureJumpBC_body;
				dTdn_on_wall_0_tmp = OforS[gb].M_circles[N_circle].dTdn_on_wall_0;
				F_TemperatureJump_tmp = OforS[gb].M_circles[N_circle].F_TemperatureJump;
			}
			else if(vol_inf_fb[ij_1] >= (Ncp_beg_b + OforS[gb].N_circles))
			{
				//Temperature Jump Boundary Conditions for polyhedrons.
				//The number of polihedron in this node is:
				N_polyhedron = vol_inf_fb[ij_1] - (Ncp_beg_b + OforS[gb].N_circles);

				T_body_tmp = OforS[gb].M_polyhedrons[N_polyhedron].T_body;

				is_TemperatureJumpBC_body_tmp = OforS[gb].M_polyhedrons[N_polyhedron].is_TemperatureJumpBC_body;
				dTdn_on_wall_0_tmp = OforS[gb].M_polyhedrons[N_polyhedron].dTdn_on_wall_0;
				F_TemperatureJump_tmp = OforS[gb].M_polyhedrons[N_polyhedron].F_TemperatureJump;
			}


			if(is_TemperatureJumpBC_body_tmp)
			{
				//If the Kn is not local this is true for all variants,
				//but F_TemperatureJump_tmp is different for all variants
				//and have to calculate this on every variant.
				if(To_Use_Kn_local_in_wall_BC)
				{
					F_TemperatureJump_Kn_for_this_BC_tmp = (F_TemperatureJump_tmp * Kn) / rho[ij];
				}
				else
				{
					F_TemperatureJump_Kn_for_this_BC_tmp = F_TemperatureJump_tmp * Kn;
				}

				//Apply BC for temperature jump according:
				a_tmp = F_TemperatureJump_Kn_for_this_BC_tmp / (0.5 * hy[j]);

				//From dimentional equation for Slip BC case NonDimentionalGiven_Kn_pch_rho_R_T
				Temper_of_fluid_on_the_wall_v[ij] = (a_tmp * Temper[ij] + T_body_tmp) / (1.0 + a_tmp);
			}
			else if(dTdn_on_wall_0_tmp)
			{
				//If normal derivative of the temperature on the wall is zero
				Temper_of_fluid_on_the_wall_v[ij] = Temper[ij];
			}
			else
			{
				//There is no temperature jump on this surface and tempetarure of the fluid on the
				//surface is equal to the temperature of the body.
				Temper_of_fluid_on_the_wall_v[ij] = T_body_tmp;
			}

		}


		//Check for body in node vol_inf_fb_bool[ij1]
		if(!vol_inf_fb_bool[ij1])
		{
			//There is body in volume ij1
			//Temperature jump BC for external vertex:
			is_TemperatureJumpBC_body_tmp = false;

			//Now check the kind of body: circle or polyhedron
			if((vol_inf_fb[ij1] >= Ncp_beg_b)
				&& (vol_inf_fb[ij1] < (Ncp_beg_b + OforS[gb].N_circles)))
			{
				//The body is circle
				//The number of circle in this node is:
				N_circle = vol_inf_fb[ij1] - Ncp_beg_b;

				T_body_tmp = OforS[gb].M_circles[N_circle].T_body;

				is_TemperatureJumpBC_body_tmp = OforS[gb].M_circles[N_circle].is_TemperatureJumpBC_body;
				dTdn_on_wall_0_tmp = OforS[gb].M_circles[N_circle].dTdn_on_wall_0;
				F_TemperatureJump_tmp = OforS[gb].M_circles[N_circle].F_TemperatureJump;
			}
			else if(vol_inf_fb[ij1] >= (Ncp_beg_b + OforS[gb].N_circles))
			{
				//Temperature Jump Boundary Conditions for polyhedrons.
				//The number of polihedron in this node is:
				N_polyhedron = vol_inf_fb[ij1] - (Ncp_beg_b + OforS[gb].N_circles);

				T_body_tmp = OforS[gb].M_polyhedrons[N_polyhedron].T_body;

				is_TemperatureJumpBC_body_tmp = OforS[gb].M_polyhedrons[N_polyhedron].is_TemperatureJumpBC_body;
				dTdn_on_wall_0_tmp = OforS[gb].M_polyhedrons[N_polyhedron].dTdn_on_wall_0;
				F_TemperatureJump_tmp = OforS[gb].M_polyhedrons[N_polyhedron].F_TemperatureJump;
			}

			if(is_TemperatureJumpBC_body_tmp)
			{
				//If the Kn is not local this is true for all variants,
				//but F_TemperatureJump_tmp is different for all variants
				//and have to calculate this on every variant.
				if(To_Use_Kn_local_in_wall_BC)
				{
					F_TemperatureJump_Kn_for_this_BC_tmp = (F_TemperatureJump_tmp * Kn) / rho[ij];
				}
				else
				{
					F_TemperatureJump_Kn_for_this_BC_tmp = F_TemperatureJump_tmp * Kn;
				}

				//Apply BC for temperature jump according:
				a_tmp = F_TemperatureJump_Kn_for_this_BC_tmp / (0.5 * hy[j]);

				//From dimentional equation for Slip BC case NonDimentionalGiven_Kn_pch_rho_R_T
				Temper_of_fluid_on_the_wall_v[ij1] = (a_tmp * Temper[ij] + T_body_tmp) / (1.0 + a_tmp);
			}
			else if(dTdn_on_wall_0_tmp)
			{
				//If normal derivative of the temperature on the wall is zero
				Temper_of_fluid_on_the_wall_v[ij1] = Temper[ij];
			}
			else
			{
				//There is no temperature jump on this surface and tempetarure of the fluid on the
				//surface is equal to the temperature of the body.
				Temper_of_fluid_on_the_wall_v[ij1] = T_body_tmp;
			}

		}



	}


}


inline void BoundaryConditions_Temperature_TemperatureJumpBC_WriteToMassiveToWriteDataToFile(void)
{
	BoundaryConditions_Temperature_TemperatureJumpBC();

	for(counter = 0; counter < Na_BoundaryConditions_Temperature_TemperatureJumpBC; counter++)
	{
		node = rule_vector_ij_BoundaryConditions_Temperature_TemperatureJumpBC[counter];
		j = rule_vector_j_BoundaryConditions_Temperature_TemperatureJumpBC[counter];
		i = node - j * Nx;

		if(Periodic_boundary_conditions_about_OX)
		{
			ij_function_PeriodicBC_i_and_j_to_the_boundaries(i, j);
		}
		else
		{
			ij_function_i_and_j_to_the_boundaries(node);
		}

		//The following assumptions are not correct for external corner, because
		//to the temperature on the control volume on the external vertex are copied
		//two temperatures on the two surfaces.
		//Here this approach is used as information only.
		if(!vol_inf_fb_bool[i_1j])
		{
			//The temperature of the surface is copied to the temperature in the body
			Temper[i_1j] = Temper_of_fluid_on_the_wall_u[ij];
		}

		if(!vol_inf_fb_bool[i1j])
		{
			//The temperature of the surface is copied to the temperature in the body
			Temper[i1j] = Temper_of_fluid_on_the_wall_u[i1j];
		}

		if(!vol_inf_fb_bool[ij_1])
		{
			//The temperature of the surface is copied to the temperature in the body
			Temper[ij_1] = Temper_of_fluid_on_the_wall_v[ij];
		}

		if(!vol_inf_fb_bool[ij1])
		{
			//The temperature of the surface is copied to the temperature in the body
			Temper[ij1] = Temper_of_fluid_on_the_wall_v[ij1];
		}


	}


}


inline void BoundaryConditions_function_of_time(void)
//This procedure change the inflow velocity of the channel.
{
	if(!Periodic_boundary_conditions_about_OX)
	{
		//Change BC for u
		if(use_interpolation_for_u_BC_by_time
#ifdef COMPILE_MPI
		//This BC is used only for the inflow of the channel.
		&& (!SubDomain[IJ].is_SubDomain_in_direction_I_1)
#endif //End of COMPILE_MPI
			)
		{
			if(it <= t_BC_b)
				u_gas_nd = u_BC_xb_t_BC_b;
			else if(t_BC_b < it && it <= t_BC_e)
			{
				//Linear interpolation of BC
				u_gas_nd = ((u_BC_xb_t_BC_e - u_BC_xb_t_BC_b) / (t_BC_e - t_BC_b)) * (it - t_BC_b) + u_BC_xb_t_BC_b;
			}
			else
				u_gas_nd = u_BC_xb_t_BC_e;


			u_given_xb = u_gas_nd;


			//Apply BC u = u_gas_nd, for x = x_b
			for (j = 0; j < Ny; j++)
			{
				for (i = 0; i < Nx_u[0]; i++)
				{
					//boundary conditions for u in xb
					if(vol_inf_fb_bool[i + j * Nx])
						u[i + j * Nx] = u_gas_nd;
				}
			}


		}
	}

}


void WriteDragCoefficient(void)
{
	CD_fr_x = 0;
	CD_p_x = 0;

	//I will count number of bodies
	unsigned int NumberOfBodies;

	NumberOfBodies = OforS[gb].N_circles + OforS[gb].N_polyhedrons;


	bool ToSolveDragCoefficientForThisBody;
	//we count all count all bodies and will check which coefficient will be solved
	for(counter = 0; counter < NumberOfBodies; counter++)
	{
		ToSolveDragCoefficientForThisBody = false;

		//check which body drag coefficent will solve
		//check for circles
		if(counter < OforS[gb].N_circles)
		{
			ToSolveDragCoefficientForThisBody =
				OforS[gb].M_circles[counter].ToSolveDragCoefficientForThisBody;
		}
		else if(OforS[gb].is_enter_M_polyhedrons)
		{
			ToSolveDragCoefficientForThisBody =
				OforS[gb].M_polyhedrons[counter - OforS[gb].N_circles].ToSolveDragCoefficientForThisBody;
		}


		//solve drag coefficent
		if(ToSolveDragCoefficientForThisBody)
		{
			using namespace std;

			//Solve Drag Coefficient CD
			CD_fr_bottom = 0;
			CD_p_bottom = 0;
			S_for_CD_bottom = 0;
			counter_bottom = 0;

			CD_fr_front = 0;
			CD_p_front = 0;
			S_for_CD_front = 0;
			counter_front = 0;

			CD_fr_top = 0;
			CD_p_top = 0;
			S_for_CD_top = 0;
			counter_top = 0;

			CD_fr_behind = 0;
			CD_p_behind = 0;
			S_for_CD_behind = 0;
			counter_behind = 0;


			for(i = Nx_others[0]; i < Nx_others[1]; i++)
			{
				for(j = Ny_others[0]; j < Ny_others[1]; j++)
				{
					node = i + j * Nx;

					//Solve CD ONLY for body with number counter.
					if((vol_inf_fb[node] - 1) == counter)
					{
						ij_function_i_and_j_to_the_boundaries(node);

						if(vol_inf_fb_bool[ij_1])
						{
							CD_fr_bottom +=
								hx[i] * sqrt(Temper_of_fluid_on_the_wall_v[ij])
										//sqrt_T[ij]
								* ((u[ij_1] - u[ij]) / hy[j - 1 * (j < (Ny - 1))]
									+ (u[i1j_1] - u[i1j]) / hy[j - 1 * (j < (Ny - 1))]);

							CD_p_bottom += p[ij_1] * hx[i];

							//Solve bottom surface of rectangular
							S_for_CD_bottom += hx[i];

							counter_bottom++;
						}

						if(vol_inf_fb_bool[i_1j])
						{
							CD_fr_front +=
								hy[j] * sqrt(Temper_of_fluid_on_the_wall_u[ij])
										//sqrt_T[ij]
								* ((v[i_1j] - v[ij]) / hx[i_1]
									+ (v[i_1j1] - v[ij1]) / hx[i_1]);

							CD_p_front += p[i_1j] * hy[j];

							//Solve front surface of rectabgular
							S_for_CD_front += hy[j];

							counter_front++;
						}

						if(vol_inf_fb_bool[ij1])
						{
							CD_fr_top +=
								hx[i] * sqrt(Temper_of_fluid_on_the_wall_v[ij1])
										//sqrt_T[ij]
								* ((u[ij1] - u[ij]) / hy[j + 1 * (0 < j)]
									+ (u[i1j1] - u[i1j]) / hy[j + 1 * (0 < j)]);

							CD_p_top += p[ij1] * hx[i];

							//Solve top surface of rectabgular
							S_for_CD_top += hx[i];

							counter_top++;
						}

						if(vol_inf_fb_bool[i1j])
						{
							CD_fr_behind +=
								hy[j] * sqrt(Temper_of_fluid_on_the_wall_u[i1j])
										//sqrt_T[ij]
								* ((v[i1j] - v[ij]) / hx[i1]
									+ (v[i1j1] - v[ij1]) / hx[i1]);

							CD_p_behind += p[i1j] * hy[j];

							//Solve behind surface of rectabgular
							S_for_CD_behind += hy[j];


							counter_behind++;
						}


					}


				}
			}

			//Multiplication coefficient for pressure term.
			double c_CD_p;

			//Multiplication coefficient for friction term.
			double c_CD_fr;

			//Use diferent type of non domentional volues.
			if(TypeOfNonDimensiolization == NonDimentionalGiven_Re_pch_rho_V2_2)
			{
				c_CD_p = 1.0 / (u_given_xb_in_it * u_given_xb_in_it * (p_BC_xb / T_BC_xb));

				c_CD_fr = 2.0 / (Re * u_given_xb_in_it);
			}
			else if(TypeOfNonDimensiolization == NonDimentionalGiven_Re_pch_rho_R_T)
			{
				c_CD_p = 2.0 / (gamma1 * Ma * Ma * (p_BC_xb / T_BC_xb));

				c_CD_fr = 2.0 / (Re * u_given_xb_in_it);
			}
			else if(TypeOfNonDimensiolization == NonDimentionalGiven_Kn_pch_rho_R_T)
			{
				c_CD_p = 1.0 / (u_given_xb_in_it * u_given_xb_in_it * (p_BC_xb / T_BC_xb));

				c_CD_fr = (5.0 / 8.0) * sqrt(M_PI) * Kn
					/ (u_given_xb_in_it * u_given_xb_in_it * (p_BC_xb / T_BC_xb));
			}


			//correct drag coefficent, becouse the area is NOT exact 1
			CD_fr_bottom = c_CD_fr * CD_fr_bottom;
			CD_p_bottom = c_CD_p * CD_p_bottom;

			CD_fr_front = c_CD_fr * CD_fr_front;
			CD_p_front = c_CD_p * CD_p_front;

			CD_fr_top = c_CD_fr * CD_fr_top;
			CD_p_top = c_CD_p * CD_p_top;

			CD_fr_behind = c_CD_fr * CD_fr_behind;
			CD_p_behind = c_CD_p * CD_p_behind;


#ifdef COMPILE_MPI
			//Exchange data for MPI program
			CollectDataForCDfromAllSubDomainsToProcessIJ0();

			if(IJ == 0) //Write data to file only from SubDomain IJ == 0
#endif //End of COMPILE_MPI
			{
				CD_fr_x = (CD_fr_top / S_for_CD_top) + (CD_fr_bottom / S_for_CD_bottom);
				CD_p_x = (CD_p_front / S_for_CD_front) - (CD_p_behind / S_for_CD_behind);

				string CD_file_name = "CD." + toString(counter) + ".txt";

				ofstream output_data;

				//if we continue the old solve the file will not be erased
				if(ToCheckForFileIsNow_CD)
				{
					if(ToReadSolvedDataFromFile)
						output_data.open(CD_file_name.c_str(), ios::out | ios::app);
					else
						output_data.open(CD_file_name.c_str());
				}
				else
					output_data.open(CD_file_name.c_str(), ios::out | ios::app);


				if(output_data.is_open())
				{
					int width_between_columns = 30;
					int Ndigits = 10;


					if(!ToReadSolvedDataFromFile && ToCheckForFileIsNow_CD)
					{
						output_data << std::setw(width_between_columns) << "//it (number of step by time) 00";
						output_data << std::setw(width_between_columns) << "ht (step by time) 01";
						output_data << std::setw(width_between_columns) << "CD_fr_x (CD from friction) 02";
						output_data << std::setw(width_between_columns) << "CD_p_x (CD from pressure) 03";

						//write exact area of surfaces
						output_data << std::setw(width_between_columns) << "S_for_CD_bottom 04";
						output_data << std::setw(width_between_columns) << "S_for_CD_front 05";
						output_data << std::setw(width_between_columns) << "S_for_CD_top 06";
						output_data << std::setw(width_between_columns) << "S_for_CD_behind 07";

						//write all coefficents
						output_data << std::setw(width_between_columns) << "CD_fr_bottom 08";
						output_data << std::setw(width_between_columns) << "CD_p_bottom 09";
						output_data << std::setw(width_between_columns) << "CD_fr_front 10";
						output_data << std::setw(width_between_columns) << "CD_p_front 11";
						output_data << std::setw(width_between_columns) << "CD_fr_top 12";
						output_data << std::setw(width_between_columns) << "CD_p_top 13";
						output_data << std::setw(width_between_columns) << "CD_fr_behind 14";
						output_data << std::setw(width_between_columns) << "CD_p_behind 15";

						output_data << endl;
					}

					output_data << std::setw(width_between_columns) << it;

					output_data
						<< std::setiosflags(std::ios::scientific)
						<< std::setprecision(Ndigits);

					output_data << std::setw(width_between_columns) << ht;
					output_data << std::setw(width_between_columns) << CD_fr_x;
					output_data << std::setw(width_between_columns) << CD_p_x;

					//write exact area of surfaces
					output_data << std::setw(width_between_columns) << S_for_CD_bottom;
					output_data << std::setw(width_between_columns) << S_for_CD_front;
					output_data << std::setw(width_between_columns) << S_for_CD_top;
					output_data << std::setw(width_between_columns) << S_for_CD_behind;

					//write all coefficents
					output_data << std::setw(width_between_columns) << CD_fr_bottom;
					output_data << std::setw(width_between_columns) << CD_p_bottom;
					output_data << std::setw(width_between_columns) << CD_fr_front;
					output_data << std::setw(width_between_columns) << CD_p_front;
					output_data << std::setw(width_between_columns) << CD_fr_top;
					output_data << std::setw(width_between_columns) << CD_p_top;
					output_data << std::setw(width_between_columns) << CD_fr_behind;
					output_data << std::setw(width_between_columns) << CD_p_behind;

					output_data << endl;


					output_data.close();


					ToCheckForFileIsNow_CD = false;

				}
				else
				{
					cout << "The program can not open file " << CD_file_name << " to write data." << endl;
				}

			}


		}
	}

#ifdef COMPILE_MPI
	if(IJ == 0) //Write data to file only from SubDomain IJ == 0
#endif //End of COMPILE_MPI
	{
		if(WriteDataInShortFormat_kind == 1)
		{
			//if WriteDataInShortFormat_kind == 1; ==> write: Re	Kn	u_gas_nd CD	CD_p

			//x_b (inflow) mass flow - q_xb
			//x_e (outflow) mass flow - q_xe
			double q_xb, q_xe;

			//solve on begin - x_b
			i = 0;//Nx_u[0];
			q_xb = 0;
			for(j = 0; j < Ny; j++)
			{
				node = i + j * Nx;
				if(vol_inf_fb_bool[node])
				{
					q_xb += rho[node]* u[node] * hy[j];
				}
			}


			//solve on exit - x_e
			i = Nx_u[1];
			q_xe = 0;
			for(j = 0; j < Ny; j++)
			{
				node = i + j * Nx;
				if(vol_inf_fb_bool[node])
				{
					q_xe += rho[node] * u[node] * hy[j];
				}
			}

			//Write data ti file DataInShortFormat.txt
			ofstream output_data;
			output_data.open("DataInShortFormat.txt");

			if(output_data.is_open())
			{
				int width_between_columns = 30;
				int Ndigits = 15;

				output_data
					<< std::setiosflags(std::ios::scientific)
					<< std::setprecision(Ndigits);

				output_data
					<< std::setw(width_between_columns) << Re
					<< std::setw(width_between_columns) << Kn
					<< std::setw(width_between_columns) << CD_p_x + CD_fr_x
					<< std::setw(width_between_columns) << CD_p_x
					<< std::setw(width_between_columns) << u_gas_nd
					<< std::setw(width_between_columns) << q_xb
					<< std::setw(width_between_columns) << q_xe
					<< std::setw(width_between_columns) << p_BC_xb
					<< std::setw(width_between_columns) << p_BC_xe
					<< endl;

				output_data.close();
			}
			else
			{
				cout << "The program can not open file DataInShortFormat.txt to write data." << endl;
			}


		}


		//write data for p_BC_xb nad fixed velocities on x_b to file
		if(Pressure_BC/* && (p_BC_xb_correction_method > 0)*/)
		{
			string p_and_u_xb_it_file_name = "p_and_u_xb_it.txt";

			ofstream output_data;

			//if we continue the old solve the file will not be erased
			if(ToCheckForFileIsNow_p_u_x_b)
			{
				if(ToReadSolvedDataFromFile)
					output_data.open(p_and_u_xb_it_file_name.c_str(), ios::out | ios::app);
				else
					output_data.open(p_and_u_xb_it_file_name.c_str());
			}
			else
				output_data.open(p_and_u_xb_it_file_name.c_str(), ios::out | ios::app);


			if(output_data.is_open())
			{
			int width_between_columns = 30;
				int Ndigits = 15;


				if(!ToReadSolvedDataFromFile && ToCheckForFileIsNow_p_u_x_b)
				{
					output_data << std::setw(width_between_columns) << "//it (number of step by time)";
					output_data << std::setw(width_between_columns) << "ht (step by time)";
					output_data << std::setw(width_between_columns) << "p_BC_xb";
					output_data << std::setw(width_between_columns) << "u_given_xb_in_it";
					output_data << std::setw(width_between_columns) << "u_given_xb_in_it_max";
					output_data << std::setw(width_between_columns) << "u_given_xb_in_it_mean";
					output_data << std::setw(width_between_columns) << "p_correction";

					output_data << endl;
				}


				output_data << std::setw(width_between_columns) << it;

				output_data
					<< std::setiosflags(std::ios::scientific)
					<< std::setprecision(Ndigits);

				output_data << std::setw(width_between_columns) << ht;

				output_data << std::setw(width_between_columns) << p_BC_xb;
				output_data << std::setw(width_between_columns) << u_given_xb_in_it;
				output_data << std::setw(width_between_columns) << u_given_xb_in_it_max;
				output_data << std::setw(width_between_columns) << u_given_xb_in_it_mean;
				output_data << std::setw(width_between_columns) << p_correction;

				output_data << endl;


				output_data.close();

				ToCheckForFileIsNow_p_u_x_b = false;
			}
			else
			{
				cout << "The program can not open file p_and_u_xb_it.txt to write data." << endl;
			}


		}

	}


}


void define_rule_vectors(void)
{
	//if rule vectors are defined they must be deleted
	if(is_rule_vectors_defined)
	{
		vol_inf_fb_bool = NULL; delete [] vol_inf_fb_bool;

		rule_vector_ij_nch_u = NULL; delete [] rule_vector_ij_nch_u;
		rule_vector_j_nch_u = NULL; delete [] rule_vector_j_nch_u;

		rule_vector_ij_bc_u = NULL; delete [] rule_vector_ij_bc_u;
		rule_vector_j_bc_u = NULL; delete [] rule_vector_j_bc_u;


		rule_vector_ij_nch_v = NULL; delete [] rule_vector_ij_nch_v;
		rule_vector_j_nch_v = NULL; delete [] rule_vector_j_nch_v;

		rule_vector_ij_bc_v = NULL; delete [] rule_vector_ij_bc_v;
		rule_vector_j_bc_v = NULL; delete [] rule_vector_j_bc_v;


		rule_vector_ij_nch_others = NULL; delete [] rule_vector_ij_nch_others;
		rule_vector_j_nch_others = NULL; delete [] rule_vector_j_nch_others;

		rule_vector_ij_bc_others = NULL; delete [] rule_vector_ij_bc_others;
		rule_vector_j_bc_others = NULL; delete [] rule_vector_j_bc_others;
	}

    vol_inf_fb_bool = new bool [Na];

	//define vol_inf_fb_bool
	for(counter = 0; counter < Na; counter++)
	{
		vol_inf_fb_bool[counter] = (vol_inf_fb[counter] == fluid);
	}


	//define vectors and variables for u
	//Define vectors and variables for area where is no check
	Na_nch_u = 0;
	for(j = Ny_u[0]; j < Ny_u[1]; j++)
	{
		for(i = Nx_u[0]; i < Nx_u[1]; i++)
		{
			node = i + j * Nx;

			if(Periodic_boundary_conditions_about_OX)
			{
				ij_function_PeriodicBC(i, j);
			}
			else
			{
				ij_function(node);
			}

			//this are all checks which is make in loop.
			//If in all we have fluid we do not need to check anything and no in edn
			//of domain.
			if(((Ny_u[0] + BlockLengthSwap_J1) < j) && (j < (Ny_u[1] - BlockLengthSwap_J_1))
			&& ((Nx_u[0] + BlockLengthSwap_I1) < i) && (i < (Nx_u[1] - BlockLengthSwap_I_1))
			&& vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[ij]
			&& vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[ij_1]
			&& vol_inf_fb_bool[i_1j1] && vol_inf_fb_bool[ij1]
			&& vol_inf_fb_bool[i1j] && vol_inf_fb_bool[i_2j])
			{
				Na_nch_u++;
			}

		}
	}

	//make rule vectors
	rule_vector_ij_nch_u = new unsigned int [Na_nch_u];
	rule_vector_j_nch_u = new unsigned int [Na_nch_u];

	//define rule vectors
	counter = 0;
	for(j = Ny_u[0]; j < Ny_u[1]; j++)
	{
		for(i = Nx_u[0]; i < Nx_u[1]; i++)
		{
			node = i + j * Nx;

			if(Periodic_boundary_conditions_about_OX)
			{
				ij_function_PeriodicBC(i, j);
			}
			else
			{
				ij_function(node);
			}

			//this are all checks which is make in loop.
			//If in all we have fluid we do not need to check anything and no in edn
			//of domain.
			if(((Ny_u[0] + BlockLengthSwap_J1) < j) && (j < (Ny_u[1] - BlockLengthSwap_J_1))
			&& ((Nx_u[0] + BlockLengthSwap_I1) < i) && (i < (Nx_u[1] - BlockLengthSwap_I_1))
			&& vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[ij]
			&& vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[ij_1]
			&& vol_inf_fb_bool[i_1j1] && vol_inf_fb_bool[ij1]
			&& vol_inf_fb_bool[i1j] && vol_inf_fb_bool[i_2j])
			{
				rule_vector_ij_nch_u[counter] = ij;
				rule_vector_j_nch_u[counter] = j;

				counter++;
			}

		}
	}


	//Define vectors and variables for area where are boundary conditions
	Na_bc_u = 0;
	for(j = Ny_u[0]; j < Ny_u[1]; j++)
	{
		for(i = Nx_u[0]; i < Nx_u[1]; i++)
		{
			node = i + j * Nx;

			if(Periodic_boundary_conditions_about_OX)
			{
				ij_function_PeriodicBC(i, j);
			}
			else
			{
				ij_function(node);
			}

			//this are all checks which is make in loop.
			//If in all we have body in one cell we need to make check.
			if(vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[ij]
			&&
			((j <= (Ny_u[0] + BlockLengthSwap_J1)) || ((Ny_u[1] - BlockLengthSwap_J_1) <= j)
			|| (i <= (Nx_u[0] + BlockLengthSwap_I1)) || ((Nx_u[1] - BlockLengthSwap_I_1) <= i)
			|| (!(vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[ij_1]
				&& vol_inf_fb_bool[i_1j1] && vol_inf_fb_bool[ij1]
				&& vol_inf_fb_bool[i1j] && vol_inf_fb_bool[i_2j]))))
			{
				Na_bc_u++;
			}

		}
	}

	//make rule vectors
	rule_vector_ij_bc_u = new unsigned int [Na_bc_u];
	rule_vector_j_bc_u = new unsigned int [Na_bc_u];

	//define rule vectors
	counter = 0;
	for(j = Ny_u[0]; j < Ny_u[1]; j++)
	{
		for(i = Nx_u[0]; i < Nx_u[1]; i++)
		{
			node = i + j * Nx;

			if(Periodic_boundary_conditions_about_OX)
			{
				ij_function_PeriodicBC(i, j);
			}
			else
			{
				ij_function(node);
			}

			//this are all checks which is make in loop.
			//If in all we have body in one cell we need to make check.
			if(vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[ij]
			&&
			((j <= (Ny_u[0] + BlockLengthSwap_J1)) || ((Ny_u[1] - BlockLengthSwap_J_1) <= j)
			|| (i <= (Nx_u[0] + BlockLengthSwap_I1)) || ((Nx_u[1] - BlockLengthSwap_I_1) <= i)
			|| (!(vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[ij_1]
				&& vol_inf_fb_bool[i_1j1] && vol_inf_fb_bool[ij1]
				&& vol_inf_fb_bool[i1j] && vol_inf_fb_bool[i_2j]))))
			{
				rule_vector_ij_bc_u[counter] = ij;
				rule_vector_j_bc_u[counter] = j;

				counter++;
			}

		}
	}


	//Define vectors and variables for v
	//Define vectors and variables for area where is no check
	Na_nch_v = 0;
	for(j = Ny_v[0]; j < Ny_v[1]; j++)
	{
		for(i = Nx_v[0]; i < Nx_v[1]; i++)
		{
			node = i + j * Nx;

			if(Periodic_boundary_conditions_about_OX)
			{
				ij_function_PeriodicBC(i, j);
			}
			else
			{
				ij_function(node);
			}

			//this are all checks which is make in loop.
			//If in all we have fluid we do not need to check anything.
			if(((Ny_u[0] + BlockLengthSwap_J1) < j) && (j < (Ny_u[1] - BlockLengthSwap_J_1))
			&& ((Nx_u[0] + BlockLengthSwap_I1) < i) && (i < (Nx_u[1] - BlockLengthSwap_I_1))
			&& vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij]
			&& vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[i_1j]
			&& vol_inf_fb_bool[i1j_1] && vol_inf_fb_bool[i1j]
			&& vol_inf_fb_bool[ij_1 - Nx * (0 <= (int)(j - 2))] && vol_inf_fb_bool[ij1])
			{
				Na_nch_v++;
			}

		}
	}

	//make rule vectors
	rule_vector_ij_nch_v = new unsigned int [Na_nch_v];
	rule_vector_j_nch_v = new unsigned int [Na_nch_v];

	//define rule vectors
	counter = 0;
	for(j = Ny_v[0]; j < Ny_v[1]; j++)
	{
		for(i = Nx_v[0]; i < Nx_v[1]; i++)
		{
			node = i + j * Nx;

			if(Periodic_boundary_conditions_about_OX)
			{
				ij_function_PeriodicBC(i, j);
			}
			else
			{
				ij_function(node);
			}

			//this are all checks which is make in loop.
			//If in all we have fluid we do not need to check anything.
			if(((Ny_u[0] + BlockLengthSwap_J1) < j) && (j < (Ny_u[1] - BlockLengthSwap_J_1))
			&& ((Nx_u[0] + BlockLengthSwap_I1) < i) && (i < (Nx_u[1] - BlockLengthSwap_I_1))
			&& vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij]
			&& vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[i_1j]
			&& vol_inf_fb_bool[i1j_1] && vol_inf_fb_bool[i1j]
			&& vol_inf_fb_bool[ij_1 - Nx * (0 <= (int)(j - 2))] && vol_inf_fb_bool[ij1])
			{
				rule_vector_ij_nch_v[counter] = ij;
				rule_vector_j_nch_v[counter] = j;

				counter++;
			}

		}
	}

	//Define vectors and variables for area where are boundary conditions
	Na_bc_v = 0;
	for(j = Ny_v[0]; j < Ny_v[1]; j++)
	{
		for(i = Nx_v[0]; i < Nx_v[1]; i++)
		{
			node = i + j * Nx;

			if(Periodic_boundary_conditions_about_OX)
			{
				ij_function_PeriodicBC(i, j);
			}
			else
			{
				ij_function(node);
			}

			if(vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij]
			&&
			((j <= (Ny_u[0] + BlockLengthSwap_J1)) || ((Ny_u[1] - BlockLengthSwap_J_1) <= j)
			|| (i <= (Nx_u[0] + BlockLengthSwap_I1)) || ((Nx_u[1] - BlockLengthSwap_I_1) <= i)
			|| (!(vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[i_1j]
				  && vol_inf_fb_bool[i1j_1] && vol_inf_fb_bool[i1j]
				  && vol_inf_fb_bool[ij_1 - Nx * (0 <= (int)(j - 2))] && vol_inf_fb_bool[ij1]))))
			{
				Na_bc_v++;
			}

		}
	}

	//make rule vectors
	rule_vector_ij_bc_v = new unsigned int [Na_bc_v];
	rule_vector_j_bc_v = new unsigned int [Na_bc_v];

	//define rule vectors
	counter = 0;
	for(j = Ny_v[0]; j < Ny_v[1]; j++)
	{
		for(i = Nx_v[0]; i < Nx_v[1]; i++)
		{
			node = i + j * Nx;

			if(Periodic_boundary_conditions_about_OX)
			{
				ij_function_PeriodicBC(i, j);
			}
			else
			{
				ij_function(node);
			}

			if(vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij]
			&&
			((j <= (Ny_u[0] + BlockLengthSwap_J1)) || ((Ny_u[1] - BlockLengthSwap_J_1) <= j)
			|| (i <= (Nx_u[0] + BlockLengthSwap_I1)) || ((Nx_u[1] - BlockLengthSwap_I_1) <= i)
			|| (!(vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[i_1j]
				  && vol_inf_fb_bool[i1j_1] && vol_inf_fb_bool[i1j]
				  && vol_inf_fb_bool[ij_1 - Nx * (0 <= (int)(j - 2))] && vol_inf_fb_bool[ij1]))))
			{
				rule_vector_ij_bc_v[counter] = ij;
				rule_vector_j_bc_v[counter] = j;

				counter++;
			}

		}
	}


	//Define vectors and variables for others
	//Define vectors and variables for area where is no check
	Na_nch_others = 0;
	for(j = Ny_others[0]; j < Ny_others[1]; j++)
	{
		for(i = Nx_others[0]; i < Nx_others[1]; i++)
		{
			node = i + j * Nx;

			if(Periodic_boundary_conditions_about_OX)
			{
				ij_function_PeriodicBC(i, j);
			}
			else
			{
				ij_function(node);
			}

			//thish are all checks which is make in loop.
			//If in all we have fluid we do not need to check anything and no in edn
			//of domain.
			if(((Ny_u[0] + BlockLengthSwap_J1) < j) && (j < (Ny_u[1] - BlockLengthSwap_J_1))
			&& ((Nx_u[0] + BlockLengthSwap_I1) < i) && (i < (Nx_u[1] - BlockLengthSwap_I_1))
			&& vol_inf_fb_bool[ij]
			&& vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[i1j]
			&& vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij1]
			&& vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[i_1j1]
			&& vol_inf_fb_bool[i1j_1] && vol_inf_fb_bool[i1j1])
			{
				Na_nch_others++;
			}

		}
	}

	//make rule vectors
	rule_vector_ij_nch_others = new unsigned int [Na_nch_others];
	rule_vector_j_nch_others = new unsigned int [Na_nch_others];

	//define rule vectors
	counter = 0;
	for(j = Ny_others[0]; j < Ny_others[1]; j++)
	{
		for(i = Nx_others[0]; i < Nx_others[1]; i++)
		{
			node = i + j * Nx;

			if(Periodic_boundary_conditions_about_OX)
			{
				ij_function_PeriodicBC(i, j);
			}
			else
			{
				ij_function(node);
			}

			//thish are all checks which is make in loop.
			//If in all we have fluid we do not need to check anything and no in edn
			//of domain.
			if(((Ny_u[0] + BlockLengthSwap_J1) < j) && (j < (Ny_u[1] - BlockLengthSwap_J_1))
			&& ((Nx_u[0] + BlockLengthSwap_I1) < i) && (i < (Nx_u[1] - BlockLengthSwap_I_1))
			&& vol_inf_fb_bool[ij]
			&& vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[i1j]
			&& vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij1]
			&& vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[i_1j1]
			&& vol_inf_fb_bool[i1j_1] && vol_inf_fb_bool[i1j1])
			{
				rule_vector_ij_nch_others[counter] = ij;
				rule_vector_j_nch_others[counter] = j;

				counter++;
			}

		}
	}


	//Define vectors and variables for area where are boundary conditions
	Na_bc_others = 0;
	for(j = Ny_others[0]; j < Ny_others[1]; j++)
	{
		for(i = Nx_others[0]; i < Nx_others[1]; i++)
		{
			node = i + j * Nx;

			if(Periodic_boundary_conditions_about_OX)
			{
				ij_function_PeriodicBC(i, j);
			}
			else
			{
				ij_function(node);
			}


			//this are all checks which is make in loop.
			//If in all we have body in one cell we need to make check.
			if(vol_inf_fb_bool[ij]
			&&
			((j <= (Ny_u[0] + BlockLengthSwap_J1)) || ((Ny_u[1] - BlockLengthSwap_J_1) <= j)
			|| (i <= (Nx_u[0] + BlockLengthSwap_I1)) || ((Nx_u[1] - BlockLengthSwap_I_1) <= i)
			|| (!(vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[i1j]
				  && vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij1]
				  && vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[i_1j1]
				  && vol_inf_fb_bool[i1j_1] && vol_inf_fb_bool[i1j1]))))
			{
				Na_bc_others++;
			}

		}
	}

	//make rule vectors
	rule_vector_ij_bc_others = new unsigned int [Na_bc_others];
	rule_vector_j_bc_others = new unsigned int [Na_bc_others];

	//define rule vectors
	counter = 0;
	for(j = Ny_others[0]; j < Ny_others[1]; j++)
	{
		for(i = Nx_others[0]; i < Nx_others[1]; i++)
		{
			node = i + j * Nx;

			if(Periodic_boundary_conditions_about_OX)
			{
				ij_function_PeriodicBC(i, j);
			}
			else
			{
				ij_function(node);
			}


			//this are all checks which is make in loop.
			//If in all we have body in one cell we need to make check.
			if(vol_inf_fb_bool[ij]
			&&
			((j <= (Ny_u[0] + BlockLengthSwap_J1)) || ((Ny_u[1] - BlockLengthSwap_J_1) <= j)
			|| (i <= (Nx_u[0] + BlockLengthSwap_I1)) || ((Nx_u[1] - BlockLengthSwap_I_1) <= i)
			|| (!(vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[i1j]
				  && vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij1]
				  && vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[i_1j1]
				  && vol_inf_fb_bool[i1j_1] && vol_inf_fb_bool[i1j1]))))
			{
				rule_vector_ij_bc_others[counter] = ij;
				rule_vector_j_bc_others[counter] = j;

				counter++;
			}

		}
	}



//OLD code
//	if(Periodic_boundary_conditions_about_OX)
//	{
//		//define vectors and variables for u
//		//Define vectors and variables for area where is no check
//		Na_nch_u = 0;
//		for(j = Ny_u[0]; j < Ny_u[1]; j++)
//		{
//			for(i = Nx_u[0]; i < Nx_u[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//				//this are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything and no in edn
//				//of domain.
//				if((0 < i) && (i < Nx - 2)
//				&& vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[ij]
//				&& vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[ij_1]
//				&& vol_inf_fb_bool[i_1j1] && vol_inf_fb_bool[ij1])
//				{
//					Na_nch_u++;
//				}
//
//			}
//		}
//
//		//make rule vectors
//		rule_vector_ij_nch_u = new unsigned int [Na_nch_u];
//		rule_vector_j_nch_u = new unsigned int [Na_nch_u];
//
//		//define rule vectors
//		counter = 0;
//		for(j = Ny_u[0]; j < Ny_u[1]; j++)
//		{
//			for(i = Nx_u[0]; i < Nx_u[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//				//thish are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything.
//				if((0 < i) && (i < Nx - 2)
//				&& vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[ij]
//				&& vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[ij_1]
//				&& vol_inf_fb_bool[i_1j1] && vol_inf_fb_bool[ij1])
//				{
//					rule_vector_ij_nch_u[counter] = ij;
//					rule_vector_j_nch_u[counter] = j;
//
//					counter++;
//				}
//
//			}
//		}
//
//
//		//Define vectors and variables for area where are boundary conditions
//		Na_bc_u = 0;
//		for(j = Ny_u[0]; j < Ny_u[1]; j++)
//		{
//			for(i = Nx_u[0]; i < Nx_u[1]; i++)
//			{
//				ij_function_PeriodicBC(i, j);
//
//
//				//this are all checks which is make in loop.
//				//If in all we have body in one cell we need to make check.
//				if(vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[ij]
//				&&
//				((i == 0) || (Nx - 2 <= i)
//				|| (!(vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[ij_1]
//					&& vol_inf_fb_bool[i_1j1] && vol_inf_fb_bool[ij1]))))
//				{
//					Na_bc_u++;
//				}
//
//			}
//		}
//
//		//make rule vectors
//		rule_vector_ij_bc_u = new unsigned int [Na_bc_u];
//		rule_vector_j_bc_u = new unsigned int [Na_bc_u];
//
//		//define rule vectors
//		counter = 0;
//		for(j = Ny_u[0]; j < Ny_u[1]; j++)
//		{
//			for(i = Nx_u[0]; i < Nx_u[1]; i++)
//			{
//				ij_function_PeriodicBC(i, j);
//
//
//				//this are all checks which is make in loop.
//				//If in all we have body in one cell we need to make check.
//				if(vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[ij]
//				&&
//				((i == 0) || (Nx - 2 <= i)
//				|| (!(vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[ij_1]
//					&& vol_inf_fb_bool[i_1j1] && vol_inf_fb_bool[ij1]))))
//				{
//					rule_vector_ij_bc_u[counter] = ij;
//					rule_vector_j_bc_u[counter] = j;
//
//					counter++;
//				}
//
//			}
//		}
//
//
//		//Define vectors and variables for v
//		//Define vectors and variables for area where is no check
//		Na_nch_v = 0;
//		for(j = Ny_v[0]; j < Ny_v[1]; j++)
//		{
//			for(i = Nx_v[0]; i < Nx_v[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//				//this are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything.
//				if((0 < i) && (i < Nx - 2)
//				&& vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij]
//				&& vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[i_1j]
//				&& vol_inf_fb_bool[i1j_1] && vol_inf_fb_bool[i1j]
//				&& vol_inf_fb_bool[ij_1 - Nx * (0 <= (int)(j - 2))] && vol_inf_fb_bool[ij1])
//				{
//					Na_nch_v++;
//				}
//
//			}
//		}
//
//		//make rule vectors
//		rule_vector_ij_nch_v = new unsigned int [Na_nch_v];
//		rule_vector_j_nch_v = new unsigned int [Na_nch_v];
//
//		//define rule vectors
//		counter = 0;
//		for(j = Ny_v[0]; j < Ny_v[1]; j++)
//		{
//			for(i = Nx_v[0]; i < Nx_v[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//				//this are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything.
//				if((0 < i) && (i < Nx - 2)
//				&& vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij]
//				&& vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[i_1j]
//				&& vol_inf_fb_bool[i1j_1] && vol_inf_fb_bool[i1j]
//				&& vol_inf_fb_bool[ij_1 - Nx * (0 <= (int)(j - 2))] && vol_inf_fb_bool[ij1])
//				{
//					rule_vector_ij_nch_v[counter] = ij;
//					rule_vector_j_nch_v[counter] = j;
//
//					counter++;
//				}
//
//			}
//		}
//
//		//Define vectors and variables for area where are boundary conditions
//		Na_bc_v = 0;
//		for(j = Ny_v[0]; j < Ny_v[1]; j++)
//		{
//			for(i = Nx_v[0]; i < Nx_v[1]; i++)
//			{
//				ij_function_PeriodicBC(i, j);
//
//
//				//this are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything.
//				if(vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij]
//				&&
//				((i == 0) || (Nx - 2 <= i)
//				|| (!(vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[i_1j]
//					  && vol_inf_fb_bool[i1j_1] && vol_inf_fb_bool[i1j]
//					  && vol_inf_fb_bool[ij_1 - Nx * (0 <= (int)(j - 2))] && vol_inf_fb_bool[ij1]))))
//				{
//					Na_bc_v++;
//				}
//
//			}
//		}
//
//		//make rule vectors
//		rule_vector_ij_bc_v = new unsigned int [Na_bc_v];
//		rule_vector_j_bc_v = new unsigned int [Na_bc_v];
//
//		//define rule vectors
//		counter = 0;
//		for(j = Ny_v[0]; j < Ny_v[1]; j++)
//		{
//			for(i = Nx_v[0]; i < Nx_v[1]; i++)
//			{
//				ij_function_PeriodicBC(i, j);
//
//
//				//thish are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything.
//				if(vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij]
//				&&
//				((i == 0) || (Nx - 2 <= i)
//				|| (!(vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[i_1j]
//					  && vol_inf_fb_bool[i1j_1] && vol_inf_fb_bool[i1j]
//					  && vol_inf_fb_bool[ij_1 - Nx * (0 <= (int)(j - 2))] && vol_inf_fb_bool[ij1]))))
//				{
//					rule_vector_ij_bc_v[counter] = ij;
//					rule_vector_j_bc_v[counter] = j;
//
//					counter++;
//				}
//
//			}
//		}
//
//
//		//Define vectors and variables for others
//		//Define vectors and variables for area where is no check
//		Na_nch_others = 0;
//		for(j = Ny_others[0]; j < Ny_others[1]; j++)
//		{
//			for(i = Nx_others[0]; i < Nx_others[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//				//thish are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything and no in edn
//				//of domain.
//				if((0 < i) && (i < Nx - 2)
//				&& vol_inf_fb_bool[ij]
//				&& vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[i1j]
//				&& vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij1])
//				{
//					Na_nch_others++;
//				}
//
//			}
//		}
//
//		//make rule vectors
//		rule_vector_ij_nch_others = new unsigned int [Na_nch_others];
//		rule_vector_j_nch_others = new unsigned int [Na_nch_others];
//
//		//define rule vectors
//		counter = 0;
//		for(j = Ny_others[0]; j < Ny_others[1]; j++)
//		{
//			for(i = Nx_others[0]; i < Nx_others[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//				//thish are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything and no in edn
//				//of domain.
//				if((0 < i) && (i < Nx - 2)
//				&& vol_inf_fb_bool[ij]
//				&& vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[i1j]
//				&& vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij1])
//				{
//					rule_vector_ij_nch_others[counter] = ij;
//					rule_vector_j_nch_others[counter] = j;
//
//					counter++;
//				}
//
//			}
//		}
//
//
//		//Define vectors and vatriables for area where are boundary conditions
//		Na_bc_others = 0;
//		for(j = Ny_others[0]; j < Ny_others[1]; j++)
//		{
//			for(i = Nx_others[0]; i < Nx_others[1]; i++)
//			{
//				ij_function_PeriodicBC(i, j);
//
//
//				//thish are all checks which is make in loop.
//				//If in all we have body in one cell we need to make check.
//				if(vol_inf_fb_bool[ij]
//				&&
//				((i == 0) || (Nx - 2 <= i)
//				|| (!(vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[i1j]
//					&& vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij1]))))
//				{
//					Na_bc_others++;
//				}
//
//			}
//		}
//
//		//make rule vectors
//		rule_vector_ij_bc_others = new unsigned int [Na_bc_others];
//		rule_vector_j_bc_others = new unsigned int [Na_bc_others];
//
//		//define rule vectors
//		counter = 0;
//		for(j = Ny_others[0]; j < Ny_others[1]; j++)
//		{
//			for(i = Nx_others[0]; i < Nx_others[1]; i++)
//			{
//				ij_function_PeriodicBC(i, j);
//
//
//				//thish are all checks which is make in loop.
//				//If in all we have body in one cell we need to make check.
//				if(vol_inf_fb_bool[ij]
//				&&
//				((i == 0) || (Nx - 2 <= i)
//				|| (!(vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[i1j]
//					&& vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij1]))))
//				{
//					rule_vector_ij_bc_others[counter] = ij;
//					rule_vector_j_bc_others[counter] = j;
//
//					counter++;
//				}
//
//			}
//		}
//
//
//	}
//	else if(Pressure_BC)
//	{
//		//define vectors and variables for u
//		//Define vectors and variables for area where is no check
//		Na_nch_u = 0;
//		for(j = Ny_u[0]; j < Ny_u[1]; j++)
//		{
//			for(i = Nx_u[0]; i < Nx_u[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//				//this are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything and no in edn
//				//of domain.
//				if((Nx_u[0] < i) && (i < (Nx_u[1] - BlockLengthSwap_I_1))
//				&& (Ny_u[0] < j) && (j < (Ny_u[1] - BlockLengthSwap_J_1))
//
//				&& vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[ij]
//				&& vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[ij_1]
//				&& vol_inf_fb_bool[i_1j1] && vol_inf_fb_bool[ij1])
//				{
//					Na_nch_u++;
//				}
//
//			}
//		}
//
//		//make rule vectors
//		rule_vector_ij_nch_u = new unsigned int [Na_nch_u];
//		rule_vector_j_nch_u = new unsigned int [Na_nch_u];
//
//		//define rule vectors
//		counter = 0;
//		for(j = Ny_u[0]; j < Ny_u[1]; j++)
//		{
//			for(i = Nx_u[0]; i < Nx_u[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//				//thish are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything.
//				if((Nx_u[0] < i) && (i < (Nx_u[1] - BlockLengthSwap_I_1))
//				&& (Ny_u[0] < j) && (j < (Ny_u[1] - BlockLengthSwap_J_1))
//
//				&& vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[ij]
//				&& vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[ij_1]
//				&& vol_inf_fb_bool[i_1j1] && vol_inf_fb_bool[ij1])
//				{
//					rule_vector_ij_nch_u[counter] = ij;
//					rule_vector_j_nch_u[counter] = j;
//
//					counter++;
//				}
//
//			}
//		}
//
//
//		//Define vectors and variables for area where are boundary conditions
//		Na_bc_u = 0;
//		for(j = Ny_u[0]; j < Ny_u[1]; j++)
//		{
//			for(i = Nx_u[0]; i < Nx_u[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//
//				//this are all checks which is make in loop.
//				//If in all we have body in one cell we need to make check.
//				if(vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[ij]
//				&&
//				//Include all area boundaries
//				((i <= Nx_u[0]) || ((Nx_u[1] - BlockLengthSwap_I_1) <= i)
//				|| (j <= Ny_u[0]) || ((Ny_u[1] - BlockLengthSwap_J_1) <= j)
//
//				|| (!(vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[ij_1]
//					&& vol_inf_fb_bool[i_1j1] && vol_inf_fb_bool[ij1]))))
//				{
//					Na_bc_u++;
//				}
//
//			}
//		}
//
//		//make rule vectors
//		rule_vector_ij_bc_u = new unsigned int [Na_bc_u];
//		rule_vector_j_bc_u = new unsigned int [Na_bc_u];
//
//		//define rule vectors
//		counter = 0;
//		for(j = Ny_u[0]; j < Ny_u[1]; j++)
//		{
//			for(i = Nx_u[0]; i < Nx_u[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//
//				//this are all checks which is make in loop.
//				//If in all we have body in one cell we need to make check.
//				if(vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[ij]
//				&&
//				//Include all area boundaries
//				((i <= Nx_u[0]) || ((Nx_u[1] - BlockLengthSwap_I_1) <= i)
//				|| (j <= Ny_u[0]) || ((Ny_u[1] - BlockLengthSwap_J_1) <= j)
//
//				|| (!(vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[ij_1]
//					&& vol_inf_fb_bool[i_1j1] && vol_inf_fb_bool[ij1]))))
//				{
//					rule_vector_ij_bc_u[counter] = ij;
//					rule_vector_j_bc_u[counter] = j;
//
//					counter++;
//				}
//
//			}
//		}
//
//
//		//define vectors and variables for v
//		//Define vectors and variables for area where is no check
//		Na_nch_v = 0;
//		for(j = Ny_v[0]; j < Ny_v[1]; j++)
//		{
//			for(i = Nx_v[0]; i < Nx_v[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//				//this are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything.
//				if((Nx_v[0] < i) && (i < (Nx_v[1] - BlockLengthSwap_I_1))
//				&& (Ny_v[0] < j) && (j < (Ny_v[1] - BlockLengthSwap_J_1))
//
//				&& vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij]
//				&& vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[i_1j]
//				&& vol_inf_fb_bool[i1j_1] && vol_inf_fb_bool[i1j])
//				{
//					Na_nch_v++;
//				}
//
//			}
//		}
//
//		//make rule vectors
//		rule_vector_ij_nch_v = new unsigned int [Na_nch_v];
//		rule_vector_j_nch_v = new unsigned int [Na_nch_v];
//
//		//define rule vectors
//		counter = 0;
//		for(j = Ny_v[0]; j < Ny_v[1]; j++)
//		{
//			for(i = Nx_v[0]; i < Nx_v[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//				//thish are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything.
//				if((Nx_v[0] < i) && (i < (Nx_v[1] - BlockLengthSwap_I_1))
//				&& (Ny_v[0] < j) && (j < (Ny_v[1] - BlockLengthSwap_J_1))
//
//				&& vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij]
//				&& vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[i_1j]
//				&& vol_inf_fb_bool[i1j_1] && vol_inf_fb_bool[i1j])
//				{
//					rule_vector_ij_nch_v[counter] = ij;
//					rule_vector_j_nch_v[counter] = j;
//
//					counter++;
//				}
//
//			}
//		}
//
//		//Define vectors and vatriables for area where are boundary conditions
//		Na_bc_v = 0;
//		for(j = Ny_v[0]; j < Ny_v[1]; j++)
//		{
//			for(i = Nx_v[0]; i < Nx_v[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//
//				//thish are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything.
//				if(vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij]
//				&&
//				//Include all area boundaries
//				((i <= Nx_v[0]) || ((Nx_v[1] - BlockLengthSwap_I_1) <= i)
//				|| (j <= Ny_v[0]) || ((Ny_v[1] - BlockLengthSwap_J_1) <= j)
//
//				|| (!(vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[i_1j]
//					&& vol_inf_fb_bool[i1j_1] && vol_inf_fb_bool[i1j]))))
//				{
//					Na_bc_v++;
//				}
//
//			}
//		}
//
//		//make rule vectors
//		rule_vector_ij_bc_v = new unsigned int [Na_bc_v];
//		rule_vector_j_bc_v = new unsigned int [Na_bc_v];
//
//		//define rule vectors
//		counter = 0;
//		for(j = Ny_v[0]; j < Ny_v[1]; j++)
//		{
//			for(i = Nx_v[0]; i < Nx_v[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//
//				//thish are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything.
//				if(vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij]
//				&&
//				//Include all area boundaries
//				((i <= Nx_v[0]) || ((Nx_v[1] - BlockLengthSwap_I_1) <= i)
//				|| (j <= Ny_v[0]) || ((Ny_v[1] - BlockLengthSwap_J_1) <= j)
//
//				|| (!(vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[i_1j]
//					&& vol_inf_fb_bool[i1j_1] && vol_inf_fb_bool[i1j]))))
//				{
//					rule_vector_ij_bc_v[counter] = ij;
//					rule_vector_j_bc_v[counter] = j;
//
//					counter++;
//				}
//
//			}
//		}
//
//
//		//define vectors and variables for others
//		//Define vectors and vatriables for area where is no check
//		Na_nch_others = 0;
//		for(j = Ny_others[0]; j < Ny_others[1]; j++)
//		{
//			for(i = Nx_others[0]; i < Nx_others[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//				//thish are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything and no in edn
//				//of domain.
//				if((Nx_others[0] < i) && (i < (Nx_others[1] - BlockLengthSwap_I_1))
//				&& (Ny_others[0] < j) && (j < (Ny_others[1] - BlockLengthSwap_J_1))
//
//				&& vol_inf_fb_bool[ij]
//				&& vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[i1j]
//				&& vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij1])
//				{
//					Na_nch_others++;
//				}
//
//			}
//		}
//
//		//make rule vectors
//		rule_vector_ij_nch_others = new unsigned int [Na_nch_others];
//		rule_vector_j_nch_others = new unsigned int [Na_nch_others];
//
//		//define rule vectors
//		counter = 0;
//		for(j = Ny_others[0]; j < Ny_others[1]; j++)
//		{
//			for(i = Nx_others[0]; i < Nx_others[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//				//thish are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything and no in edn
//				//of domain.
//				if((Nx_others[0] < i) && (i < (Nx_others[1] - BlockLengthSwap_I_1))
//				&& (Ny_others[0] < j) && (j < (Ny_others[1] - BlockLengthSwap_J_1))
//
//				&& vol_inf_fb_bool[ij]
//				&& vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[i1j]
//				&& vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij1])
//				{
//					rule_vector_ij_nch_others[counter] = ij;
//					rule_vector_j_nch_others[counter] = j;
//
//					counter++;
//				}
//
//			}
//		}
//
//
//		//Define vectors and vatriables for area where are boundary conditions
//		Na_bc_others = 0;
//		for(j = Ny_others[0]; j < Ny_others[1]; j++)
//		{
//			for(i = Nx_others[0]; i < Nx_others[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//
//				//thish are all checks which is make in loop.
//				//If in all we have body in one cell we need to make check.
//				if(vol_inf_fb_bool[ij]
//				&&
//				//Include all area boundaries
//				((i <= Nx_others[0]) || ((Nx_others[1] - BlockLengthSwap_I_1) <= i)
//				|| (j <= Ny_others[0]) || ((Ny_others[1] - BlockLengthSwap_J_1) <= j)
//
//				|| (!(vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[i1j]
//					&& vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij1]))))
//				{
//					Na_bc_others++;
//				}
//
//			}
//		}
//
//		//make rule vectors
//		rule_vector_ij_bc_others = new unsigned int [Na_bc_others];
//		rule_vector_j_bc_others = new unsigned int [Na_bc_others];
//
//		//define rule vectors
//		counter = 0;
//		for(j = Ny_others[0]; j < Ny_others[1]; j++)
//		{
//			for(i = Nx_others[0]; i < Nx_others[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//
//				//thish are all checks which is make in loop.
//				//If in all we have body in one cell we need to make check.
//				if(vol_inf_fb_bool[ij]
//				&&
//				//Include all area boundaries
//				((i <= Nx_others[0]) || ((Nx_others[1] - BlockLengthSwap_I_1) <= i)
//				|| (j <= Ny_others[0]) || ((Ny_others[1] - BlockLengthSwap_J_1) <= j)
//
//				|| (!(vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[i1j]
//					&& vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij1]))))
//				{
//					rule_vector_ij_bc_others[counter] = ij;
//					rule_vector_j_bc_others[counter] = j;
//
//					counter++;
//				}
//
//			}
//		}
//
//
//	}
//	else
//	{
//		//define vectors and variables for u
//		//Define vectors and variables for area where is no check
//		Na_nch_u = 0;
//		for(j = Ny_u[0]; j < Ny_u[1]; j++)
//		{
//			for(i = Nx_u[0]; i < Nx_u[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//				//this are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything and no in edn
//				//of domain.
//				if((Nx_u[0] < i) && (i < (Nx_u[1] - BlockLengthSwap_I_1))
//				&& (Ny_u[0] < j) && (j < (Ny_u[1] - BlockLengthSwap_J_1))
//
//				&& vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[ij]
//				&& vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[ij_1]
//				&& vol_inf_fb_bool[i_1j1] && vol_inf_fb_bool[ij1])
//				{
//					Na_nch_u++;
//				}
//
//			}
//		}
//
//		//make rule vectors
//		rule_vector_ij_nch_u = new unsigned int [Na_nch_u];
//		rule_vector_j_nch_u = new unsigned int [Na_nch_u];
//
//		//define rule vectors
//		counter = 0;
//		for(j = Ny_u[0]; j < Ny_u[1]; j++)
//		{
//			for(i = Nx_u[0]; i < Nx_u[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//				//thish are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything.
//				if((Nx_u[0] < i) && (i < (Nx_u[1] - BlockLengthSwap_I_1))
//				&& (Ny_u[0] < j) && (j < (Ny_u[1] - BlockLengthSwap_J_1))
//
//				&& vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[ij]
//				&& vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[ij_1]
//				&& vol_inf_fb_bool[i_1j1] && vol_inf_fb_bool[ij1])
//				{
//					rule_vector_ij_nch_u[counter] = ij;
//					rule_vector_j_nch_u[counter] = j;
//
//					counter++;
//				}
//
//			}
//		}
//
//
//		//Define vectors and vatriables for area where are boundary conditions
//		Na_bc_u = 0;
//		for(j = Ny_u[0]; j < Ny_u[1]; j++)
//		{
//			for(i = Nx_u[0]; i < Nx_u[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//
//				//thish are all checks which is make in loop.
//				//If in all we have body in one cell we need to make check.
//				if(vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[ij]
//				&&
//				//Include all area boundaries
//				((i <= Nx_u[0]) || ((Nx_u[1] - BlockLengthSwap_I_1) <= i)
//				|| (j <= Ny_u[0]) || ((Ny_u[1] - BlockLengthSwap_J_1) <= j)
//
//				|| (!(vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[ij_1]
//					&& vol_inf_fb_bool[i_1j1] && vol_inf_fb_bool[ij1]))))
//				{
//					Na_bc_u++;
//				}
//
//			}
//		}
//
//		//make rule vectors
//		rule_vector_ij_bc_u = new unsigned int [Na_bc_u];
//		rule_vector_j_bc_u = new unsigned int [Na_bc_u];
//
//		//define rule vectors
//		counter = 0;
//		for(j = Ny_u[0]; j < Ny_u[1]; j++)
//		{
//			for(i = Nx_u[0]; i < Nx_u[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//
//				//thish are all checks which is make in loop.
//				//If in all we have body in one cell we need to make check.
//				if(vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[ij]
//				&&
//				//Include all area boundaries
//				((i <= Nx_u[0]) || ((Nx_u[1] - BlockLengthSwap_I_1) <= i)
//				|| (j <= Ny_u[0]) || ((Ny_u[1] - BlockLengthSwap_J_1) <= j)
//
//				|| (!(vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[ij_1]
//					&& vol_inf_fb_bool[i_1j1] && vol_inf_fb_bool[ij1]))))
//				{
//					rule_vector_ij_bc_u[counter] = ij;
//					rule_vector_j_bc_u[counter] = j;
//
//					counter++;
//				}
//
//			}
//		}
//
//
//		//define vectors and variables for v
//		//Define vectors and vatriables for area where is no check
//		Na_nch_v = 0;
//		for(j = Ny_v[0]; j < Ny_v[1]; j++)
//		{
//			for(i = Nx_v[0]; i < Nx_v[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//				//thish are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything.
//				if((Nx_v[0] < i) && (i < (Nx_v[1] - BlockLengthSwap_I_1))
//				&& (Ny_v[0] < j) && (j < (Ny_v[1] - BlockLengthSwap_J_1))
//
//				&& vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij]
//				&& vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[i_1j]
//				&& vol_inf_fb_bool[i1j_1] && vol_inf_fb_bool[i1j])
//				{
//					Na_nch_v++;
//				}
//
//			}
//		}
//
//		//make rule vectors
//		rule_vector_ij_nch_v = new unsigned int [Na_nch_v];
//		rule_vector_j_nch_v = new unsigned int [Na_nch_v];
//
//		//define rule vectors
//		counter = 0;
//		for(j = Ny_v[0]; j < Ny_v[1]; j++)
//		{
//			for(i = Nx_v[0]; i < Nx_v[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//				//thish are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything.
//				if((Nx_v[0] < i) && (i < (Nx_v[1] - BlockLengthSwap_I_1))
//				&& (Ny_v[0] < j) && (j < (Ny_v[1] - BlockLengthSwap_J_1))
//
//				&& vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij]
//				&& vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[i_1j]
//				&& vol_inf_fb_bool[i1j_1] && vol_inf_fb_bool[i1j])
//				{
//					rule_vector_ij_nch_v[counter] = ij;
//					rule_vector_j_nch_v[counter] = j;
//
//					counter++;
//				}
//
//			}
//		}
//
//		//Define vectors and vatriables for area where are boundary conditions
//		Na_bc_v = 0;
//		for(j = Ny_v[0]; j < Ny_v[1]; j++)
//		{
//			for(i = Nx_v[0]; i < Nx_v[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//
//				//thish are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything.
//				if(vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij]
//				&&
//				//Include all area boundaries
//				((i <= Nx_v[0]) || ((Nx_v[1] - BlockLengthSwap_I_1) <= i)
//				|| (j <= Ny_v[0]) || ((Ny_v[1] - BlockLengthSwap_J_1) <= j)
//
//				|| (!(vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[i_1j]
//					&& vol_inf_fb_bool[i1j_1] && vol_inf_fb_bool[i1j]))))
//				{
//					Na_bc_v++;
//				}
//
//			}
//		}
//
//		//make rule vectors
//		rule_vector_ij_bc_v = new unsigned int [Na_bc_v];
//		rule_vector_j_bc_v = new unsigned int [Na_bc_v];
//
//		//define rule vectors
//		counter = 0;
//		for(j = Ny_v[0]; j < Ny_v[1]; j++)
//		{
//			for(i = Nx_v[0]; i < Nx_v[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//
//				//thish are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything.
//				if(vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij]
//				&&
//				//Include all area boundaries
//				((i <= Nx_v[0]) || ((Nx_v[1] - BlockLengthSwap_I_1) <= i)
//				|| (j <= Ny_v[0]) || ((Ny_v[1] - BlockLengthSwap_J_1) <= j)
//
//				|| (!(vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[i_1j]
//					&& vol_inf_fb_bool[i1j_1] && vol_inf_fb_bool[i1j]))))
//				{
//					rule_vector_ij_bc_v[counter] = ij;
//					rule_vector_j_bc_v[counter] = j;
//
//					counter++;
//				}
//
//			}
//		}
//
//
//		//define vectors and variables for others
//		//Define vectors and vatriables for area where is no check
//		Na_nch_others = 0;
//		for(j = Ny_others[0]; j < Ny_others[1]; j++)
//		{
//			for(i = Nx_others[0]; i < Nx_others[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//				//thish are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything and no in edn
//				//of domain.
//				if((Nx_others[0] < i) && (i < (Nx_others[1] - BlockLengthSwap_I_1))
//				&& (Ny_others[0] < j) && (j < (Ny_others[1] - BlockLengthSwap_J_1))
//
//				&& vol_inf_fb_bool[ij]
//				&& vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[i1j]
//				&& vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij1])
//				{
//					Na_nch_others++;
//				}
//
//			}
//		}
//
//		//make rule vectors
//		rule_vector_ij_nch_others = new unsigned int [Na_nch_others];
//		rule_vector_j_nch_others = new unsigned int [Na_nch_others];
//
//		//define rule vectors
//		counter = 0;
//		for(j = Ny_others[0]; j < Ny_others[1]; j++)
//		{
//			for(i = Nx_others[0]; i < Nx_others[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//				//thish are all checks which is make in loop.
//				//If in all we have fluid we do not need to check anything and no in edn
//				//of domain.
//				if((Nx_others[0] < i) && (i < (Nx_others[1] - BlockLengthSwap_I_1))
//				&& (Ny_others[0] < j) && (j < (Ny_others[1] - BlockLengthSwap_J_1))
//
//				&& vol_inf_fb_bool[ij]
//				&& vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[i1j]
//				&& vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij1])
//				{
//					rule_vector_ij_nch_others[counter] = ij;
//					rule_vector_j_nch_others[counter] = j;
//
//					counter++;
//				}
//
//			}
//		}
//
//
//		//Define vectors and vatriables for area where are boundary conditions
//		Na_bc_others = 0;
//		for(j = Ny_others[0]; j < Ny_others[1]; j++)
//		{
//			for(i = Nx_others[0]; i < Nx_others[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//
//				//thish are all checks which is make in loop.
//				//If in all we have body in one cell we need to make check.
//				if(vol_inf_fb_bool[ij]
//				&&
//				//Include all area boundaries
//				((i <= Nx_others[0]) || ((Nx_others[1] - BlockLengthSwap_I_1) <= i)
//				|| (j <= Ny_others[0]) || ((Ny_others[1] - BlockLengthSwap_J_1) <= j)
//
//				|| (!(vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[i1j]
//					&& vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij1]))))
//				{
//					Na_bc_others++;
//				}
//
//			}
//		}
//
//		//make rule vectors
//		rule_vector_ij_bc_others = new unsigned int [Na_bc_others];
//		rule_vector_j_bc_others = new unsigned int [Na_bc_others];
//
//		//define rule vectors
//		counter = 0;
//		for(j = Ny_others[0]; j < Ny_others[1]; j++)
//		{
//			for(i = Nx_others[0]; i < Nx_others[1]; i++)
//			{
//				node = i + j * Nx;
//
//				ij_function(node);
//
//
//				//thish are all checks which is make in loop.
//				//If in all we have body in one cell we need to make check.
//				if(vol_inf_fb_bool[ij]
//				&&
//				//Include all area boundaries
//				((i <= Nx_others[0]) || ((Nx_others[1] - BlockLengthSwap_I_1) <= i)
//				|| (j <= Ny_others[0]) || ((Ny_others[1] - BlockLengthSwap_J_1) <= j)
//
//				|| (!(vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[i1j]
//					&& vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij1]))))
//				{
//					rule_vector_ij_bc_others[counter] = ij;
//					rule_vector_j_bc_others[counter] = j;
//
//					counter++;
//				}
//
//			}
//		}
//
//
//	}



	is_rule_vectors_defined = true;

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

	//Define the rule vectors for velocity slip boundary conditions for u and v
	//Define rule vector for velocity slip boundary conditions for u
	if(is_rule_vectors_BoundaryConditions_Velocities_SlipBC_u_defined)
	{
		delete [] rule_vector_ij_BoundaryConditions_Velocities_SlipBC_u;
		delete [] rule_vector_j_BoundaryConditions_Velocities_SlipBC_u;
	}

	//Initiate the number of elements of rule_vector_ij_BoundaryConditions_Velocities_SlipBC_u
	Na_BoundaryConditions_Velocities_SlipBC_u = 0;
	for(j = 0; j < Ny; j++)
	{
		for(i = 0; i < Nx; i++)
		{
			node = i + j * Nx;

			if(Periodic_boundary_conditions_about_OX)
			{
				ij_function_PeriodicBC_i_and_j_to_the_boundaries(i, j);
			}
			else
			{
				ij_function_i_and_j_to_the_boundaries(node);
			}

			if((fluid < vol_inf_fb[ij])
			&& ((fluid == vol_inf_fb[ij_1]) || (fluid == vol_inf_fb[ij1])))
			{
				Na_BoundaryConditions_Velocities_SlipBC_u++;
			}

		}
	}

	//Make rule vectors
	rule_vector_ij_BoundaryConditions_Velocities_SlipBC_u = new int [Na_BoundaryConditions_Velocities_SlipBC_u];
	rule_vector_j_BoundaryConditions_Velocities_SlipBC_u = new int [Na_BoundaryConditions_Velocities_SlipBC_u];

	//define elements in rule vector
	counter = 0;
	for(j = 0; j < Ny; j++)
	{
		for(i = 0; i < Nx; i++)
		{
			node = i + j * Nx;

			if(Periodic_boundary_conditions_about_OX)
			{
				ij_function_PeriodicBC_i_and_j_to_the_boundaries(i, j);
			}
			else
			{
				ij_function_i_and_j_to_the_boundaries(node);
			}

			if((fluid < vol_inf_fb[ij])
			&& ((fluid == vol_inf_fb[ij_1]) || (fluid == vol_inf_fb[ij1])))
			{
				rule_vector_ij_BoundaryConditions_Velocities_SlipBC_u[counter] = ij;
				rule_vector_j_BoundaryConditions_Velocities_SlipBC_u[counter] = j;
				counter++;
			}

		}
	}


	is_rule_vectors_BoundaryConditions_Velocities_SlipBC_u_defined = true;


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

	//Define rule vector for velocity slip boundary conditions for v
	if(is_rule_vectors_BoundaryConditions_Velocities_SlipBC_v_defined)
	{
		delete [] rule_vector_ij_BoundaryConditions_Velocities_SlipBC_v;
		delete [] rule_vector_j_BoundaryConditions_Velocities_SlipBC_v;
	}

	//Initiate the number of elements of rule_vector_ij_BoundaryConditions_Velocities_SlipBC_v
	Na_BoundaryConditions_Velocities_SlipBC_v = 0;
	for(j = 0; j < Ny; j++)
	{
		for(i = 0; i < Nx; i++)
		{
			node = i + j * Nx;

			if(Periodic_boundary_conditions_about_OX)
			{
				ij_function_PeriodicBC_i_and_j_to_the_boundaries(i, j);
			}
			else
			{
				ij_function_i_and_j_to_the_boundaries(node);
			}

			//Check is body in this node vol_inf_fb[node] and is fluid up or down of this volume.
			//Here are two velocities v[ij] and v[ij1].
			if((fluid < vol_inf_fb[ij])
			&& ((fluid == vol_inf_fb[i_1j]) || (fluid == vol_inf_fb[i1j])))
			{
				Na_BoundaryConditions_Velocities_SlipBC_v++;
			}

		}
	}

	//Make rule vectors
	rule_vector_ij_BoundaryConditions_Velocities_SlipBC_v = new int [Na_BoundaryConditions_Velocities_SlipBC_v];
	rule_vector_j_BoundaryConditions_Velocities_SlipBC_v = new int [Na_BoundaryConditions_Velocities_SlipBC_v];

	//define elements in rule vector
	counter = 0;
	for(j = 0; j < Ny; j++)
	{
		for(i = 0; i < Nx; i++)
		{
			node = i + j * Nx;

			if(Periodic_boundary_conditions_about_OX)
			{
				ij_function_PeriodicBC_i_and_j_to_the_boundaries(i, j);
			}
			else
			{
				ij_function_i_and_j_to_the_boundaries(node);
			}

			//Check is body in this node vol_inf_fb[node] and is fluid up or down of this volume.
			//Here are two velocities v[ij] and v[ij1].
			if((fluid < vol_inf_fb[ij])
			&& ((fluid == vol_inf_fb[i_1j]) || (fluid == vol_inf_fb[i1j])))
			{
				rule_vector_ij_BoundaryConditions_Velocities_SlipBC_v[counter] = ij;
				rule_vector_j_BoundaryConditions_Velocities_SlipBC_v[counter] = j;
				counter++;
			}

		}
	}

	is_rule_vectors_BoundaryConditions_Velocities_SlipBC_v_defined = true;


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

	//Define rule vector for temperature jump boundary condition on surface of the body
	//Define rule vector for temperature jump boundary conditions
	if(is_rule_vectors_BoundaryConditions_Temperature_TemperatureJumpBC_defined)
	{
		delete [] rule_vector_ij_BoundaryConditions_Temperature_TemperatureJumpBC;
		delete [] rule_vector_j_BoundaryConditions_Temperature_TemperatureJumpBC;
	}

	//Initiate the number of elements of rule_vector_ij_BoundaryConditions_Temperature_TemperatureJumpBC
	Na_BoundaryConditions_Temperature_TemperatureJumpBC = 0;
	for(j = 0; j < Ny; j++)
	{
		for(i = 0; i < Nx; i++)
		{
			node = i + j * Nx;

			if(Periodic_boundary_conditions_about_OX)
			{
				ij_function_PeriodicBC_i_and_j_to_the_boundaries(i, j);
			}
			else
			{
				ij_function_i_and_j_to_the_boundaries(node);
			}

			//Check is control volume ij fluid and is in neighbors control volumes i_1j, i1j, ij_1 and ij1 body.
			if(vol_inf_fb_bool[ij]
			&& ((!vol_inf_fb_bool[i_1j]) || (!vol_inf_fb_bool[i1j])
			    || (!vol_inf_fb_bool[ij_1]) || (!vol_inf_fb_bool[ij1])))
			{
				Na_BoundaryConditions_Temperature_TemperatureJumpBC++;
			}

		}
	}

	//Make rule vectors
	rule_vector_ij_BoundaryConditions_Temperature_TemperatureJumpBC = new int [Na_BoundaryConditions_Temperature_TemperatureJumpBC];
	rule_vector_j_BoundaryConditions_Temperature_TemperatureJumpBC = new int [Na_BoundaryConditions_Temperature_TemperatureJumpBC];

	//define elements in rule vector
	counter = 0;
	for(j = 0; j < Ny; j++)
	{
		for(i = 0; i < Nx; i++)
		{
			node = i + j * Nx;

			if(Periodic_boundary_conditions_about_OX)
			{
				ij_function_PeriodicBC_i_and_j_to_the_boundaries(i, j);
			}
			else
			{
				ij_function_i_and_j_to_the_boundaries(node);
			}

			//Check is control volume ij fluid and is in neighbors control volumes i_1j, i1j, ij_1 and ij1 body.
			if(vol_inf_fb_bool[ij]
			&& ((!vol_inf_fb_bool[i_1j]) || (!vol_inf_fb_bool[i1j])
			    || (!vol_inf_fb_bool[ij_1]) || (!vol_inf_fb_bool[ij1])))
			{
				rule_vector_ij_BoundaryConditions_Temperature_TemperatureJumpBC[counter] = ij;
				rule_vector_j_BoundaryConditions_Temperature_TemperatureJumpBC[counter] = j;
				counter++;
			}

		}
	}

	is_rule_vectors_BoundaryConditions_Temperature_TemperatureJumpBC_defined = true;


}


inline void Solve_rho_in_middle_points(void)
{
	for(j = 0; j < Ny; j++)
	{
		for(i = 0; i < Nx; i++)
		{
			node = i + j * Nx;

			if(Periodic_boundary_conditions_about_OX)
			{
				ij_function_PeriodicBC_i_and_j_to_the_boundaries(i, j);
			}
			else
			{
				ij_function_i_and_j_to_the_boundaries(node);
			}



			if(vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[ij])
			{
				/*Fluid is in two neighbour control volumes*/
				rho_u[ij] = upwind_rho(rho[i_1j], rho[ij], u[ij]);
			}
			else if(vol_inf_fb_bool[i_1j] != vol_inf_fb_bool[ij])
			{
				/*one control volume contains fluid other is body*/
				rho_u[ij] = ((vol_inf_fb_bool[i_1j]) ? (rho[i_1j]) : (rho[ij]));
			}
			else
			{
				/*both control volume is body*/
				rho_u[ij] = 0.0;
			}


			if(vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij])
			{
				/*Fluid is in two neighbor control volumes*/
				rho_v[ij] = upwind_rho(rho[ij_1], rho[ij], v[ij]);
			}
			else if(vol_inf_fb_bool[ij_1] != vol_inf_fb_bool[ij])
			{
				/*one control volume contains fluid other is body*/
				rho_v[ij] = ((vol_inf_fb_bool[ij_1]) ? (rho[ij_1]) : (rho[ij]));
			}
			else
			{
				/*both control volume is body*/
				rho_v[ij] = 0.0;
			}
		}
	}


}


inline void SolveFluxes(void)
{
	//Fluxes in OX and OY.
	for(j = 0; j < Ny; j++)
	{
		for(i = 0; i < Nx; i++)
		{
			node = i + j * Nx;

			Fx[node] = hy[j] * rho_u[node] * u[node];
			Fy[node] = hx[i] * rho_v[node] * v[node];
		}
	}

}


inline void SolveDiffisionTerms(void)
{
	for(j = 0; j < Ny; j++)
	{
		for(i = 0; i < Nx; i++)
		{
			node = i + j * Nx;

			if(Periodic_boundary_conditions_about_OX)
			{
				ij_function_PeriodicBC_i_and_j_to_the_boundaries(i, j);
			}
			else
			{
				ij_function_i_and_j_to_the_boundaries(node);
			}


			//Solve diffusion terms for velocities
			//Solve diffusion terms for u
			if(vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[ij_1]
			&& vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[ij])
			{
				Gamma_yf[ij] = //1.0;
								0.5 * (Li(sqrt_T[i_1j_1], sqrt_T[i_1j], hy[j - 1 * (0 < j)], hy[j])
										+ Li(sqrt_T[ij_1], sqrt_T[ij], hy[j - 1 * (0 < j)], hy[j]));
			}
			else if(vol_inf_fb_bool[i_1j_1] || vol_inf_fb_bool[ij_1]
			|| vol_inf_fb_bool[i_1j] || vol_inf_fb_bool[ij])
			{
				//If there is body, the result is linear interpolation between sqrt from TWO possible
				//temperatures of fluid on the surface of the wall
				Gamma_yf[ij] = //1.0;
								0.5 * ((vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[i_1j]) * Li(sqrt_T[i_1j_1], sqrt_T[i_1j], hy[j - 1 * (0 < j)], hy[j])
										+ (vol_inf_fb_bool[i_1j_1] != vol_inf_fb_bool[i_1j]) * sqrt(Temper_of_fluid_on_the_wall_v[i_1j])

										+ (vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij]) * Li(sqrt_T[ij_1], sqrt_T[ij], hy[j - 1 * (0 < j)], hy[j])
										+ (vol_inf_fb_bool[ij_1] != vol_inf_fb_bool[ij]) * sqrt(Temper_of_fluid_on_the_wall_v[ij])

										+ (  (!(vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[i_1j]))
										  || (!(vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[ij])))
										  * Li(sqrt(Temper_of_fluid_on_the_wall_u[ij_1]), sqrt(Temper_of_fluid_on_the_wall_u[ij]), hy[j - 1 * (0 < j)], hy[j]));
			}


			D_ux[ij] = cT * sqrt_T[i_1j] * (hy[j] / hx[i_1]);


			if((vol_inf_fb_bool[ij] && vol_inf_fb_bool[i_1j])
			|| (vol_inf_fb_bool[ij_1] && vol_inf_fb_bool[i_1j_1]))
			{
				D_uy[ij] = cT * Gamma_yf[ij] * ((hx[i] + hx[i_1])
						/ (vol_inf_fb_bool[ij] * vol_inf_fb_bool[i_1j] * hy[j]
							+ vol_inf_fb_bool[ij_1] * vol_inf_fb_bool[i_1j_1] * hy[j - 1 * (0 < j)]));
			}
			else
			{
				D_uy[ij] = 0;
			}




			//Solve diffusion terms for v
			if(vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[ij_1]
			&& vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[ij])
			{
				Gamma_xf[ij] = //1.0;
								0.5 * (Li(sqrt_T[i_1j_1], sqrt_T[ij_1], hx[i_1], hx[i])
										+ Li(sqrt_T[i_1j], sqrt_T[ij], hx[i_1], hx[i]));
			}
			else if(vol_inf_fb_bool[i_1j_1] || vol_inf_fb_bool[ij_1]
			|| vol_inf_fb_bool[i_1j] || vol_inf_fb_bool[ij])
			{
				//If there is body, the result is linear interpolation between sqrt from TWO possible
				//temperatures of fluid on the surface of the wall
				Gamma_xf[ij] = //1.0;
								0.5 * ((vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[ij_1]) * Li(sqrt_T[i_1j_1], sqrt_T[ij_1], hx[i_1], hx[i])
										+ (vol_inf_fb_bool[i_1j_1] != vol_inf_fb_bool[ij_1]) * sqrt(Temper_of_fluid_on_the_wall_u[ij_1])

										+ (vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[ij]) * Li(sqrt_T[i_1j], sqrt_T[ij], hx[i_1], hx[i])
										+ (vol_inf_fb_bool[i_1j] != vol_inf_fb_bool[ij]) * sqrt(Temper_of_fluid_on_the_wall_u[ij])

										+ (  (!(vol_inf_fb_bool[i_1j_1] && vol_inf_fb_bool[ij_1])) \
										  || (!(vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[ij])))	\
											* Li(sqrt(Temper_of_fluid_on_the_wall_v[i_1j]), sqrt(Temper_of_fluid_on_the_wall_v[ij]), hx[i_1], hx[i]));
			}

			//Solve diffusion terms for v
			if((vol_inf_fb_bool[ij] && vol_inf_fb_bool[ij_1])
			|| (vol_inf_fb_bool[i_1j] && vol_inf_fb_bool[i_1j_1]))
			{
				D_vx[ij] = cT * Gamma_xf[ij] * ((hy[j] + hy[j - 1 * (0 < j)])
					/ (vol_inf_fb_bool[ij] * vol_inf_fb_bool[ij_1] * hx[i]
						+ vol_inf_fb_bool[i_1j] * vol_inf_fb_bool[i_1j_1] * hx[i_1]));
			}
			else
			{
				D_vx[ij] = 0;
			}


			D_vy[ij] = cT * sqrt_T[ij_1] * (hx[i] / hy[j - 1 * (0 < j)]);


			//Solve diffusion terms for Temper
			//Solve DTx
			if(vol_inf_fb_bool[ij] || vol_inf_fb_bool[i_1j])
			{
				DTx[ij] = cT1 * (hy[j]
					/ (0.5 * (vol_inf_fb_bool[ij] * hx[i] + vol_inf_fb_bool[i_1j] * hx[i_1])))

					* (
						//If the volume is to the wall
						((!vol_inf_fb_bool[ij]) || (!vol_inf_fb_bool[i_1j]))* sqrt(Temper_of_fluid_on_the_wall_u[ij])

						//If volume is to the fluid
						+ vol_inf_fb_bool[ij] * vol_inf_fb_bool[i_1j]
							////use Li() or Hi() function
							//* Li(sqrt_T[i_1j], sqrt_T[ij], hx[i_1], hx[i]) //linear interpolation for sqrt_T on boundaries

							* Hi(sqrt_T[i_1j], sqrt_T[ij], hx[i_1], hx[i]) //average harmonic between two given values and steps in mesh
					);
			}
			else
			{
				DTx[ij] = 0;
			}




			//Solve DTy
			if(vol_inf_fb_bool[ij] || vol_inf_fb_bool[ij_1])
			{
				DTy[ij] = cT1 * (hx[i]
					/ (0.5 * (vol_inf_fb_bool[ij] * hy[j] + vol_inf_fb_bool[ij_1] * hy[j - 1 * (0 < j)])))

					* (
						//If the volume is to the wall
							((!vol_inf_fb_bool[ij]) || (!vol_inf_fb_bool[ij_1])) * sqrt(Temper_of_fluid_on_the_wall_v[ij])

						//If volume is to the fluid
						+ vol_inf_fb_bool[ij] * vol_inf_fb_bool[ij_1]
							////use Li() or Hi() function
							//* Li(sqrt_T[ij_1], sqrt_T[ij], hy[j - 1 * (0 < j)], hy[j]) //linear interpolation for sqrt_T on boundaries

							* Hi(sqrt_T[ij_1], sqrt_T[ij], hy[j - 1 * (0 < j)], hy[j]) //average harmonic between two given values and steps in mesh
					);
			}
			else
			{
				DTy[ij] = 0;
			}


		}

	}


}


inline void SolveContinuityEquation(void)
//This procedure solve Continuity Equation and write data to files.
{
	SolveFluxes();

	for(j = 0; j < (Ny - 1); j++)
	{
		for(i = 0; i < (Nx - 1); i++)
		{
			node = i + j * Nx;

			if(Periodic_boundary_conditions_about_OX)
			{
				ij_function_PeriodicBC(i, j);
			}
			else
			{
				ij_function(node);
			}

			ContinuityEquation_stationary[ij] = Fx[i1j] - Fx[ij] + Fy[ij1] - Fy[ij];
			ContinuityEquation[ij] = (rho[ij] - rho_pr[ij]) * hx[i] * hy[j] / ht
				+ ContinuityEquation_stationary[ij];

		}
	}

	if(ToWriteSolvedDataToBinaryFile)
	{
#ifndef COMPILE_MPI
		string ContinuityEquation_stationary_filename = "ContinuityEquation_stationary.bin";
		string ContinuityEquation_filename = "ContinuityEquation.bin";
#endif //End of COMPILE_MPI

#ifdef COMPILE_MPI
		string ContinuityEquation_stationary_filename = "ContinuityEquation_stationary." + toString(IJ) + ".bin";
		string ContinuityEquation_filename = "ContinuityEquation." + toString(IJ) + ".bin";
#endif //End of COMPILE_MPI

		WriteMassiveToBinaryFile(ContinuityEquation_stationary_filename.c_str(), ContinuityEquation_stationary, Na);
		WriteMassiveToBinaryFile(ContinuityEquation_filename.c_str(), ContinuityEquation, Na);


		if(ToWriteSolvedDataToNewFiles)
		{
#ifndef COMPILE_MPI
			string ContinuityEquation_stationary_filename_new = "ContinuityEquation_stationary." + toString(it) + ".bin";
			string ContinuityEquation_filename_new = "ContinuityEquation." + toString(it) + ".bin";
#endif //End of COMPILE_MPI

#ifdef COMPILE_MPI
			string ContinuityEquation_stationary_filename_new = "ContinuityEquation_stationary." + toString(IJ) + "." + toString(it) + ".bin";
			string ContinuityEquation_filename_new = "ContinuityEquation." + toString(IJ) + "." + toString(it) + ".bin";
#endif //End of COMPILE_MPI

			WriteMassiveToBinaryFile(ContinuityEquation_stationary_filename_new.c_str(), ContinuityEquation_stationary, Na);
			WriteMassiveToBinaryFile(ContinuityEquation_filename_new.c_str(), ContinuityEquation, Na);
		}
	}
	else
	{
#ifndef COMPILE_MPI
		string ContinuityEquation_stationary_filename = "ContinuityEquation_stationary.txt";
		string ContinuityEquation_filename = "ContinuityEquation.txt";
#endif //End of COMPILE_MPI

#ifdef COMPILE_MPI
		string ContinuityEquation_stationary_filename = "ContinuityEquation_stationary." + toString(IJ) + ".txt";
		string ContinuityEquation_filename = "ContinuityEquation." + toString(IJ) + ".txt";
#endif //End of COMPILE_MPI

		WriteMassiveToFile_1D(ContinuityEquation_stationary_filename.c_str(), ContinuityEquation_stationary, Na);
		WriteMassiveToFile_1D(ContinuityEquation_filename.c_str(), ContinuityEquation, Na);


		if(ToWriteSolvedDataToNewFiles)
		{
#ifndef COMPILE_MPI
			string ContinuityEquation_stationary_filename_new = "ContinuityEquation_stationary." + toString(it) + ".txt";
			string ContinuityEquation_filename_new = "ContinuityEquation." + toString(it) + ".txt";
#endif //End of COMPILE_MPI

#ifdef COMPILE_MPI
			string ContinuityEquation_stationary_filename_new = "ContinuityEquation_stationary." + toString(IJ) + "." + toString(it) + ".txt";
			string ContinuityEquation_filename_new = "ContinuityEquation." + toString(IJ) + "." + toString(it) + ".txt";
#endif //End of COMPILE_MPI

			WriteMassiveToFile_1D(ContinuityEquation_stationary_filename_new.c_str(), ContinuityEquation_stationary, Na);
			WriteMassiveToFile_1D(ContinuityEquation_filename_new.c_str(), ContinuityEquation, Na);
		}

	}

}




//Here are placed fincions and procedures for MPI
#ifdef COMPILE_MPI

void DefineSubDomains_MPI(void)
{
	//This is made to read data from EnteredData.txt from the same way
	Nx_entire_computational_domain = Nx;
	Ny_entire_computational_domain = Ny;

	Na_entire_computational_domain = Nx_entire_computational_domain * Ny_entire_computational_domain;

	//Calculate N_SubDomains_a for solving nodes in subdomains
	N_SubDomains_a = N_SubDomains_x * N_SubDomains_y;

	SubDomain = new DomainDecomposition2D[N_SubDomains_a];

	//Target of this separation of subdomains is to make subdomains with
	//almost equal nodes.

	//Calculate number of nodes in all subdomains exept last subdomain on OX
	double Nx_all_N_SubDomains_x_tmp;

	//Calculate number of nodes in all subdomains exept last subdomain on OY
	double Ny_all_N_SubDomains_y_tmp;

	//Separate SubDomains in OX and OY direction:
	//Mantissa of the nodes Nx and Ny
	double mantissa_x, mantissa_y;
	mantissa_x = 0;
	mantissa_y = 0;

	//Sum of number of nodes, which are calculated in all subdomains
	//exept last one. In this number are not included points
	//in boundary of neighbourhood domain, if this one exist.
	unsigned int Nx_all_exept_last, Ny_all_exept_last;
	Nx_all_exept_last = 0;
	Ny_all_exept_last = 0;


	//If the calculation continue Nx and Ny for all subdomains have to be readed
	//from file Nx_Ny_ForAllSubDomains.txt
	if(ToReadSolvedDataFromFile)
	{
		ReadDataFor_Nx_Ny_ForAllSubDomains();
	}
	else
	{
		if(is_defined_arrays_for_Nx_Ny_for_all_SubDomains)
		{
			Nx_all_SubDomains = NULL; delete [] Nx_all_SubDomains;
			Ny_all_SubDomains = NULL; delete [] Ny_all_SubDomains;
		}

		Nx_all_SubDomains = new int [N_SubDomains_a];
		Ny_all_SubDomains = new int [N_SubDomains_a];

		is_defined_arrays_for_Nx_Ny_for_all_SubDomains = true;
	}


	//Define nodes on every subdomains on OX and OY
	for(J = 0; J < N_SubDomains_y; J++)
	{
		for(I = 0; I < N_SubDomains_x; I++)
		{
			//Calculate indexe
			IJ = I + J * N_SubDomains_x;

			//Define variable in SubDomain[IJ]
			SubDomain[IJ].N_SubDomains_a = N_SubDomains_a;
			SubDomain[IJ].N_SubDomains_x = N_SubDomains_x;
			SubDomain[IJ].N_SubDomains_y = N_SubDomains_y;
			SubDomain[IJ].N_neighbouhoods = N_neighbouhoods;
			SubDomain[IJ].N_exchanged_variables_with_neighbouhoods = N_exchanged_variables_with_neighbouhoods;
			SubDomain[IJ].N_send_variables_to_IJ0 = N_send_variables_to_IJ0;
			SubDomain[IJ].N_send_from_IJ0 = N_send_from_IJ0;

			SubDomain[IJ].DefineIndexesAndTags(I, J, Periodic_boundary_conditions_about_OX);

			I_1J = SubDomain[IJ].I_1J;
			I1J = SubDomain[IJ].I1J;
			IJ_1 = SubDomain[IJ].IJ_1;
			IJ1 = SubDomain[IJ].IJ1;
			//End of Define variable in SubDomain[IJ]

			//Define nodes in SubDomain[IJ] on OX
			if(J == 0)
			{
				//The nodes will be calculated when J == 0.
				//For 0 < J the result will be copied.
				if(I == 0)
				{
					//The first computational domain IJ == 0 collect information for residuals for all computational domains,
					//therefore the number cells for calculations have to smaller
					//Here the multiplication by 1.0 is to make divide type double of int variables

					//The reduction of Nx of subdomain is linear interpolation.
					//For N_processes_min the part of all is  Part_of_all_min
					int N_processes_min = 2;
					double Part_of_all_min = 0.95;

					//For N_processes_max the part of all is  Part_of_all_max
					int N_processes_max = 256;
					double Part_of_all_max = 0.5;
					
					//Munimum value of Nx on subdomain
					int Nx_subdimain_min = 10;

					if(N_SubDomains_a < N_processes_max)
					{
						//Number of processes N_SubDomains_a belongs to [N_processes_min, N_processes_max]
						//Reduce coefficient - reduce the number of cells on OX (Nx) of subdomain in function of number of processes (N_SubDomains_a) - linear interpolation
						double reduce_coefficient;
						reduce_coefficient = Part_of_all_min + (N_SubDomains_a - N_processes_min) * (Part_of_all_max - Part_of_all_min) / (N_processes_max - N_processes_min);
						//The coefficient is limited beteen [Part_of_all_min; Part_of_all_max]
						if(reduce_coefficient < Part_of_all_min) reduce_coefficient = maximum(Part_of_all_max, reduce_coefficient);
						if(Part_of_all_max < reduce_coefficient) reduce_coefficient = minimum(Part_of_all_min, reduce_coefficient);


						Nx_all_N_SubDomains_x_tmp = (int)((reduce_coefficient * Nx_entire_computational_domain / N_SubDomains_x));
						
						if(Nx_all_N_SubDomains_x_tmp < Nx_subdimain_min)
						{
							//Nx_all_N_SubDomains_x_tmp have to be larger then Nx_subdimain_min
							Nx_all_N_SubDomains_x_tmp = Nx_subdimain_min;
						}
					}
					else
					{
						//Number of processes is larger then N_processes_max < N_SubDomains_a
						//if N_processes_max == 100, then 1 process is less then a percent,
						//therefore for the I == 0, Nx = Nx_subdimain_min
						//The target of this to releas process 0 from calculations to be able to start non-blocking communications
						//as faster as possible.
						//The runs show that the speedup decrease, when are used more then 200-256 processes
						Nx_all_N_SubDomains_x_tmp = Nx_subdimain_min;
					}
					
					

					if(ToReadSolvedDataFromFile && (!ToOnlyWriteDataFor_Nx_Ny_ForAllSubDomains_AndExit))
					{
						SubDomain[IJ].Nx = Nx_all_SubDomains[IJ];
					}
					else
					{
						//In nodes have to be included and nodes for calculated in
						//neighbourhood subdomains.
						SubDomain[IJ].Nx = int(Nx_all_N_SubDomains_x_tmp)
							+ BlockLengthSwap_I_1 * SubDomain[IJ].is_SubDomain_in_direction_I_1
							+ BlockLengthSwap_I1 * SubDomain[IJ].is_SubDomain_in_direction_I1;

						Nx_all_SubDomains[IJ] = SubDomain[IJ].Nx;
					}


					Nx_all_exept_last += SubDomain[IJ].Nx
						- (BlockLengthSwap_I_1 * SubDomain[IJ].is_SubDomain_in_direction_I_1
							+ BlockLengthSwap_I1 * SubDomain[IJ].is_SubDomain_in_direction_I1);

				}
				else if(I < (N_SubDomains_x - 1))
				{
					//Here the multiplication by 1.0 is to make divide type double of int variables
					Nx_all_N_SubDomains_x_tmp = 1.0 * (Nx_entire_computational_domain - SubDomain[0].Nx) / (N_SubDomains_x - 1);

					mantissa_x += Nx_all_N_SubDomains_x_tmp - floor(Nx_all_N_SubDomains_x_tmp);

					if(ToReadSolvedDataFromFile && (!ToOnlyWriteDataFor_Nx_Ny_ForAllSubDomains_AndExit))
					{
						SubDomain[IJ].Nx = Nx_all_SubDomains[IJ];
					}
					else
					{
						//In nodes have to be included and nodes for calculated in
						//neighbourhood subdomains.
						SubDomain[IJ].Nx = int(Nx_all_N_SubDomains_x_tmp) + int(mantissa_x)
							+ BlockLengthSwap_I_1 * SubDomain[IJ].is_SubDomain_in_direction_I_1
							+ BlockLengthSwap_I1 * SubDomain[IJ].is_SubDomain_in_direction_I1;

						Nx_all_SubDomains[IJ] = SubDomain[IJ].Nx;
					}


					if(1 <= int(mantissa_x))
					{
						mantissa_x--;
					}

					Nx_all_exept_last += SubDomain[IJ].Nx
						- (BlockLengthSwap_I_1 * SubDomain[IJ].is_SubDomain_in_direction_I_1
							+ BlockLengthSwap_I1 * SubDomain[IJ].is_SubDomain_in_direction_I1);
				}
				else
				{
					if(ToReadSolvedDataFromFile && (!ToOnlyWriteDataFor_Nx_Ny_ForAllSubDomains_AndExit))
					{
						SubDomain[IJ].Nx = Nx_all_SubDomains[IJ];
					}
					else
					{
						//Calculate number of nodes in last subdomain on OX
						SubDomain[IJ].Nx = Nx_entire_computational_domain - Nx_all_exept_last
							+ BlockLengthSwap_I_1 * SubDomain[IJ].is_SubDomain_in_direction_I_1
							+ BlockLengthSwap_I1 * SubDomain[IJ].is_SubDomain_in_direction_I1;

						Nx_all_SubDomains[IJ] = SubDomain[IJ].Nx;
					}

				}


				//Define mesh for this subdomain will be made later,
				//here will be defined begin and end index on OX.
				if(I == 0)
				{
					SubDomain[IJ].i_b = 0;
					SubDomain[IJ].i_e = SubDomain[IJ].Nx - 1
										- BlockLengthSwap_I1 * SubDomain[IJ].is_SubDomain_in_direction_I1;
				}
				else
				{
					////Calculate number of first node of this subdomain, when counting
					////number of nodes in all subdomains on the same exept this
					//int Nx_all_to_this_subdomain;
					//Nx_all_to_this_subdomain = 0;
					//for(counter = 0; counter < I; counter++)
					//{
					//	Nx_all_to_this_subdomain +=	SubDomain[counter].Nx
					//									- (BlockLengthSwap_I_1 * SubDomain[counter].is_SubDomain_in_direction_I_1
					//										+ BlockLengthSwap_I1 * SubDomain[counter].is_SubDomain_in_direction_I1);
					//}


					SubDomain[IJ].i_b = SubDomain[IJ - 1].i_e + 1;

					SubDomain[IJ].i_e = SubDomain[IJ].i_b
										+ SubDomain[IJ].Nx - 1
											- (BlockLengthSwap_I_1 * SubDomain[IJ].is_SubDomain_in_direction_I_1
												+ BlockLengthSwap_I1 * SubDomain[IJ].is_SubDomain_in_direction_I1);

				}

			}
			else
			{
				//The nodes are already calculated, when J == 0.
				//For 0 < J the result will be copied.
				SubDomain[IJ].Nx = SubDomain[IJ_1].Nx;
				SubDomain[IJ].i_b = SubDomain[IJ_1].i_b;
				SubDomain[IJ].i_e = SubDomain[IJ_1].i_e;

				Nx_all_SubDomains[IJ] = SubDomain[IJ].Nx;
			}
			//End of Define nodes in SubDomain[IJ] on OX


			//Define nodes in SubDomain[IJ] on OY
			if(I == 0)
			{
				//The nodes will be calculated when I == 0.
				//For 0 < I the result will be copied.
				if(J < (N_SubDomains_y - 1))
				{
					//Here the multiplication by 1.0 is to make divide type double of int variables
					Ny_all_N_SubDomains_y_tmp = 1.0 * Ny_entire_computational_domain / N_SubDomains_y;

					mantissa_y += Ny_all_N_SubDomains_y_tmp - floor(Ny_all_N_SubDomains_y_tmp);

					if(ToReadSolvedDataFromFile && (!ToOnlyWriteDataFor_Nx_Ny_ForAllSubDomains_AndExit))
					{
						SubDomain[IJ].Ny = Ny_all_SubDomains[IJ];
					}
					else
					{
						//In nodes have to be included and nodes for calculated in
						//neighbourhood subdomains.
						SubDomain[IJ].Ny = int(Ny_all_N_SubDomains_y_tmp) + int(mantissa_y)
							+ BlockLengthSwap_J_1 * SubDomain[IJ].is_SubDomain_in_direction_J_1
							+ BlockLengthSwap_J1 * SubDomain[IJ].is_SubDomain_in_direction_J1;

						Ny_all_SubDomains[IJ] = SubDomain[IJ].Ny;
					}


					if(1 <= int(mantissa_y))
					{
						mantissa_y--;
					}

					Ny_all_exept_last += SubDomain[IJ].Ny
						- (BlockLengthSwap_J_1 * SubDomain[IJ].is_SubDomain_in_direction_J_1
							+ BlockLengthSwap_J1 * SubDomain[IJ].is_SubDomain_in_direction_J1);
				}
				else
				{
					if(ToReadSolvedDataFromFile && (!ToOnlyWriteDataFor_Nx_Ny_ForAllSubDomains_AndExit))
					{
						SubDomain[IJ].Ny = Ny_all_SubDomains[IJ];
					}
					else
					{
						//Calculate number of nodes in last subdomain on OY
						SubDomain[IJ].Ny = Ny_entire_computational_domain - Ny_all_exept_last
							+ BlockLengthSwap_J_1 * SubDomain[IJ].is_SubDomain_in_direction_J_1
							+ BlockLengthSwap_J1 * SubDomain[IJ].is_SubDomain_in_direction_J1;

						Ny_all_SubDomains[IJ] = SubDomain[IJ].Ny;
					}

				}

				//Define mesh for this subdomain will be made later,
				//here will be defined begin and end index on OY.
				if(J == 0)
				{
					SubDomain[IJ].j_b = 0;
					SubDomain[IJ].j_e = SubDomain[IJ].Ny - 1
										- BlockLengthSwap_J1 * SubDomain[IJ].is_SubDomain_in_direction_J1;
				}
				else
				{
					////Calculate number of first node of this subdomain, when counting
					////number of nodes in all subdomains on the same exept this
					//int Ny_all_to_this_subdomain;
					//Ny_all_to_this_subdomain = 0;
					//for(counter = 0; counter < J; counter++)
					//{
					//	Ny_all_to_this_subdomain +=	SubDomain[counter].Ny
					//									- (BlockLengthSwap_J_1 * SubDomain[counter].is_SubDomain_in_direction_J_1
					//										+ BlockLengthSwap_J1 * SubDomain[counter].is_SubDomain_in_direction_J1);
					//}


					SubDomain[IJ].j_b = SubDomain[IJ - N_SubDomains_x].j_e + 1;

					SubDomain[IJ].j_e = SubDomain[IJ].j_b
										+ SubDomain[IJ].Ny - 1
											- (BlockLengthSwap_J_1 * SubDomain[IJ].is_SubDomain_in_direction_J_1
												+ BlockLengthSwap_J1 * SubDomain[IJ].is_SubDomain_in_direction_J1);

				}
			}
			else
			{
				//The nodes are already calculated, when I == 0.
				//For 0 < I the result will be copied.
				SubDomain[IJ].Ny = SubDomain[I_1J].Ny;
				SubDomain[IJ].j_b = SubDomain[I_1J].j_b;
				SubDomain[IJ].j_e = SubDomain[I_1J].j_e;

				Ny_all_SubDomains[IJ] = SubDomain[IJ].Ny;
			}
			//End of Define nodes in SubDomain[IJ] on OY

			//Calculate all nodes in subdomain IJ
			SubDomain[IJ].Na = SubDomain[IJ].Nx * SubDomain[IJ].Ny;


		}


	}

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Comm_rank(MPI_COMM_WORLD, &IJ);


	//This is option to write only number of subdomains if there is need to make change.
	//The program have to be made to start from begining just for this this change.
	if((IJ == 0) && ToOnlyWriteDataFor_Nx_Ny_ForAllSubDomains_AndExit)
	{
		WriteDataFor_Nx_Ny_ForAllSubDomains();
		std::cout << "Program write data for Nx Ny for all SubDomains and now will exit."
			<< "\n";
	}

	MPI_Comm_size(MPI_COMM_WORLD, &N_SubDomains_a);


	//Define MPI variables
	N_max_requests =
		N_exchanged_variables_with_neighbouhoods * N_actions * N_neighbouhoods
		+ N_send_variables_to_IJ0
		+ N_send_from_IJ0;

	request = new MPI_Request[N_max_requests];
	status = new MPI_Status[N_max_requests];



	//Get version of MPI
	if(IJ == 0)
	{
		int version, subversion;
		MPI_Get_version(&version, &subversion);

		std::cout << "MPI version.subversion: " << version << "." << subversion << "\n";
	}

	//Indexes if this SubDomain are:
	I = SubDomain[IJ].I;
	J = SubDomain[IJ].J;

	//Neighbourhoods are
	I_1J = SubDomain[IJ].I_1J;
	I1J = SubDomain[IJ].I1J;
	IJ_1 = SubDomain[IJ].IJ_1;
	IJ1 = SubDomain[IJ].IJ1;

	//Make arrayes for IJ == 0
	if(IJ == 0)
	{
		max_residual_in_u_all_subdomains = new double [N_SubDomains_a];
		max_residual_in_v_all_subdomains = new double [N_SubDomains_a];
		max_residual_in_p_all_subdomains = new double [N_SubDomains_a];
		max_residual_in_T_all_subdomains = new double [N_SubDomains_a];

		request_max_residual_in_u = new MPI_Request[N_SubDomains_a];
		status_max_residual_in_u = new MPI_Status[N_SubDomains_a];

		request_max_residual_in_v = new MPI_Request[N_SubDomains_a];
		status_max_residual_in_v = new MPI_Status[N_SubDomains_a];

		request_max_residual_in_p = new MPI_Request[N_SubDomains_a];
		status_max_residual_in_p = new MPI_Status[N_SubDomains_a];

		request_max_residual_in_T = new MPI_Request[N_SubDomains_a];
		status_max_residual_in_T = new MPI_Status[N_SubDomains_a];
	}


	//Copy mesh fa all computational domain and create mesh for this subdomain
	//Copy mesh of all computational domain
	//Make arrays for mesh for all computational domain
	if(!is_defined_array_for_mesh_entire_computational_domain)
	{
		x_f_entire_computational_domain = new double [Nx + 1];
		y_f_entire_computational_domain = new double [Ny + 1];

		hx_entire_computational_domain = new double [Nx];
		hy_entire_computational_domain = new double [Ny];

		x_v_entire_computational_domain = new double [Nx];
		y_v_entire_computational_domain = new double [Ny];

		is_defined_array_for_mesh_entire_computational_domain = true;
	}

	//Copy data from x_f, y_f, hx, hy, x_v and y_v
	//to x_f_entire_computational_domain, y_f_entire_computational_domain, hx_entire_computational_domain, hy_entire_computational_domain, x_v_entire_computational_domain and y_v_entire_computational_domain
	for(i = 0; i < Nx_entire_computational_domain; i++)
	{
		x_f_entire_computational_domain[i] = x_f[i];
		hx_entire_computational_domain[i] = hx[i];
		x_v_entire_computational_domain[i] = x_v[i];
	}

	//copy data for the last node on OX
	x_f_entire_computational_domain[Nx_entire_computational_domain] = x_f[Nx_entire_computational_domain];

	for(j = 0; j < Ny_entire_computational_domain; j++)
	{
		y_f_entire_computational_domain[j] = y_f[j];
		hy_entire_computational_domain[j] = hy[j];
		y_v_entire_computational_domain[j] = y_v[j];
	}

	//Copy data for the last node on OY
	y_f_entire_computational_domain[Ny_entire_computational_domain] = y_f[Ny_entire_computational_domain];


	//Number of nodes in this calculation for this subdomain
	Nx = SubDomain[IJ].Nx;
	Ny = SubDomain[IJ].Ny;
	Na = SubDomain[IJ].Na;

	//Make mesh arraies for this subdomain
	if(is_defined_array_for_mesh)
	{
		x_f = NULL; delete [] x_f;
		y_f = NULL; delete [] y_f;

		hx = NULL; delete [] hx;
		hy = NULL; delete [] hy;

		x_v = NULL; delete [] x_v;
		y_v = NULL; delete [] y_v;

		is_defined_array_for_mesh = false;
	}


	if(!is_defined_array_for_mesh)
	{
		x_f = new double [Nx + 1];
		y_f = new double [Ny + 1];

		hx = new double [Nx];
		hy = new double [Ny];

		x_v = new double [Nx];
		y_v = new double [Ny];

		is_defined_array_for_mesh = true;
	}



	//Copy data for mesh of this subdomain
	//from x_f_entire_computational_domain, y_f_entire_computational_domain, hx_entire_computational_domain, hy_entire_computational_domain, x_v_entire_computational_domain and y_v_entire_computational_domain
	//Copy data from x_f_entire_computational_domain, y_f_entire_computational_domain, hx_entire_computational_domain, hy_entire_computational_domain, x_v_entire_computational_domain and y_v_entire_computational_domain
	//to x_f, y_f, hx, hy, x_v and y_v
	if(Periodic_boundary_conditions_about_OX)
	{
		for(i = 1; i < (SubDomain[IJ].Nx - 1); i++)
		{
			x_f[i] = x_f_entire_computational_domain[SubDomain[IJ].i_b + i];
			hx[i] = hx_entire_computational_domain[SubDomain[IJ].i_b + i];
			x_v[i] = x_v_entire_computational_domain[SubDomain[IJ].i_b + i];
		}

		if(N_SubDomains_x == 1)
		{
			//If on OX there only one subdomain
			x_f[0] = x_f_entire_computational_domain[0];
			hx[0] = hx_entire_computational_domain[0];
			x_v[0] = x_v_entire_computational_domain[0];


			x_f[Nx - 1] = x_f_entire_computational_domain[Nx - 1];
			hx[Nx - 1] = hx_entire_computational_domain[Nx - 1];
			x_v[Nx - 1] = x_v_entire_computational_domain[Nx - 1];

			x_f[Nx] = x_f_entire_computational_domain[Nx];
		}
		else
		{
			//If on OX there are more then one subdomains on OX
			if(I == 0)
			{
				//This is zero subdomain on OX
				hx[0] = hx_entire_computational_domain[Nx_entire_computational_domain - 1];
				x_f[0] = x_f_entire_computational_domain[1] - hx[0];
				x_v[0] = x_v_entire_computational_domain[0] - 0.5 * (hx[1] + hx[0]);


				hx[Nx - 1] = hx_entire_computational_domain[SubDomain[IJ].i_b + Nx - 1];
				x_f[Nx - 1] = x_f_entire_computational_domain[SubDomain[IJ].i_b + Nx - 1];
				x_v[Nx - 1] = x_v_entire_computational_domain[SubDomain[IJ].i_b + Nx - 1];

				x_f[Nx] = x_f_entire_computational_domain[SubDomain[IJ].i_b + Nx];
			}
			else if(I == (N_SubDomains_x - 1))
			{
				//This is last subdomain on OX
				hx[0] = hx_entire_computational_domain[SubDomain[IJ].i_b];
				x_f[0] = x_f_entire_computational_domain[SubDomain[IJ].i_b];
				x_v[0] = x_v_entire_computational_domain[SubDomain[IJ].i_b];


				hx[Nx - 1] = hx_entire_computational_domain[0];
				x_f[Nx - 1] = x_f[Nx - 2] + hx[Nx - 2];
				x_v[Nx - 1] = x_v[Nx - 2] + 0.5 * (hx[Nx - 2] + hx[Nx - 1]);

				x_f[Nx] = x_f[Nx - 1] + hx[Nx - 1];
			}
			else
			{
				//This subdomain is between zero and last on OX.
				//(No zero or last on OX.)
				x_f[0] = x_f_entire_computational_domain[SubDomain[IJ].i_b];
				hx[0] = hx_entire_computational_domain[SubDomain[IJ].i_b];
				x_v[0] = x_v_entire_computational_domain[SubDomain[IJ].i_b];


				x_f[Nx - 1] = x_f_entire_computational_domain[SubDomain[IJ].i_b + Nx - 1];
				hx[Nx - 1] = hx_entire_computational_domain[SubDomain[IJ].i_b + Nx - 1];
				x_v[Nx - 1] = x_v_entire_computational_domain[SubDomain[IJ].i_b + Nx - 1];

				x_f[Nx] = x_f_entire_computational_domain[SubDomain[IJ].i_b + Nx];
			}

		}

	}
	else
	{
		for(i = 0; i < SubDomain[IJ].Nx; i++)
		{
			x_f[i] = x_f_entire_computational_domain[SubDomain[IJ].i_b - BlockLengthSwap_I_1 * SubDomain[IJ].is_SubDomain_in_direction_I_1+ i];
			hx[i] = hx_entire_computational_domain[SubDomain[IJ].i_b - BlockLengthSwap_I_1 * SubDomain[IJ].is_SubDomain_in_direction_I_1 + i];
			x_v[i] = x_v_entire_computational_domain[SubDomain[IJ].i_b - BlockLengthSwap_I_1 * SubDomain[IJ].is_SubDomain_in_direction_I_1 + i];
		}

		//copy data for the last node on OX
		x_f[Nx] = x_f_entire_computational_domain[SubDomain[IJ].i_b - BlockLengthSwap_I_1 * SubDomain[IJ].is_SubDomain_in_direction_I_1 + Nx];
	}


	for(j = 0; j < Ny; j++)
	{
		y_f[j] = y_f_entire_computational_domain[SubDomain[IJ].j_b - BlockLengthSwap_J_1 * SubDomain[IJ].is_SubDomain_in_direction_J_1 + j];
		hy[j] = hy_entire_computational_domain[SubDomain[IJ].j_b - BlockLengthSwap_J_1 * SubDomain[IJ].is_SubDomain_in_direction_J_1 + j];
		y_v[j] = y_v_entire_computational_domain[SubDomain[IJ].j_b - BlockLengthSwap_J_1 * SubDomain[IJ].is_SubDomain_in_direction_J_1 + j];
	}

	//copy data for the last node on OY
	y_f[Ny] = y_f_entire_computational_domain[SubDomain[IJ].j_b - BlockLengthSwap_J_1 * SubDomain[IJ].is_SubDomain_in_direction_J_1 + Ny];



	//Save mesh for this subdomain
	string x_f_filename = "x_f." + toString(IJ) + ".txt";
	string hx_filename = "hx." + toString(IJ) + ".txt";
	string x_v_filename = "x_v." + toString(IJ) + ".txt";

	string y_f_filename = "y_f." + toString(IJ) + ".txt";
	string hy_filename = "hy." + toString(IJ) + ".txt";
	string y_v_filename = "y_v." + toString(IJ) + ".txt";

	//write mesh to files
	WriteMassiveToFile_1D(x_f_filename.c_str(), x_f, Nx+1);
	WriteMassiveToFile_1D(hx_filename.c_str(), hx, Nx);
	WriteMassiveToFile_1D(x_v_filename.c_str(), x_v, Nx);

	WriteMassiveToFile_1D(y_f_filename.c_str(), y_f, Ny+1);
	WriteMassiveToFile_1D(hy_filename.c_str(), hy, Ny);
	WriteMassiveToFile_1D(y_v_filename.c_str(), y_v, Ny);


	//Save mesh for entire computational domain
	if(IJ == 0)
	{
		//write mesh to files
		WriteMassiveToFile_1D("x_f.txt", x_f_entire_computational_domain, Nx_entire_computational_domain+1);
		WriteMassiveToFile_1D("hx.txt", hx_entire_computational_domain, Nx_entire_computational_domain);
		WriteMassiveToFile_1D("x_v.txt", x_v_entire_computational_domain, Nx_entire_computational_domain);

		WriteMassiveToFile_1D("y_f.txt", y_f_entire_computational_domain, Ny_entire_computational_domain+1);
		WriteMassiveToFile_1D("hy.txt", hy_entire_computational_domain, Ny_entire_computational_domain);
		WriteMassiveToFile_1D("y_v.txt", y_v_entire_computational_domain, Ny_entire_computational_domain);
	}


	//Make MPI_Datatypes to exchange data with neighbourhoods
	//Exchange data with subdomain I_1J, if exist.
	if(SubDomain[IJ].is_SubDomain_in_direction_I_1)
	{
		if(!is_created_MPI_Datatype_send_to_I_1J)
		{
			//Number of elements, blocklength and stride which will be send to subdomain I_1J
			N_send_to_I_1J = Ny;

			//Send BlockLengthSwap_I1, because this is swap I1J of subdomain I_1J.
			blocklength_send_to_I_1J = BlockLengthSwap_I1;

			stride_send_to_I_1J = Nx;


			MPI_Type_vector(N_send_to_I_1J, blocklength_send_to_I_1J, stride_send_to_I_1J,
				MPI_DOUBLE, &send_to_I_1J);

			MPI_Type_commit(&send_to_I_1J);
			//End of Create derived datatype MPI_TYPE_VECTOR to send data to I_1J

			is_created_MPI_Datatype_send_to_I_1J = true;
		}


		if(!is_created_MPI_Datatype_recv_from_I_1J)
		{
			//Number of elements, blocklength and stride which will be receive from subdomain I_1J
			N_recv_from_I_1J = Ny;
			blocklength_recv_from_I_1J = BlockLengthSwap_I_1;
			stride_recv_from_I_1J = Nx;

			MPI_Type_vector(N_recv_from_I_1J, blocklength_recv_from_I_1J, stride_recv_from_I_1J,
				MPI_DOUBLE, &recv_from_I_1J);

			MPI_Type_commit(&recv_from_I_1J);
			//End of Create derived datatype MPI_TYPE_VECTOR to receive data from I_1J

			is_created_MPI_Datatype_recv_from_I_1J = true;
		}

	}


	//Exchange data with subdomain IJ1, if exist.
	if(SubDomain[IJ].is_SubDomain_in_direction_I1)
	{
		if(!is_created_MPI_Datatype_send_to_I1J)
		{
			//Number of elements, blocklength and stride which will be send to subdomain I1J
			N_send_to_I1J = Ny;

			//Send BlockLengthSwap_I_1, because this is swap I_1J of subdomain I1J.
			blocklength_send_to_I1J = BlockLengthSwap_I_1;

			stride_send_to_I1J = Nx;


			MPI_Type_vector(N_send_to_I1J, blocklength_send_to_I1J, stride_send_to_I1J,
				MPI_DOUBLE, &send_to_I1J);

			MPI_Type_commit(&send_to_I1J);
			//End of Create derived datatype MPI_TYPE_VECTOR to send data to I1J

			is_created_MPI_Datatype_send_to_I1J = true;
		}

		if(!is_created_MPI_Datatype_recv_from_I1J)
		{
			//Number of elements, blocklength and stride which will be receive from subdomain I1J
			N_recv_from_I1J = Ny;
			blocklength_recv_from_I1J = BlockLengthSwap_I1;
			stride_recv_from_I1J = Nx;

			MPI_Type_vector(N_recv_from_I1J, blocklength_recv_from_I1J, stride_recv_from_I1J,
				MPI_DOUBLE, &recv_from_I1J);

			MPI_Type_commit(&recv_from_I1J);
			//End of Create derived datatype MPI_TYPE_VECTOR to receive data from I1J

			is_created_MPI_Datatype_recv_from_I1J = true;
		}

	}



}


template <typename T_massive, typename T_Nx_limits_massive, typename T_Ny_limits_massive>
inline void StartDataExchangeWithNeighbourhoodsSubDomains_MPI(
	T_massive massive[],
	T_Nx_limits_massive Nx_massive_limits[],
	T_Ny_limits_massive Ny_massive_limits[])
{
	//Define tags and index of requests for data exchange.
	//Define tags for data exchange.
	//tag_send_fi_to_I1J_tmp - this mean that is send to I1J from IJ
	//tag_send_fi_to_I1J_tmp_from_I_1J - this mean the is send to I1J from SubDomain I_1J
	int tag_send_fi_to_I_1J_tmp, tag_send_fi_to_I1J_tmp_from_I_1J,
		tag_send_fi_to_I1J_tmp, tag_send_fi_to_I_1J_tmp_from_I1J,
		tag_send_fi_to_IJ_1_tmp, tag_send_fi_to_IJ1_tmp_from_IJ_1,
		tag_send_fi_to_IJ1_tmp, tag_send_fi_to_IJ_1_tmp_from_IJ1;

	//Define index of requests for data exchange.
	int index_request_send_fi_to_I_1J_tmp, index_request_recv_fi_from_I_1J_tmp,
		index_request_send_fi_to_I1J_tmp, index_request_recv_fi_from_I1J_tmp,
		index_request_send_fi_to_IJ_1_tmp, index_request_recv_fi_from_IJ_1_tmp,
		index_request_send_fi_to_IJ1_tmp, index_request_recv_fi_from_IJ1_tmp;

	int variable_send_tmp;


	//Ininicialize tags
	tag_send_fi_to_I_1J_tmp = 0;
	tag_send_fi_to_I1J_tmp_from_I_1J = 0;
	tag_send_fi_to_I1J_tmp = 0;
	tag_send_fi_to_I_1J_tmp_from_I1J = 0;
	tag_send_fi_to_IJ_1_tmp = 0;
	tag_send_fi_to_IJ1_tmp_from_IJ_1 = 0;
	tag_send_fi_to_IJ1_tmp = 0;
	tag_send_fi_to_IJ_1_tmp_from_IJ1 = 0;



	//Ininicialize index of requests
	index_request_send_fi_to_I_1J_tmp = 0;
	index_request_recv_fi_from_I_1J_tmp = 0;
	index_request_send_fi_to_I1J_tmp = 0;
	index_request_recv_fi_from_I1J_tmp = 0;
	index_request_send_fi_to_IJ_1_tmp = 0;
	index_request_recv_fi_from_IJ_1_tmp = 0;
	index_request_send_fi_to_IJ1_tmp = 0;
	index_request_recv_fi_from_IJ1_tmp = 0;


	//Calculate tags and index of requests in depend of exchanged variable.
	//The variable will be recognized in depend of address of massive.
	if(massive == d_p_c_12)
	{
		//Exchange array is d_p_c_12
		variable_send_tmp = d_p_c_12_send;
	}
	else if(massive == apx)
	{
		//Exchange array is apx
		variable_send_tmp = apx_send;
	}
	else if(massive == bpx)
	{
		//Exchange array is bpx
		variable_send_tmp = bpx_send;
	}
	else if(massive == d_p_c_34)
	{
		//Exchange array is d_p_c_34
		variable_send_tmp = d_p_c_34_send;
	}
	else if(massive == apy)
	{
		//Exchange array is apy
		variable_send_tmp = apy_send;
	}
	else if(massive == bpy)
	{
		//Exchange array is bpy
		variable_send_tmp = bpy_send;
	}
	else if(massive == p)
	{
		//Exchange array is p
		variable_send_tmp = p_send;
	}
	else if(massive == Temper)
	{
		//Exchange array is Temper
		variable_send_tmp = Temper_send;
	}
	else if(massive == u)
	{
		//Exchange array is u
		variable_send_tmp = u_send;
	}
	else if(massive == v)
	{
		//Exchange array is v
		variable_send_tmp = v_send;
	}



	//Exchange variable_send_tmp
	//Exchange data with subdomain I_1J, if exist.
	if(SubDomain[IJ].is_SubDomain_in_direction_I_1)
	{
		tag_send_fi_to_I_1J_tmp = SubDomain[IJ].CalculateTagForSendVariableToNeighbour(
			variable_send_tmp, to_I_1J);

		tag_send_fi_to_I1J_tmp_from_I_1J = SubDomain[I_1J].CalculateTagForSendVariableToNeighbour(
			variable_send_tmp, to_I1J);


		index_request_send_fi_to_I_1J_tmp = CalculateIndexRequestForSendOrRecvVariableNeighbour(
			variable_send_tmp, to_I_1J, index_send);

		index_request_recv_fi_from_I_1J_tmp = CalculateIndexRequestForSendOrRecvVariableNeighbour(
			variable_send_tmp, to_I_1J, index_recv);
	}

	//Exchange data with subdomain I1J, if exist.
	if(SubDomain[IJ].is_SubDomain_in_direction_I1)
	{
		tag_send_fi_to_I1J_tmp = SubDomain[IJ].CalculateTagForSendVariableToNeighbour(
			variable_send_tmp, to_I1J);

		tag_send_fi_to_I_1J_tmp_from_I1J = SubDomain[I1J].CalculateTagForSendVariableToNeighbour(
			variable_send_tmp, to_I_1J);


		index_request_send_fi_to_I1J_tmp = CalculateIndexRequestForSendOrRecvVariableNeighbour(
			variable_send_tmp, to_I1J, index_send);

		index_request_recv_fi_from_I1J_tmp = CalculateIndexRequestForSendOrRecvVariableNeighbour(
			variable_send_tmp, to_I1J, index_recv);
	}

	//Exchange data with subdomain IJ_1, if exist.
	if(SubDomain[IJ].is_SubDomain_in_direction_J_1)
	{
		tag_send_fi_to_IJ_1_tmp = SubDomain[IJ].CalculateTagForSendVariableToNeighbour(
			variable_send_tmp, to_IJ_1);

		tag_send_fi_to_IJ1_tmp_from_IJ_1 = SubDomain[IJ_1].CalculateTagForSendVariableToNeighbour(
			variable_send_tmp, to_IJ1);


		index_request_send_fi_to_IJ_1_tmp = CalculateIndexRequestForSendOrRecvVariableNeighbour(
			variable_send_tmp, to_IJ_1, index_send);

		index_request_recv_fi_from_IJ_1_tmp = CalculateIndexRequestForSendOrRecvVariableNeighbour(
			variable_send_tmp, to_IJ_1, index_recv);
	}

	//Exchange data with subdomain IJ1, if exist.
	if(SubDomain[IJ].is_SubDomain_in_direction_J1)
	{
		tag_send_fi_to_IJ1_tmp = SubDomain[IJ].CalculateTagForSendVariableToNeighbour(
			variable_send_tmp, to_IJ1);

		tag_send_fi_to_IJ_1_tmp_from_IJ1 = SubDomain[IJ1].CalculateTagForSendVariableToNeighbour(
			variable_send_tmp, to_IJ_1);


		index_request_send_fi_to_IJ1_tmp = CalculateIndexRequestForSendOrRecvVariableNeighbour(
			variable_send_tmp, to_IJ1, index_send);

		index_request_recv_fi_from_IJ1_tmp = CalculateIndexRequestForSendOrRecvVariableNeighbour(
			variable_send_tmp, to_IJ1, index_recv);
	}
	//End of Define tags and index of requests for data exchange.


	//Exchange data with neighbourhoods subdomains.
	//Exchange data with subdomain I_1J, if exist.
	if(SubDomain[IJ].is_SubDomain_in_direction_I_1)
	{
		//Here is count is 1, because is send only one
		//element of type send_to_I_1J, wich contain all information.
		MPI_Isend(&massive[Nx_massive_limits[0] + 0 * Nx], 1, send_to_I_1J,
			SubDomain[IJ].I_1J, tag_send_fi_to_I_1J_tmp,
			MPI_COMM_WORLD, &request[index_request_send_fi_to_I_1J_tmp]);


		//if(is_created_MPI_Datatype_send_to_I_1J)
		//{
		//	//This is not make, because the MPI_Datatypes for this data exchange
		//	//are the same for all computations and is defined once.
		//	//Finally free datatype - to not make error from too many MPI objects.
		//	MPI_Type_free(&send_to_I_1J);
		//	is_created_MPI_Datatype_send_to_I_1J = false;
		//}



		//Here is count is 1, because is receive only one
		//element of type recv_from_I_1J, wich contain all information.
		MPI_Irecv(&massive[Nx_massive_limits[0] - blocklength_recv_from_I_1J + 0 * Nx],
			1, recv_from_I_1J,
			SubDomain[IJ].I_1J, tag_send_fi_to_I1J_tmp_from_I_1J,
			MPI_COMM_WORLD, &request[index_request_recv_fi_from_I_1J_tmp]);


		//if(is_created_MPI_Datatype_recv_from_I_1J)
		//{
		//	//This is not make, because the MPI_Datatypes for this data exchange
		//	//are the same for all computations and is defined once.
		//	//Finally free datatype - to not make error from too many MPI objects.
		//	MPI_Type_free(&recv_from_I_1J);
		//	is_created_MPI_Datatype_recv_from_I_1J = false;
		//}

	}


	//Exchange data with subdomain I1J, if exist.
	if(SubDomain[IJ].is_SubDomain_in_direction_I1)
	{
		//Here is count is 1, because is send only one
		//element of type send_to_I_1J, wich contain all information.
		MPI_Isend(&massive[Nx_massive_limits[1] - blocklength_send_to_I1J + 0 * Nx],
			1, send_to_I1J,
			SubDomain[IJ].I1J, tag_send_fi_to_I1J_tmp,
			MPI_COMM_WORLD, &request[index_request_send_fi_to_I1J_tmp]);


		//if(is_created_MPI_Datatype_send_to_I1J)
		//{
		//	//This is not make, because the MPI_Datatypes for this data exchange
		//	//are the same for all computations and is defined once.
		//	//Finally free datatype - to not make error from too many MPI objects.
		//	MPI_Type_free(&send_to_I1J);
		//	is_created_MPI_Datatype_send_to_I1J = false;
		//}



		//Here is count is 1, because is receive only one
		//element of type recv_from_I1J, wich contain all information.
		MPI_Irecv(&massive[Nx_massive_limits[1] + 0 * Nx], 1, recv_from_I1J,
			SubDomain[IJ].I1J, tag_send_fi_to_I_1J_tmp_from_I1J,
			MPI_COMM_WORLD, &request[index_request_recv_fi_from_I1J_tmp]);

		//if(is_created_MPI_Datatype_recv_from_I1J)
		//{
		//	//This is not make, because the MPI_Datatypes for this data exchange
		//	//are the same for all computations and is defined once.
		//	//Finally free datatype - to not make error from too many MPI objects.
		//	MPI_Type_free(&recv_from_I1J);
		//	is_created_MPI_Datatype_recv_from_I1J = false;
		//}
	}


	//Exchange data with subdomain IJ_1, if exist.
	if(SubDomain[IJ].is_SubDomain_in_direction_J_1)
	{
		//Here the data are continues and is not need to make derived datatype
		//Number of ellement, which will be send to subdomain IJ_1
		N_send_to_IJ_1 = BlockLengthSwap_J1 * Nx;

		MPI_Isend(&massive[0 + Ny_massive_limits[0] * Nx], N_send_to_IJ_1, MPI_DOUBLE,
			SubDomain[IJ].IJ_1, tag_send_fi_to_IJ_1_tmp,
			MPI_COMM_WORLD, &request[index_request_send_fi_to_IJ_1_tmp]);

		//Number of ellement, which will be recv from subdomain IJ_1
		N_recv_from_IJ_1 = BlockLengthSwap_J_1 * Nx;

		MPI_Irecv(&massive[0 + (Ny_massive_limits[0] - BlockLengthSwap_J_1) * Nx],
			N_recv_from_IJ_1, MPI_DOUBLE,
			SubDomain[IJ].IJ_1, tag_send_fi_to_IJ1_tmp_from_IJ_1,
			MPI_COMM_WORLD, &request[index_request_recv_fi_from_IJ_1_tmp]);
	}


	//Exchange data with subdomain IJ1, if exist.
	if(SubDomain[IJ].is_SubDomain_in_direction_J1)
	{
		//Number of ellement, which will be send to subdomain IJ1
		N_send_to_IJ1 = BlockLengthSwap_J_1 * Nx;

		MPI_Isend(&massive[0 + (Ny_massive_limits[1] - BlockLengthSwap_J_1) * Nx],
			N_send_to_IJ1, MPI_DOUBLE,
			SubDomain[IJ].IJ1, tag_send_fi_to_IJ1_tmp,
			MPI_COMM_WORLD, &request[index_request_send_fi_to_IJ1_tmp]);

		//Number of ellement, which will be recv from subdomain IJ1
		N_recv_from_IJ1 = BlockLengthSwap_J1 * Nx;

		MPI_Irecv(&massive[0 + Ny_massive_limits[1] * Nx],
			N_recv_from_IJ1, MPI_DOUBLE,
			SubDomain[IJ].IJ1, tag_send_fi_to_IJ_1_tmp_from_IJ1,
			MPI_COMM_WORLD, &request[index_request_recv_fi_from_IJ1_tmp]);
	}
	//End of Exchange data with neighbourhoods subdomains.

}

template <typename T_massive>
inline void WaitToCompleteDataExchangeWithNeighbourhoodsSubDomains_MPI(
	T_massive massive[])
{
	//Define index of requests for data exchange.
	int index_request_send_fi_to_I_1J_tmp, index_request_recv_fi_from_I_1J_tmp,
		index_request_send_fi_to_I1J_tmp, index_request_recv_fi_from_I1J_tmp,
		index_request_send_fi_to_IJ_1_tmp, index_request_recv_fi_from_IJ_1_tmp,
		index_request_send_fi_to_IJ1_tmp, index_request_recv_fi_from_IJ1_tmp;

	int variable_send_tmp;

	//Ininicialize index of requests
	index_request_send_fi_to_I_1J_tmp = 0;
	index_request_recv_fi_from_I_1J_tmp = 0;
	index_request_send_fi_to_I1J_tmp = 0;
	index_request_recv_fi_from_I1J_tmp = 0;
	index_request_send_fi_to_IJ_1_tmp = 0;
	index_request_recv_fi_from_IJ_1_tmp = 0;
	index_request_send_fi_to_IJ1_tmp = 0;
	index_request_recv_fi_from_IJ1_tmp = 0;


	//Calculate index of requests in depend of exchanged variable.
	//The variable will be recognized in depend of address of massive.
	if(massive == d_p_c_12)
	{
		//Exchange array is d_p_c_12
		variable_send_tmp = d_p_c_12_send;
	}
	else if(massive == apx)
	{
		//Exchange array is apx
		variable_send_tmp = apx_send;
	}
	else if(massive == bpx)
	{
		//Exchange array is bpx
		variable_send_tmp = bpx_send;
	}
	else if(massive == d_p_c_34)
	{
		//Exchange array is d_p_c_34
		variable_send_tmp = d_p_c_34_send;
	}
	else if(massive == apy)
	{
		//Exchange array is apy
		variable_send_tmp = apy_send;
	}
	else if(massive == bpy)
	{
		//Exchange array is bpy
		variable_send_tmp = bpy_send;
	}
	else if(massive == p)
	{
		//Exchange array is p
		variable_send_tmp = p_send;
	}
	else if(massive == Temper)
	{
		//Exchange array is Temper
		variable_send_tmp = Temper_send;
	}
	else if(massive == u)
	{
		//Exchange array is u
		variable_send_tmp = u_send;
	}
	else if(massive == v)
	{
		//Exchange array is v
		variable_send_tmp = v_send;
	}



	//Exchange variable_send_tmp
	//Exchange data with subdomain I_1J, if exist.
	if(SubDomain[IJ].is_SubDomain_in_direction_I_1)
	{
		index_request_send_fi_to_I_1J_tmp = CalculateIndexRequestForSendOrRecvVariableNeighbour(
			variable_send_tmp, to_I_1J, index_send);

		index_request_recv_fi_from_I_1J_tmp = CalculateIndexRequestForSendOrRecvVariableNeighbour(
			variable_send_tmp, to_I_1J, index_recv);
	}

	//Exchange data with subdomain I1J, if exist.
	if(SubDomain[IJ].is_SubDomain_in_direction_I1)
	{
		index_request_send_fi_to_I1J_tmp = CalculateIndexRequestForSendOrRecvVariableNeighbour(
			variable_send_tmp, to_I1J, index_send);

		index_request_recv_fi_from_I1J_tmp = CalculateIndexRequestForSendOrRecvVariableNeighbour(
			variable_send_tmp, to_I1J, index_recv);
	}

	//Exchange data with subdomain IJ_1, if exist.
	if(SubDomain[IJ].is_SubDomain_in_direction_J_1)
	{
		index_request_send_fi_to_IJ_1_tmp = CalculateIndexRequestForSendOrRecvVariableNeighbour(
			variable_send_tmp, to_IJ_1, index_send);

		index_request_recv_fi_from_IJ_1_tmp = CalculateIndexRequestForSendOrRecvVariableNeighbour(
			variable_send_tmp, to_IJ_1, index_recv);
	}

	//Exchange data with subdomain IJ1, if exist.
	if(SubDomain[IJ].is_SubDomain_in_direction_J1)
	{
		index_request_send_fi_to_IJ1_tmp = CalculateIndexRequestForSendOrRecvVariableNeighbour(
			variable_send_tmp, to_IJ1, index_send);

		index_request_recv_fi_from_IJ1_tmp = CalculateIndexRequestForSendOrRecvVariableNeighbour(
			variable_send_tmp, to_IJ1, index_recv);
	}
	//End of Define index of requests for data exchange.


	//Wait to recieve data frrom neighbourhoods subdomains.
	//Wait to complete exchange data with subdomain I_1J, if exist.
	if(SubDomain[IJ].is_SubDomain_in_direction_I_1)
	{
		MPI_Wait(&request[index_request_send_fi_to_I_1J_tmp],
			&status[index_request_send_fi_to_I_1J_tmp]);

		MPI_Wait(&request[index_request_recv_fi_from_I_1J_tmp],
			&status[index_request_recv_fi_from_I_1J_tmp]);
	}

	//Wait to complete exchange data with domain I1J, if exist.
	if(SubDomain[IJ].is_SubDomain_in_direction_I1)
	{
		MPI_Wait(&request[index_request_send_fi_to_I1J_tmp],
			&status[index_request_send_fi_to_I1J_tmp]);

		MPI_Wait(&request[index_request_recv_fi_from_I1J_tmp],
			&status[index_request_recv_fi_from_I1J_tmp]);
	}

	//Wait to complete exchange data with subdomain IJ_1, if exist.
	if(SubDomain[IJ].is_SubDomain_in_direction_J_1)
	{
		MPI_Wait(&request[index_request_send_fi_to_IJ_1_tmp],
			&status[index_request_send_fi_to_IJ_1_tmp]);

		MPI_Wait(&request[index_request_recv_fi_from_IJ_1_tmp],
			&status[index_request_recv_fi_from_IJ_1_tmp]);
	}

	//Wait to complete exchange data with subdomain IJ1, if exist.
	if(SubDomain[IJ].is_SubDomain_in_direction_J1)
	{
		MPI_Wait(&request[index_request_send_fi_to_IJ1_tmp],
			&status[index_request_send_fi_to_IJ1_tmp]);

		MPI_Wait(&request[index_request_recv_fi_from_IJ1_tmp],
			&status[index_request_recv_fi_from_IJ1_tmp]);
	}
	//End of Wait to recieve data frrom neighbourhood subdomains.


}


inline void StartSendDataForConvergenceCriterionsFromSubDomainsFor_pT_MPI(void)
{
	if(IJ == 0)
	{
		//Receive max_residual_in_u, max_residual_in_v, max_residual_in_p
		//and max_residual_in_T from all subdomains.
		//For this process (IJ == 0)
		max_residual_in_p_all_subdomains[0] = max_residual_in_p[Iter];
		max_residual_in_T_all_subdomains[0] = max_residual_in_T[Iter];

		//Recv data from all subdomains
		//Must be defined internal counter, because the part of calculation
		//continiue parrallel to this procces (exchanging data).
		int counter_tmp;
		for(counter_tmp = 1; counter_tmp < N_SubDomains_a; counter_tmp++)
		{
			//Recv max_residual_in_p from process counter_tmp
			MPI_Irecv(&max_residual_in_p_all_subdomains[counter_tmp], 1, MPI_DOUBLE,
				counter_tmp,
				SubDomain[counter_tmp].CalculateTagForSendVariableToIJ0(max_residual_in_p_send_to_IJ0),
				MPI_COMM_WORLD,
				&request_max_residual_in_p[counter_tmp]);

			//Recv max_residual_in_T from process counter_tmp
			MPI_Irecv(&max_residual_in_T_all_subdomains[counter_tmp], 1, MPI_DOUBLE,
				counter_tmp,
				SubDomain[counter_tmp].CalculateTagForSendVariableToIJ0(max_residual_in_T_send_to_IJ0),
				MPI_COMM_WORLD,
				&request_max_residual_in_T[counter_tmp]);
		}

	}
	else
	{
		//Send max_residual_in_p and max_residual_in_T for this subdomain to
		//process IJ == 0.
		//Send max_residual_in_p to IJ == 0
		MPI_Isend(&max_residual_in_p[Iter], 1, MPI_DOUBLE,
			0, SubDomain[IJ].CalculateTagForSendVariableToIJ0(max_residual_in_p_send_to_IJ0),
			MPI_COMM_WORLD,
			&request[CalculateIndexRequestForSendVariableToIJ0(max_residual_in_p_send_to_IJ0)]);

		//Send max_residual_in_T to IJ == 0
		MPI_Isend(&max_residual_in_T[Iter], 1, MPI_DOUBLE,
			0, SubDomain[IJ].CalculateTagForSendVariableToIJ0(max_residual_in_T_send_to_IJ0),
			MPI_COMM_WORLD,
			&request[CalculateIndexRequestForSendVariableToIJ0(max_residual_in_T_send_to_IJ0)]);
	}

}


inline void StartSendDataForConvergenceCriterionsFromSubDomainsFor_u_MPI(void)
{
	if(IJ == 0)
	{
		//Receive max_residual_in_u all subdomains.
		//For this process (IJ == 0)
		max_residual_in_u_all_subdomains[0] = max_residual_in_u;

		//Recv data from all subdomains
		//Must be defined internal counter, because the part of calculation
		//continiue parrallel to this procces (exchanging data).
		int counter_tmp;
		for(counter_tmp = 1; counter_tmp < N_SubDomains_a; counter_tmp++)
		{
			//Recv max_residual_in_u from process counter_tmp
			MPI_Irecv(&max_residual_in_u_all_subdomains[counter_tmp], 1, MPI_DOUBLE,
				counter_tmp,
				SubDomain[counter_tmp].CalculateTagForSendVariableToIJ0(max_residual_in_u_send_to_IJ0),
				MPI_COMM_WORLD,
				&request_max_residual_in_u[counter_tmp]);
		}

	}
	else
	{
		//Send max_residual_in_ufor this subdomain to process IJ == 0.
		//Send max_residual_in_u to IJ == 0
		MPI_Isend(&max_residual_in_u, 1, MPI_DOUBLE,
			0, SubDomain[IJ].CalculateTagForSendVariableToIJ0(max_residual_in_u_send_to_IJ0),
			MPI_COMM_WORLD,
			&request[CalculateIndexRequestForSendVariableToIJ0(max_residual_in_u_send_to_IJ0)]);
	}

}


inline void StartSendDataForConvergenceCriterionsFromSubDomainsFor_v_MPI(void)
{
	if(IJ == 0)
	{
		//Receive max_residual_in_v from all subdomains.
		//For this process (IJ == 0)
		max_residual_in_v_all_subdomains[0] = max_residual_in_v;

		//Recv data from all subdomains
		//Must be defined internal counter, because the part of calculation
		//continiue parrallel to this procces (exchanging data).
		int counter_tmp;
		for(counter_tmp = 1; counter_tmp < N_SubDomains_a; counter_tmp++)
		{
			//Recv max_residual_in_v from process counter_tmp
			MPI_Irecv(&max_residual_in_v_all_subdomains[counter_tmp], 1, MPI_DOUBLE,
				counter_tmp,
				SubDomain[counter_tmp].CalculateTagForSendVariableToIJ0(max_residual_in_v_send_to_IJ0),
				MPI_COMM_WORLD,
				&request_max_residual_in_v[counter_tmp]);
		}

	}
	else
	{
		//Send max_residual_in_v for this subdomain to process IJ == 0.
		//Send max_residual_in_v to IJ == 0
		MPI_Isend(&max_residual_in_v, 1, MPI_DOUBLE,
			0, SubDomain[IJ].CalculateTagForSendVariableToIJ0(max_residual_in_v_send_to_IJ0),
			MPI_COMM_WORLD,
			&request[CalculateIndexRequestForSendVariableToIJ0(max_residual_in_v_send_to_IJ0)]);
	}

}


inline void WaitToCompleteRecvDataForConvergenceCriterionsFromSubDomains_for_pT_MPI(void)
{
	if(IJ == 0)
	{
		//Must be defined internal counter, because the part of calculation
		//continiue parrallel to this procces
		int counter_tmp;
		for(counter_tmp = 1; counter_tmp < N_SubDomains_a; counter_tmp++)
		{
			//Wait to complete receive data for max_residual_in_p from process counter_tmp
			MPI_Wait(&request_max_residual_in_p[counter_tmp],
				&status_max_residual_in_p[counter_tmp]);

			//Wait to complete receive data for max_residual_in_T from process counter_tmp
			MPI_Wait(&request_max_residual_in_T[counter_tmp],
				&status_max_residual_in_T[counter_tmp]);
		}

	}
	else
	{
		//Wait to complete sending of data for max_residual_in_p to IJ == 0
		MPI_Wait(&request[CalculateIndexRequestForSendVariableToIJ0(max_residual_in_p_send_to_IJ0)],
			&status[CalculateIndexRequestForSendVariableToIJ0(max_residual_in_p_send_to_IJ0)]);

		//Wait to complete sending of data for max_residual_in_T to IJ == 0
		MPI_Wait(&request[CalculateIndexRequestForSendVariableToIJ0(max_residual_in_T_send_to_IJ0)],
			&status[CalculateIndexRequestForSendVariableToIJ0(max_residual_in_T_send_to_IJ0)]);
	}

}


inline void WaitToCompleteRecvDataForConvergenceCriterionsFromSubDomains_for_uv_MPI(void)
{
	if(IJ == 0)
	{
		//Must be defined internal counter, because the part of calculation
		//continiue parrallel to this procces
		int counter_tmp;
		for(counter_tmp = 1; counter_tmp < N_SubDomains_a; counter_tmp++)
		{
			//Wait to complete receive data for max_residual_in_u from process counter_tmp
			MPI_Wait(&request_max_residual_in_u[counter_tmp],
				&status_max_residual_in_u[counter_tmp]);

			//Wait to complete receive data for max_residual_in_v from process counter_tmp
			MPI_Wait(&request_max_residual_in_v[counter_tmp],
				&status_max_residual_in_v[counter_tmp]);
		}

	}
	else
	{
		//Wait to complete sending of data for max_residual_in_u to IJ == 0
		MPI_Wait(&request[CalculateIndexRequestForSendVariableToIJ0(max_residual_in_u_send_to_IJ0)],
			&status[CalculateIndexRequestForSendVariableToIJ0(max_residual_in_u_send_to_IJ0)]);

		//Wait to complete sending of data for max_residual_in_v to IJ == 0
		MPI_Wait(&request[CalculateIndexRequestForSendVariableToIJ0(max_residual_in_v_send_to_IJ0)],
			&status[CalculateIndexRequestForSendVariableToIJ0(max_residual_in_v_send_to_IJ0)]);
	}

}


inline void	FindMaximumResudualsFromAllSubDomains_for_pT_MPI(void)
{
	//Find maximun residual from all SubDomains
	//Must be defined internal counter, because the part of calculation
	//continiue parrallel to this procces
	int counter_tmp;
	for(counter_tmp = 0; counter_tmp < N_SubDomains_a; counter_tmp++)
	{
		//Search for maximum value of max_residual_in_p from all subdomains
		if(max_residual_in_p[Iter] < max_residual_in_p_all_subdomains[counter_tmp])
		{
			//There is already value for max_residual_in_p from this subdomain IJ == 0
			max_residual_in_p[Iter] = max_residual_in_p_all_subdomains[counter_tmp];
		}

		//Search for maximum value of max_residual_in_T from all subdomains
		if(max_residual_in_T[Iter] < max_residual_in_T_all_subdomains[counter_tmp])
		{
			//There is already value for max_residual_in_T from this subdomain IJ == 0
			max_residual_in_T[Iter] = max_residual_in_T_all_subdomains[counter_tmp];
		}

	}

}

inline void	FindMaximumResudualsFromAllSubDomains_for_uv_MPI(void)
{
	//Find maximun residual from all SubDomains
	//Must be defined internal counter, because the part of calculation
	//continiue parrallel to this procces
	int counter_tmp;
	for(counter_tmp = 0; counter_tmp < N_SubDomains_a; counter_tmp++)
	{
		//Even in one subdomain the calculation are not completed the calculations
		//are continue in all sundomains.
		//Search for maximum value of max_residual_in_u from all subdomains
		if(max_residual_in_u < max_residual_in_u_all_subdomains[counter_tmp])
		{
			//There is already value for max_residual_in_u from this subdomain IJ == 0
			max_residual_in_u = max_residual_in_u_all_subdomains[counter_tmp];
		}

		//Search for maximum value of max_residual_in_v from all subdomains
		if(max_residual_in_v < max_residual_in_v_all_subdomains[counter_tmp])
		{
			//There is already value for max_residual_in_v from this subdomain IJ == 0
			max_residual_in_v = max_residual_in_v_all_subdomains[counter_tmp];
		}
	}

}

inline void StartSendDataForConvergenceCriterionsToSubDomains_for_pT_MPI(void)
{
	if(IJ == 0)
	{
		continue_iter_p_c_int = continue_iter_p_c;

		//continue_iter_p_c_int will be send to all subdomains.
		//Must be defined internal counter, because the part of calculation
		//continiue parrallel to this procces
		int counter_tmp;
		for(counter_tmp = 1; counter_tmp < N_SubDomains_a; counter_tmp++)
		{
			//Send continue_iter_int to subdomain counter_tmp
			//Here is used request_max_residual_in_u to check for completing data send.
			MPI_Isend(&continue_iter_p_c_int, 1, MPI_INT,
				counter_tmp,
				SubDomain[counter_tmp].CalculateTagForSendVariableFromIJ0(continue_iter_p_c_int_send_from_IJ0),
				MPI_COMM_WORLD, &request_max_residual_in_p[counter_tmp]);
		}

	}
	else
	{
		//Receive data for continue_iter_int from IJ == 0
		MPI_Irecv(&continue_iter_p_c_int, 1, MPI_INT,
			0, SubDomain[IJ].CalculateTagForSendVariableFromIJ0(continue_iter_p_c_int_send_from_IJ0),
			MPI_COMM_WORLD,
			&request[CalculateIndexRequestForRecvVariableFromIJ0(continue_iter_p_c_int_send_from_IJ0)]);
	}

}


inline void WaitToCompleteRecvDataForConvergenceCriterionsFromIJ_0_for_pT_MPI(void)
{
	if(IJ == 0)
	{
		//Must be defined internal counter, because the part of calculation
		//continiue parrallel to this procces
		int counter_tmp;
		for(counter_tmp = 1; counter_tmp < N_SubDomains_a; counter_tmp++)
		{
			//Wait to complete send data for continue_iter_int to process counter_tmp
			MPI_Wait(&request_max_residual_in_p[counter_tmp],
				&status_max_residual_in_p[counter_tmp]);
		}

	}
	else
	{
		//Wait to complete receiving data for continue_iter_int from process IJ == 0
		MPI_Wait(&request[CalculateIndexRequestForRecvVariableFromIJ0(continue_iter_p_c_int_send_from_IJ0)],
			&status[CalculateIndexRequestForRecvVariableFromIJ0(continue_iter_p_c_int_send_from_IJ0)]);

		//After recv is complete, can be solved continue_iter_p_c.
		continue_iter_p_c = continue_iter_p_c_int;
	}
}


inline void StartSendDataForConvergenceCriterionsToSubDomains_for_it_loop_MPI(void)
{
	if(IJ == 0)
	{
		continue_iter_int = continue_iter;
		countinue_time_step_int = countinue_time_step;

		//continue_iter_tmp_int and countinue_time_step_int_tmp will be send to
		//all subdomains.
		//Must be defined internal counter, because the part of calculation
		//continiue parrallel to this procces
		int counter_tmp;
		for(counter_tmp = 1; counter_tmp < N_SubDomains_a; counter_tmp++)
		{
			//Send continue_iter_int to subdomain counter_tmp
			//Here is used request_max_residual_in_u to check for completing data send.
			MPI_Isend(&continue_iter_int, 1, MPI_INT,
				counter_tmp,
				SubDomain[counter_tmp].CalculateTagForSendVariableFromIJ0(continue_iter_int_send_from_IJ0),
				MPI_COMM_WORLD, &request_max_residual_in_u[counter_tmp]);

			//Send countinue_time_step_int to subdomain counter_tmp
			//Here is used request_max_residual_in_v to check for completing data send.
			MPI_Isend(&countinue_time_step_int, 1, MPI_INT,
				counter_tmp,
				SubDomain[counter_tmp].CalculateTagForSendVariableFromIJ0(continue_time_step_int_send_from_IJ0),
				MPI_COMM_WORLD, &request_max_residual_in_v[counter_tmp]);
		}

	}
	else
	{
		//Receive data for continue_iter_int from IJ == 0
		MPI_Irecv(&continue_iter_int, 1, MPI_INT,
			0, SubDomain[IJ].CalculateTagForSendVariableFromIJ0(continue_iter_int_send_from_IJ0),
			MPI_COMM_WORLD,
			&request[CalculateIndexRequestForRecvVariableFromIJ0(continue_iter_int_send_from_IJ0)]);

		//Receive data for countinue_time_step_int from IJ == 0
		MPI_Irecv(&countinue_time_step_int, 1, MPI_INT,
			0, SubDomain[IJ].CalculateTagForSendVariableFromIJ0(continue_time_step_int_send_from_IJ0),
			MPI_COMM_WORLD,
			&request[CalculateIndexRequestForRecvVariableFromIJ0(continue_time_step_int_send_from_IJ0)]);
	}

}

inline void WaitToCompleteRecvDataForConvergenceCriterionsFromIJ_0_for_it_loop_MPI(void)
{
	if(IJ == 0)
	{
		//Must be defined internal counter, because the part of calculation
		//continiue parrallel to this procces
		int counter_tmp;
		for(counter_tmp = 1; counter_tmp < N_SubDomains_a; counter_tmp++)
		{
			//Wait to complete send data for continue_iter_int to process counter_tmp
			MPI_Wait(&request_max_residual_in_u[counter_tmp],
				&status_max_residual_in_u[counter_tmp]);

			//Wait to complete send data for countinue_time_step_int to process counter_tmp
			MPI_Wait(&request_max_residual_in_v[counter_tmp],
				&status_max_residual_in_v[counter_tmp]);
		}

	}
	else
	{
		//Wait to complete receiving data for continue_iter_int from process IJ == 0
		MPI_Wait(&request[CalculateIndexRequestForRecvVariableFromIJ0(continue_iter_int_send_from_IJ0)],
			&status[CalculateIndexRequestForRecvVariableFromIJ0(continue_iter_int_send_from_IJ0)]);

		//Wait to complete receiving data for countinue_time_step_int from process IJ == 0
		MPI_Wait(&request[CalculateIndexRequestForRecvVariableFromIJ0(continue_time_step_int_send_from_IJ0)],
			&status[CalculateIndexRequestForRecvVariableFromIJ0(continue_time_step_int_send_from_IJ0)]);

		//After recv is complete, can be solved continue_iter and countinue_time_step.
		continue_iter = continue_iter_int;
		countinue_time_step = countinue_time_step_int;
	}
}


template <typename T_FileName, typename T_massive>
void Write_massive_entire_computational_domain_MPI(
	const T_FileName& FileNameWithoutExtension,
	T_massive massive[], const unsigned int& massive_type)
{
	if(IJ == 0)
	{
		//Make temporal massive, where data will be stored and
		//after that will be written into the file.
		T_massive * massive_entire_computational_domain_tmp;
		massive_entire_computational_domain_tmp = new T_massive [Na_entire_computational_domain];

		//copy data from massive, because here IJ == 0
		for(j = 0; j < Ny; j++)
		{
			for(i = 0; i < Nx; i++)
			{
				massive_entire_computational_domain_tmp[i + j * Nx_entire_computational_domain] = massive[i + j * Nx];
			}
		}

		MPI_Status status_tmp;

		////Index_of_addition is index from, where will be added information
		////from other threads.
		//unsigned int Index_of_addition;
		//Index_of_addition = Nx_others[1] * Ny_others[1];

		//Define variables needed to make derived data type
		unsigned int N_recv_from_counter,
			blocklength_recv_from_counter,
			stride_recv_from_counter;


		//Number of elements, blocklength and stride which will be receive from subdomain I_1J
		int N_recv_from_counter_tmp, blocklength_recv_from_counter_tmp,
			stride_recv_from_counter_tmp;


		//Must be defined internal counter, because the part of calculation
		//continue parallel to this process
		int counter_tmp;
		//Recv data from SubDomains 0 < IJ
		for(counter_tmp = 1; counter_tmp < N_SubDomains_a; counter_tmp++)
		{
			////Receive all array and put it in correct place
			////Here have to be created Derived Datatipe to transfer the data.
			////Create derived datatype MPI_TYPE_VECTOR to receive data from subdomain counter
			//MPI_Datatype recv_from_counter_tmp;

			//Number of all elements to receive
			N_recv_from_counter_tmp = SubDomain[counter_tmp].Na;

			//Calculate blocklength
			blocklength_recv_from_counter_tmp = SubDomain[counter_tmp].Nx;

			////Stide is Nx_entire_computational_domain
			//stride_recv_from_counter_tmp = Nx_entire_computational_domain;

			//if(massive_type == double_type)
			//{
			//	//If recv data are from type double
			//	MPI_Type_vector(
			//		N_recv_from_counter_tmp, blocklength_recv_from_counter_tmp,
			//		stride_recv_from_counter_tmp,
			//		MPI_DOUBLE, &recv_from_counter_tmp);
			//}
			//else if(massive_type == unsigned_int_type)
			//{
			//	//If recv data are from type unsigned int
			//	MPI_Type_vector(
			//		N_recv_from_counter_tmp, blocklength_recv_from_counter_tmp,
			//		stride_recv_from_counter_tmp,
			//		MPI_UNSIGNED, &recv_from_counter_tmp);
			//}


			//MPI_Type_commit(&recv_from_counter_tmp);
			////End of Create derived datatype MPI_TYPE_VECTOR to receive data from subdomain counter


			////Start write to massive fi_all from index Index_of_addition.

			////Calculate Index_of_addition
			//Index_of_addition = SubDomain[counter_tmp].i_b - BlockLengthSwap_I_1 * SubDomain[counter_tmp].is_SubDomain_in_direction_I_1
			//	+ (SubDomain[counter_tmp].j_b - BlockLengthSwap_J_1 * SubDomain[counter_tmp].is_SubDomain_in_direction_J_1) * Nx_entire_computational_domain;


			////When from subdomains the data are send with MPI_Isend:
			//MPI_Irecv(&massive_entire_computational_domain_tmp[Index_of_addition],
			//	1, recv_from_counter_tmp,
			//	counter_tmp, SubDomain[counter_tmp].tag_send_fi0_to_I_1J,
			//	MPI_COMM_WORLD, &request[index_request_send_fi0_to_I_1J]);

			////Wait to complete MPI_Irecv
			//MPI_Wait(&request[index_request_send_fi0_to_I_1J],
			//	&status[index_request_send_fi0_to_I_1J]);

			//
			////MPI_Status status_tmp;
			//////When from subdomains the data are send with MPI_Isend:
			////MPI_Recv(&massive_entire_computational_domain_tmp[Index_of_addition],
			////	1, recv_from_counter_tmp,
			////	counter_tmp, SubDomain[counter_tmp].tag_send_fi0_to_I_1J,
			////	MPI_COMM_WORLD, &status_tmp);


			////Finally free datatype - to not make error from too many MPI objects.
			//MPI_Type_free(&recv_from_counter_tmp);



			//The above code works only on Windows. It can not put data in right place
			//on Linux. May be this is, because for Windows I use MPI 2.0
			//for Linux MPI 1.2.
			//The important is that this is make slow only saving of data.
			//
			//Make temporal array where data will be received.
			T_massive * massive_recv_from_SubDomain_tmp;
			//This array will be freed in the end of the loop, because is locally defined.
			massive_recv_from_SubDomain_tmp = new T_massive [N_recv_from_counter_tmp];

			//Data are copied to massive_entire_computational_domain_tmp
			//and temporral massive will be deleted.
			//If recv data are from type double
			if(massive_type == double_type)
			{
				//When from subdomains the data are send with MPI_Isend:
				MPI_Irecv(&massive_recv_from_SubDomain_tmp[0],
					N_recv_from_counter_tmp, MPI_DOUBLE,
					counter_tmp,
					SubDomain[counter_tmp].CalculateTagForSendVariableToNeighbour(u_send, to_I_1J),
					MPI_COMM_WORLD,
					&request[CalculateIndexRequestForSendOrRecvVariableNeighbour(u_send, to_I_1J, index_recv)]);
			}
			else if(massive_type == unsigned_int_type)
			{
				//When from subdomains the data are send with MPI_Isend:
				MPI_Irecv(&massive_recv_from_SubDomain_tmp[0],
					N_recv_from_counter_tmp, MPI_UNSIGNED,
					counter_tmp,
					SubDomain[counter_tmp].CalculateTagForSendVariableToNeighbour(u_send, to_I_1J),
					MPI_COMM_WORLD,
					&request[CalculateIndexRequestForSendOrRecvVariableNeighbour(u_send, to_I_1J, index_recv)]);
			}


			//Wait to complete MPI_Irecv
			MPI_Wait(&request[CalculateIndexRequestForSendOrRecvVariableNeighbour(u_send, to_I_1J, index_recv)],
				&status[CalculateIndexRequestForSendOrRecvVariableNeighbour(u_send, to_I_1J, index_recv)]);

			//Temporal local variables
			int i_tmp, j_tmp, node_entire_compuattional_domain_tmp;

			//Where will be copied massive from SunDomain[counter_tmp], index
			//index in axis OX
			int i_addition;
			i_addition = SubDomain[counter_tmp].i_b - BlockLengthSwap_I_1 * SubDomain[counter_tmp].is_SubDomain_in_direction_I_1;

			//index in axis OY
			int j_addition;
			j_addition = SubDomain[counter_tmp].j_b - BlockLengthSwap_J_1 * SubDomain[counter_tmp].is_SubDomain_in_direction_J_1;

			//copy data from massive counter_tmp
			for(j_tmp = BlockLengthSwap_J_1 * SubDomain[counter_tmp].is_SubDomain_in_direction_J_1;
				j_tmp < SubDomain[counter_tmp].Ny - BlockLengthSwap_J1 * SubDomain[counter_tmp].is_SubDomain_in_direction_J1;
				j_tmp++)
			{
				for(i_tmp = BlockLengthSwap_I_1 * SubDomain[counter_tmp].is_SubDomain_in_direction_I_1;
					i_tmp < SubDomain[counter_tmp].Nx - BlockLengthSwap_I1 * SubDomain[counter_tmp].is_SubDomain_in_direction_I1;
					i_tmp++)
				{
					//Node in massive_entire_computational_domain_tmp,
					//where node will be copied.
					node_entire_compuattional_domain_tmp =
						(i_addition + i_tmp)
						+ (j_addition + j_tmp) * Nx_entire_computational_domain;

					massive_entire_computational_domain_tmp[node_entire_compuattional_domain_tmp] =
						massive_recv_from_SubDomain_tmp[i_tmp + j_tmp * SubDomain[counter_tmp].Nx];
				}
			}


			//This array MUST be deleted, because is defined in function.
			//The function not delete arrayes defined in in automaticaly.
			//Delete only local variables (this is make by compiler).
			delete [] massive_recv_from_SubDomain_tmp;

		}

		//WriteMassiveToFile_2D("massive_2D.txt",
		//	massive_entire_computational_domain_tmp,
		//	Nx_entire_computational_domain, Ny_entire_computational_domain);

		string FileNameWithoutExtension_string_tmp = FileNameWithoutExtension;

		//Write massive_entire_computational_domain_tmp to file
		if(ToWriteSolvedDataToBinaryFile)
		{
			string filename_entire = FileNameWithoutExtension_string_tmp + ".bin";

			WriteMassiveToBinaryFile(filename_entire.c_str(),
				massive_entire_computational_domain_tmp, Na_entire_computational_domain);

			if(ToWriteSolvedDataToNewFiles)
			{
				string filename_entire_new = FileNameWithoutExtension_string_tmp + "." + toString(it) + ".bin";

				WriteMassiveToBinaryFile(filename_entire_new.c_str(),
					massive_entire_computational_domain_tmp, Na_entire_computational_domain);

			}


		}
		else
		{
			string filename_entire = FileNameWithoutExtension_string_tmp + ".txt";

			WriteMassiveToFile_1D(filename_entire.c_str(),
				massive_entire_computational_domain_tmp, Na_entire_computational_domain);

			if(ToWriteSolvedDataToNewFiles)
			{
				string filename_entire_new = FileNameWithoutExtension_string_tmp + "." + toString(it) + ".txt";

				WriteMassiveToFile_1D(filename_entire_new.c_str(),
					massive_entire_computational_domain_tmp, Na_entire_computational_domain);

			}


		}


		//This array MUST be deleted, because is defined in function.
		//The function not delete arrayes defined in in automaticaly.
		//Delete only local variables (this is make by compiler).
		delete [] massive_entire_computational_domain_tmp;

	}
	else
	{
		//Isend all array massive
		if(massive_type == double_type)
		{
			//If recv data are from type double
			MPI_Isend(&massive[0], SubDomain[IJ].Na, MPI_DOUBLE,
				0, SubDomain[IJ].CalculateTagForSendVariableToNeighbour(u_send, to_I_1J),
				MPI_COMM_WORLD,
				&request[CalculateIndexRequestForSendOrRecvVariableNeighbour(u_send, to_I_1J, index_send)]);

			////If recv data are from type double
			//MPI_Ssend(&massive[0], SubDomain[IJ].Na, MPI_DOUBLE,
			//	0, SubDomain[IJ].tag_send_fi0_to_I_1J,
			//	MPI_COMM_WORLD);

		}
		else if(massive_type == unsigned_int_type)
		{
			//If recv data are from type unsigned int
			MPI_Isend(&massive[0], SubDomain[IJ].Na, MPI_UNSIGNED,
				0, SubDomain[IJ].CalculateTagForSendVariableToNeighbour(u_send, to_I_1J),
				MPI_COMM_WORLD,
				&request[CalculateIndexRequestForSendOrRecvVariableNeighbour(u_send, to_I_1J, index_send)]);

			////If recv data are from type unsigned int
			//MPI_Ssend(&massive[0], SubDomain[IJ].Na, MPI_UNSIGNED,
			//	0, SubDomain[IJ].tag_send_fi0_to_I_1J,
			//	MPI_COMM_WORLD);
		}


		//Wait to complete MPI_Isend
		MPI_Wait(&request[CalculateIndexRequestForSendOrRecvVariableNeighbour(u_send, to_I_1J, index_send)],
			&status[CalculateIndexRequestForSendOrRecvVariableNeighbour(u_send, to_I_1J, index_send)]);

	}



}

template <typename T_FileName, typename T_massive>
void Read_massive_entire_computational_domain_MPI(
	const T_FileName& FileNameWithoutExtension,
	T_massive massive[], const unsigned int& massive_type)
{
	if(IJ == 0)
	{
		//Make temporal massive, where data will be stored, when is readed from file and
		//after that will be send to other subdomains.
		T_massive * massive_entire_computational_domain_tmp;
		massive_entire_computational_domain_tmp = new T_massive [Na_entire_computational_domain];


		string FileNameWithoutExtension_string_tmp = FileNameWithoutExtension;

		//Read data from file
		if(ToReadSolvedDataFromBinaryFile)
		{
			string filename_entire = FileNameWithoutExtension_string_tmp + ".bin";

			ReadMassiveFromBinaryFile(filename_entire.c_str(),
				massive_entire_computational_domain_tmp, Na_entire_computational_domain);
		}
		else
		{
			string filename_entire = FileNameWithoutExtension_string_tmp + ".txt";

			ReadMassiveFromFile_1D(filename_entire.c_str(),
				massive_entire_computational_domain_tmp, Na_entire_computational_domain);
		}

		//WriteMassiveToFile_2D("massive_entire_computational_domain_2D.txt", massive_entire_computational_domain_tmp, Nx_entire_computational_domain, Ny_entire_computational_domain);

		//Copy data to massive, becouse here IJ == 0
		for(j = 0; j < Ny; j++)
		{
			for(i = 0; i < Nx; i++)
			{
				massive[i + j * Nx] = massive_entire_computational_domain_tmp[i + j * Nx_entire_computational_domain];
			}
		}

		////Write massive to file
		//WriteMassiveToFile_2D("massive_2D.txt", massive, Nx, Ny);

		MPI_Status status_tmp;

		//Index_of_addition is index from, where will be send information to other subdomains.
		unsigned int Index_of_addition;

		Index_of_addition = Nx_others[1] * Ny_others[1];

		//Number of elements, blocklength and stride, which will be send to subdomain I_1J
		int N_send_to_counter_tmp, blocklength_send_to_counter_tmp,
			stride_send_to_counter_tmp;


		//Must be defined internal counter, because the part of calculation
		//continiue parrallel to this procces
		int counter_tmp;
		//Send data to SubDomains 0 < IJ
		for(counter_tmp = 1; counter_tmp < N_SubDomains_a; counter_tmp++)
		{
			//Send all array and put it in correct place

			//Number of all elements to send
			N_send_to_counter_tmp = SubDomain[counter_tmp].Na;

			//Calculate blocklength
			blocklength_send_to_counter_tmp = SubDomain[counter_tmp].Nx;

			//Stide is Nx_entire_computational_domain
			stride_send_to_counter_tmp = Nx_entire_computational_domain;


			//Make temporal array, which will be send to subdomain counter_tmp.
			T_massive * massive_send_to_SubDomain_tmp;
			//This array will be freed in the end of loop, because is locally defined.
			massive_send_to_SubDomain_tmp = new T_massive [N_send_to_counter_tmp];

			//Temporal local variables
			int i_tmp, j_tmp, node_entire_compuattional_domain_tmp;

			//Where will be copied part of massive_entire_computational_domain_tmp to massive
			//for SunDomain[counter_tmp], index index in axis OX
			int i_addition;
			i_addition = SubDomain[counter_tmp].i_b - BlockLengthSwap_I_1 * SubDomain[counter_tmp].is_SubDomain_in_direction_I_1;

			//index in axis OY
			int j_addition;
			j_addition = SubDomain[counter_tmp].j_b - BlockLengthSwap_J_1 * SubDomain[counter_tmp].is_SubDomain_in_direction_J_1;

			//copy data from massive_entire_computational_domain_tmp to massive
			//for SunDomain[counter_tmp].
			for(j_tmp = 0; j_tmp < SubDomain[counter_tmp].Ny; j_tmp++)
			{
				for(i_tmp = 0; i_tmp < SubDomain[counter_tmp].Nx; i_tmp++)
				{
					//Node in massive_entire_computational_domain_tmp,
					//from where node wil be copied.
					node_entire_compuattional_domain_tmp =
						(i_tmp + i_addition)
						+ (j_tmp + j_addition) * Nx_entire_computational_domain;

					massive_send_to_SubDomain_tmp[i_tmp + j_tmp * SubDomain[counter_tmp].Nx] =
						massive_entire_computational_domain_tmp[node_entire_compuattional_domain_tmp];
				}
			}

			////Write massive to file
			//WriteMassiveToFile_2D("massive_2D.txt", massive_send_to_SubDomain_tmp,
			//	SubDomain[counter_tmp].Nx, SubDomain[counter_tmp].Ny);

			//Send massive_send_to_SubDomain_tmp to SubDomain[counter_tmp]
			if(massive_type == double_type)
			{
				//If data massive is from type double
				MPI_Isend(&massive_send_to_SubDomain_tmp[0],
					N_send_to_counter_tmp, MPI_DOUBLE,
					counter_tmp,
					SubDomain[counter_tmp].CalculateTagForSendVariableToNeighbour(u_send, to_I_1J),
					MPI_COMM_WORLD,
					&request[CalculateIndexRequestForSendOrRecvVariableNeighbour(u_send, to_I_1J, index_send)]);
			}
			else if(massive_type == unsigned_int_type)
			{
				//If data massive is from type unsigned int
				MPI_Isend(&massive_send_to_SubDomain_tmp[0],
					N_send_to_counter_tmp, MPI_UNSIGNED,
					counter_tmp,
					SubDomain[counter_tmp].CalculateTagForSendVariableToNeighbour(u_send, to_I_1J),
					MPI_COMM_WORLD,
					&request[CalculateIndexRequestForSendOrRecvVariableNeighbour(u_send, to_I_1J, index_send)]);
			}


			//Wait to complete MPI_Isend
			MPI_Wait(&request[CalculateIndexRequestForSendOrRecvVariableNeighbour(u_send, to_I_1J, index_send)],
				&status[CalculateIndexRequestForSendOrRecvVariableNeighbour(u_send, to_I_1J, index_send)]);


			//This array MUST be deleted, because is defined in function.
			//The function not delete arrayes defined in in automaticaly.
			//Delete only local variables (this is make by compiler).
			delete [] massive_send_to_SubDomain_tmp;
		}



		//This array MUST be deleted, because is defined in function.
		//The function not delete arrayes defined in in automaticaly.
		//Delete only local variables (this is make by compiler).
		delete [] massive_entire_computational_domain_tmp;
	}
	else
	{
		//Irecv all array massive
		if(massive_type == double_type)
		{
			//If recv data are from type double
			MPI_Irecv(&massive[0], SubDomain[IJ].Na, MPI_DOUBLE,
				0, SubDomain[IJ].CalculateTagForSendVariableToNeighbour(u_send, to_I_1J),
				MPI_COMM_WORLD,
				&request[CalculateIndexRequestForSendOrRecvVariableNeighbour(u_send, to_I_1J, index_recv)]);
		}
		else if(massive_type == unsigned_int_type)
		{
			//If recv data are from type unsigned int
			MPI_Irecv(&massive[0], SubDomain[IJ].Na, MPI_UNSIGNED,
				0, SubDomain[IJ].CalculateTagForSendVariableToNeighbour(u_send, to_I_1J),
				MPI_COMM_WORLD,
				&request[CalculateIndexRequestForSendOrRecvVariableNeighbour(u_send, to_I_1J, index_recv)]);
		}


		//Wait to complete MPI_Irecv
		MPI_Wait(&request[CalculateIndexRequestForSendOrRecvVariableNeighbour(u_send, to_I_1J, index_recv)],
			&status[CalculateIndexRequestForSendOrRecvVariableNeighbour(u_send, to_I_1J, index_recv)]);

	}
}



void WriteDataFor_Nx_Ny_ForAllSubDomains(void)
{
	if(is_defined_arrays_for_Nx_Ny_for_all_SubDomains)
	{
		//Define temporal variables
		string FileName_tmp = "Nx_Ny_ForAllSubDomains.txt";

		using namespace std;

		ofstream output_data;
		output_data.open(FileName_tmp.c_str());

		if(output_data.is_open())
		{
			int counter_tmp;
			for (counter_tmp = 0; counter_tmp < N_SubDomains_a; counter_tmp++)
			{
				output_data
					<< Nx_all_SubDomains[counter_tmp] << "	"
					<< Ny_all_SubDomains[counter_tmp] << "	"
					<< "//Nx_all_SubDomains	Ny_all_SubDomains for subdomain IJ = " << counter_tmp
					<< endl;
			}

			output_data.close();
		}
		else
		{
			cout << "The program can not open file Nx_Ny_ForAllSubDomains.txt to write data." << endl;
		}


	}
	else
	{
		cout << "You try to write arrays Nx_all_SubDomains and Ny_all_SubDomains, \
				but they are not defined. Repair the source code!" << endl;
	}

}


void ReadDataFor_Nx_Ny_ForAllSubDomains(void)
{
	//Define temporal variables
	string FileName_tmp = "Nx_Ny_ForAllSubDomains.txt";

	if(is_defined_arrays_for_Nx_Ny_for_all_SubDomains)
	{
		Nx_all_SubDomains = NULL; delete [] Nx_all_SubDomains;
		Ny_all_SubDomains = NULL; delete [] Ny_all_SubDomains;
	}

	Nx_all_SubDomains = new int [N_SubDomains_a];
	Ny_all_SubDomains = new int [N_SubDomains_a];

	is_defined_arrays_for_Nx_Ny_for_all_SubDomains = true;

	using namespace std;

	const int MaxCharRead = 5000;
	char buffer[MaxCharRead];


	ifstream input_data;
	input_data.open(FileName_tmp.c_str());

	if(input_data.is_open())
	{
		int counter_tmp;
		for(counter_tmp = 0; counter_tmp < N_SubDomains_a; counter_tmp++)
		{
			input_data >> Nx_all_SubDomains[counter_tmp];
			input_data >> Ny_all_SubDomains[counter_tmp];
			input_data.getline(&buffer[0], MaxCharRead, '\n');
		}

		input_data.close();
	}
	else
	{
		cout << "The program can't open file " << FileName_tmp
			<< " to read data for Nx and Ny for SubDomains"
			<< " process rank = " << rank
			<< " Number_of_processes = " << Number_of_processes
			<< endl;
	}

}


void CollectDataForCDfromAllSubDomainsToProcessIJ0(void)
{
	//Make temporal variables
	double CD_fr_bottom_tmp, CD_p_bottom_tmp, S_for_CD_bottom_tmp,
		CD_fr_front_tmp, CD_p_front_tmp, S_for_CD_front_tmp,
		CD_fr_top_tmp, CD_p_top_tmp, S_for_CD_top_tmp,
		CD_fr_behind_tmp, CD_p_behind_tmp, S_for_CD_behind_tmp;

	//Initialize temporal variables
	CD_fr_bottom_tmp = 0;
	CD_p_bottom_tmp = 0;
	S_for_CD_bottom_tmp = 0;

	CD_fr_front_tmp = 0;
	CD_p_front_tmp = 0;
	S_for_CD_front_tmp = 0;

	CD_fr_top_tmp = 0;
	CD_p_top_tmp = 0;
	S_for_CD_top_tmp = 0;

	CD_fr_behind_tmp = 0;
	CD_p_behind_tmp = 0;
	S_for_CD_behind_tmp = 0;


	//Sum of CD from all subdomains
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&CD_fr_bottom, &CD_fr_bottom_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&CD_p_bottom, &CD_p_bottom_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&S_for_CD_bottom, &S_for_CD_bottom_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Reduce(&CD_fr_front, &CD_fr_front_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&CD_p_front, &CD_p_front_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&S_for_CD_front, &S_for_CD_front_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Reduce(&CD_fr_top, &CD_fr_top_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&CD_p_top, &CD_p_top_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&S_for_CD_top, &S_for_CD_top_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Reduce(&CD_fr_behind, &CD_fr_behind_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&CD_p_behind, &CD_p_behind_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&S_for_CD_behind, &S_for_CD_behind_tmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);


	if(IJ == 0)
	{
		//Copy to the variable
		CD_fr_bottom = CD_fr_bottom_tmp;
		CD_p_bottom = CD_p_bottom_tmp;
		S_for_CD_bottom = S_for_CD_bottom_tmp;

		CD_fr_front = CD_fr_front_tmp;
		CD_p_front = CD_p_front_tmp;
		S_for_CD_front = S_for_CD_front_tmp;

		CD_fr_top = CD_fr_top_tmp;
		CD_p_top = CD_p_top_tmp;
		S_for_CD_top = S_for_CD_top_tmp;

		CD_fr_behind = CD_fr_behind_tmp;
		CD_p_behind = CD_p_behind_tmp;
		S_for_CD_behind = S_for_CD_behind_tmp;
	}

	MPI_Barrier(MPI_COMM_WORLD);

}



#endif //End of COMPILE_MPI




//This is begin of file for class SIMPLEST_ANL_compr_T_h_var_NoDim_Algorithm

