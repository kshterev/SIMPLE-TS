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

//In this file DomainDecomposition.h
//is defined struct to define data for all SubDomains.
//This struct conain data which SubDomain is current and
//who are the neighbourhood SubDomains, which is importantwhen
//when data are exchanged.


#ifndef DomainDecomposition2D_H
#define DomainDecomposition2D_H

//#pragma once


struct DomainDecomposition2D
{
	DomainDecomposition2D(void);
	~DomainDecomposition2D(void);

	//Variables for DomainDecomposition
	//Number of SubDomains on OX
	int N_SubDomains_x;
	//Number of SubDomains on OY
	int N_SubDomains_y;
	
	//N_SubDomains_a - number of all SubDomains
	int N_SubDomains_a;

	//Number of nodes on OX in this SubDomain
	int Nx;
	//Number of nodes on OY in this SubDomain
	int Ny;
	//Number of all nodes in this SubDomain
	int Na;

	//This SubDomain is IJ
	//The left neighbourhood is I_1J
	//The right neighbourhood is I1J
	//The bottom neighbourhood is IJ_1
	//The top neighbourhood is IJ1
	int I, J, IJ, I_1J, I1J, IJ_1, IJ1;

	//Information about neighbourhood SubDomans
	bool is_SubDomain_in_direction_I_1, is_SubDomain_in_direction_I1,
		is_SubDomain_in_direction_J_1, is_SubDomain_in_direction_J1;


	//Begin and end index of SubDomain on OX, according to entire computational domain,
	//note that begin may be 0, end may be Nx - 1.
	int i_b, i_e;
	//Begin and end index of SubDomain on OY, according to entire computational domain,
	//note that begin may be 0, end may be Ny - 1.
	int j_b, j_e;

	//maximum numbet of tag for this kind of problems
	//are N_tag_max = N_SubDomains_a * N_neighbouhoods * N_exchanged_variables;
	//where:
	//N_Domains_a - number of all domains
	//N_neighbouhoods - number of neighbouhoods
	//		- for 1D case -> N_neighbouhoods = 2
	//		- for 2D case -> N_neighbouhoods = 4
	//		- for 3D case -> N_neighbouhoods = 6
	//
	//N_exchanged_variables_with_neighbouhoods - number of axchanged variabled with neighbouhoods
	//		- in this case N_exchanged_variables = 10
	//			(d_p_c_12, apx, bpx, d_p_c_34, apy, bpy, u, v, p, Temper)
	//
	//N_send_variables_to_IJ0 - number of send variables to process IJ == 0
	//		- in this case N_send_variables_to_IJ0 = 4
	//			(max_residual_in_p, max_residual_in_T, max_residual_in_u, max_residual_in_v)
	//
	//N_send_from_IJ0 - number of send variables from process IJ == 0 to all SubDomains.
	//		- in this case N_send_from_IJ0 = 2
	//			(continue_iter_int, countinue_time_step_int)
	//
	//The equation to calculate N_tag_max is:
	//N_tag_max = (N_exchanged_variables_with_neighbouhoods * N_neighbouhoods
	//				+ N_send_variables_to_IJ0) * N_Domains_a
	//			+ N_send_from_IJ0;
	//
	//For simplification the tags will be calculated by equation
	//N_tag_max = (N_exchanged_variables_with_neighbouhoods * N_neighbouhoods
	//				+ N_send_variables_to_IJ0 + N_send_from_IJ0) * N_Domains_a;
	//
	//Where the tag with indexes continue_iter_int_send_from_IJ0 and
	//countinue_time_step_int_send_from_IJ0 always will be number of tag
	//of process IJ == 0, respectively 44 and 45.
	//
	//For this case we have:
	//N_exchanged_variables_with_neighbouhoods = 10;
	//N_neighbouhoods = 4;
	//N_send_variables_to_IJ0 = 4;
	//N_send_from_IJ0 = 2;
	//N_tag_max = (10 * 4 + 4 + 2) * N_Domains_a;
	int N_tag_max, N_neighbouhoods, N_exchanged_variables_with_neighbouhoods,
		N_send_variables_to_IJ0, N_send_from_IJ0;


	//zero_for_this_tag_for_this_IJ - zero tag for this process.
	//It is used the equation to calculate this variable:
	//zero_for_this_tag_for_this_IJ = (N_exchanged_variables_with_neighbouhoods * N_neighbouhoods
	//										+ N_send_variables_to_IJ0 + N_send_from_IJ0) * IJ;
	int zero_tag_for_this_for_this_IJ;



	/*
	Target:		Calculate tag for send variable to neighbour. This function return
				value for tag, which is calculated using equation in depend of variable
				neighbour and process from, where will be send data. The equation just
				arrange the tags. This tag is unique.
	Receive:	variable_tmp - variable which will be send. The variable for which can
				be used this function are: d_p_c_12_send, apx_send, bpx_send,
				d_p_c_34_send, apy_send, bpy_send, p_send, Temper_send, u_send and
				v_send.
				
				neighbour_tmp - neighbour to which will be send variable. The variables
				are to_I_1J, to_I1J, to_IJ_1 and to_IJ1.
	*/
	inline int CalculateTagForSendVariableToNeighbour(
		const int& variable_tmp, const int& neighbour_tmp)
	{
		return (zero_tag_for_this_for_this_IJ
			+ variable_tmp * N_neighbouhoods + neighbour_tmp);
	};

	/*
	Target:		Calculate tag for send variable to process IJ == 0. This function return
				value for tag, which is calculated using equation in depend of variable
				and process from, where will be send data. The equation just arrange
				the tags. This tag is unique.
	Receive:	variable_tmp - variable which will be send. The variable for which can
				be used this function are: max_residual_in_p_send_to_IJ0,
				max_residual_in_T_send_to_IJ0, max_residual_in_u_send_to_IJ0 and
				max_residual_in_v_send_to_IJ0.
	*/
	inline int CalculateTagForSendVariableToIJ0(
		const int& variable_tmp)
	{
		return (zero_tag_for_this_for_this_IJ
			+ N_exchanged_variables_with_neighbouhoods * N_neighbouhoods
			+ variable_tmp);
	};

	/*
	Target:		Calculate tag for send variable from process IJ == 0 to all processes
				(subdomains). This function return value for tag, which is calculated
				using equation in depend of variable and process, where will be send
				data. The equation just arrange the tags. This tag is unique.
	Receive:	variable_tmp - variable which will be send. The variable for which can
				be used this function are: continue_iter_int_send_from_IJ0 and
				countinue_time_step_int_send_from_IJ0.
	*/
	inline int CalculateTagForSendVariableFromIJ0(
		const int& variable_tmp)
	{
		return (zero_tag_for_this_for_this_IJ
			+ N_exchanged_variables_with_neighbouhoods * N_neighbouhoods
			+ N_send_variables_to_IJ0
			+ variable_tmp);
	};


	//Periodic_boundary_conditions_about_OX for all computational domain
	bool Periodic_boundary_conditions_about_OX;
	
	//Define Indexes
	//if I_1J == IJ that mean that there is no SubDomain on left (is_SubDomain_in_direction_I_1 == false)
	//if I1J == IJ that mean that there is no SubDomain on right (is_SubDomain_in_direction_I1 == false) 
	//if IJ_1 == IJ that mean that there is no SubDomain on bottom (is_SubDomain_in_direction_J_1 == false)
	//if IJ1 == IJ that mean that there is no SubDomain on top (is_SubDomain_in_direction_J1 == false)
	//Also calculate tags.
	//The I, J and Periodic_boundary_conditions_about_OX must be defined before use thie function.
	void DefineIndexesAndTags(
		const int& I_tmp, const int& J_tmp, const bool& Periodic_boundary_conditions_about_OX_tmp);

};


#endif //#ifndef DomainDecomposition2D_H

