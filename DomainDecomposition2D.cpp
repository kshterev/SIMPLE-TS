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

#include "DomainDecomposition2D.h"


////contain informatin about all SubDomains.
//extern DomainDecomposition2D * SubDomain;


DomainDecomposition2D::DomainDecomposition2D(void)
{
	N_SubDomains_x = 0;
	N_SubDomains_y = 0;
	N_SubDomains_a = 0;

	Nx = 0;
	Ny = 0;
	Na = 0;

	I = 0;
	J = 0;
	IJ = 0;
	I_1J = 0;
	I1J = 0;
	IJ_1 = 0;
	IJ1 = 0;

	is_SubDomain_in_direction_I_1 = false;
	is_SubDomain_in_direction_I1 = false;
	is_SubDomain_in_direction_J_1 = false;
	is_SubDomain_in_direction_J1 = false;

	i_b = 0;
	i_e = 0;
	j_b = 0;
	j_e = 0;

	N_tag_max = 0;
	N_neighbouhoods = 0;
	N_exchanged_variables_with_neighbouhoods = 0;
	N_send_variables_to_IJ0 = 0;
	N_send_from_IJ0 = 0;

	zero_tag_for_this_for_this_IJ = 0;

	Periodic_boundary_conditions_about_OX = false;
}

DomainDecomposition2D::~DomainDecomposition2D(void)
{
}

void DomainDecomposition2D::DefineIndexesAndTags(
		const int& I_tmp, const int& J_tmp, const bool& Periodic_boundary_conditions_about_OX_tmp)
{
	I = I_tmp;
	J = J_tmp;
	Periodic_boundary_conditions_about_OX = Periodic_boundary_conditions_about_OX_tmp;

	//Calculate indexes
	IJ = I + J * N_SubDomains_x;

	if(Periodic_boundary_conditions_about_OX)
	{
		if(I == 0)
		{
			//The I_1 neighbouhood is last subdomain, if(I == 0)
			I_1J = (N_SubDomains_x - 1) + J * N_SubDomains_x;
		}
		else
		{
			I_1J = IJ - 1 * (0 < I);
		}


		if(I == (N_SubDomains_x - 1))
		{
			//The I1 neighbouhood is first (zero) subdomain, if(I == (N_SubDomains_x - 1))
			I1J = 0 + J * N_SubDomains_x;
		}
		else
		{
			I1J = IJ + 1 * (I < (N_SubDomains_x - 1));
		}


		//If there are more then 1 subdomains on OX, there is neighbouhood
		is_SubDomain_in_direction_I_1 = (1 < N_SubDomains_x);

		//If there are more then 1 subdomains on OX, there is neighbouhood
		is_SubDomain_in_direction_I1 = (1 < N_SubDomains_x);
	}
	else
	{
		I_1J = IJ - 1 * (0 < I);
		I1J = IJ + 1 * (I < (N_SubDomains_x - 1));

		is_SubDomain_in_direction_I_1 = (0 < I);
		is_SubDomain_in_direction_I1 = (I < (N_SubDomains_x - 1));
	}

	IJ_1 = IJ - N_SubDomains_x * (0 < J);
	IJ1 = IJ + N_SubDomains_x * (J < (N_SubDomains_y - 1));

	is_SubDomain_in_direction_J_1 = (0 < J);
	is_SubDomain_in_direction_J1 = (J < (N_SubDomains_y - 1));



	//Calculate maximum naumber of tags:
	N_tag_max = (N_exchanged_variables_with_neighbouhoods * N_neighbouhoods
					+ N_send_variables_to_IJ0 + N_send_from_IJ0) * N_SubDomains_a;

	zero_tag_for_this_for_this_IJ = (N_exchanged_variables_with_neighbouhoods * N_neighbouhoods
									+ N_send_variables_to_IJ0 + N_send_from_IJ0) * IJ;

}



