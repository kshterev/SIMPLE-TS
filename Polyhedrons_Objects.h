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


//File name: Polyhedons_Objects.h
#ifndef Polyhedons_Objects_H
#define Polyhedons_Objects_H

struct /*class*/ Polyhedrons_Objects
{
//public:
	Polyhedrons_Objects(void);
	~Polyhedrons_Objects(void);


	//Matrix which contain data for vertex of polyhedron
	unsigned int N_vertex_body;

	//Vector of x coordinate of vertexes
	double * x_coordinate;
	//Vector of y coordinate of vertexes
	double * y_coordinate;


	//Coordinate for first point inside polyhedron
	double point1_inside_x;
	double point1_inside_y;
	//Coordinate for second point inside polyhedron
	double point2_inside_x;
	double point2_inside_y;
	//Coefficients for line between points inside (point1 end point2): y = kpi * x - bpi
	double kpi;
	double bpi;
	//This is min and max of point inside polyhedron - this is nesesarry when we deside where is intersection point
	double x_p_min, x_p_max, y_p_min, y_p_max;


	//Horizontal velocity of the body polihedron
	double u_body;
	//Vertical velocity of the body polihedron
	double v_body;

	//Pressure of the body polihedron
    double p_body;

	//Temperature of the body polihedron
    double T_body;


	//Variable for Velocity Slip Boundary conditions
	//VelocitySlipBC_body == false; ==> that mean that body have no slip boundary conditions
	//VelocitySlipBC_body == true; ==> that mean that body have velocity slip boundary conditions
	//with coefficent F_VelocitySlip.
	//Condition is: v_gas_on_the_wall = F_VelocitySlip * Kn_local * dvdn_wall
	//where:
	//Kn_local - local knudsen number, depend of way of non-duimensialization
	//dvdn_wall - this is the derivation of velocity again normal vektor to the wall
	bool is_VelocitySlipBC_body;

	//F_VelocitySlip - slip coefficients
	//F_VelocitySlip = (2.0 - sigma_V) / sigma_V,
	//where: sigma_V represent streamwise momentum accommodation coefficient
	//sigma_V = [0.0, 1.0]
	//Note: F_VelocitySlip, can be 1.1466 or something else
	double F_VelocitySlip;

	//w_VelocitySlipBC - under-relacation coefficient for Velocity Slip BC.
	//Depend from Kn_local and from step of mesh on space.
	double w_VelocitySlipBC;

	//Variable for Temperature Jump Boundary conditions
	//TemperatureJumpBC_body == false; ==> that mean that body have no temperature jump boundary conditions
	//TemperatureJumpBC_body == true; ==> that mean that body have temperature jump boundary conditions
	//with coefficent sigma.
	//Condition is: T_gas_on_the_wall = F_VelocitySlip * Kn_local * dvdn_wall
	//where:
	//Kn_local - local knudsen number, depend of way of non-duimensialization
	//dvdn_wall - this is the derivation of velocity again normal vektor to the wall
	bool is_TemperatureJumpBC_body;

	//F_TemperatureJump - Temperature Jump coefficients
	//F_TemperatureJump = (2.0 - sigma_T) / sigma_T,
	//where: sigma_V represent energy accommodation coefficient
	//sigma_T = [0.0, 1.0]
	//Note: F_TemperatureJump, can be 0.29205333333333 or something else
	double F_TemperatureJump;

	//dTdn_on_wall_0 == true - normal derivative to the surface of the body
	//is zero, i.e. the temperature on the wall is equal to the temperature
	//of the gas to the wall.
	//dTdn_on_wall_0 == false - the temperature on the body is fixed.
	double dTdn_on_wall_0;

	//ToSolveDragCoefficientForThisBody - if that is true the drag coefficient will be solved
	//and data will be write to file with name CD_ + number of body.
	bool ToSolveDragCoefficientForThisBody;


	//Variable for rectangular around polyhedron
	double x_rect_min, x_rect_max, y_rect_min, y_rect_max;

	//variable to define is already solved ractangular around polyhedron
	bool rectangular_around_polyhedron_is_defined;

	//Define rectangular around polyhedron
	void define_rectangular_around_polyhedron(void);


	//Functions for polyhedron
	void define_vector_for_vertexes(void);

	bool is_data_for_polyhedron_correct(void);

	bool is_point_inside_polyhedron(const double& x, const double& y);


	bool is_two_lines_intersections(const double& x11, const double& y11, const double& x12, const double& y12,
								const double& x21, const double& y21, const double& x22, const double& y22);


	template <typename T_a_max, typename T_b_max>
	inline double maximum(const T_a_max& a, const T_b_max& b)
	{
		return (a > b ? a : b);
	};

	template <typename T_a_min, typename T_b_min>
	inline double minimum(const T_a_min& a, const T_b_min& b)
	{
		return (a < b ? a : b);
	};
};


#endif //ifndef Polyhedons_Objects_H
