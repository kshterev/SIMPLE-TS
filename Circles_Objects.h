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


//File name: Circles_Objects.h
#ifndef Circles_Objects_H
#define Circles_Objects_H

//#pragma once

#define _USE_MATH_DEFINES
#include <math.h> //This is for mathematical calculation


struct /*class*/ CCircles_Objects
{
//public:
	CCircles_Objects(void);
	~CCircles_Objects(void);


	//Variables for circle
	//Center of the circle by OX
	double x_c;
	//Center of the circle by OY
	double y_c;
	//Radius of the circle
	double r_c;

	//Horizontal velocity of the body circle
	double u_body;
	//Vertical velocity of the body circle
	double v_body;

	//Pressure of the body circle
    double p_body;

	//Temperature of the body circle
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
	//where: sigma_V represent energy accomodation coefficient
	//sigma_T = [0.0, 1.0]
	//Note: F_TemperatureJump, can be 0.292053333333334 or something else
	double F_TemperatureJump;

	//dTdn_on_wall_0 == true - normal derivative to the surface of the body
	//is zero, i.e. the temperature on the wall is equal to the temperature
	//of the gas to the wall.
	//dTdn_on_wall_0 == false - the temperature on the body is fixed.
	double dTdn_on_wall_0;


	//ToSolveDragCoefficientForThisBody - if that is true the drag coefficient will be solved
	//and data will be write to file with name CD_ + number of body.
	bool ToSolveDragCoefficientForThisBody;


	//Functions for circle
    template <typename T_x, typename T_y>
	bool is_point_inside_circle(const T_x& x, const T_y& y)
	/*	Target:		to return true if point(x, y) is in circle
		Recieve:	coordinates of point which we check - x and y
		Return:		true if point is inside and false if point is outside circle
	*/
	{
		 return (((pow((x - x_c), 2) + pow((y - y_c), 2))
			 <= pow(r_c, 2)));
	};


};


#endif //#ifndef Circles_Objects_H
