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


#include "Polyhedrons_Objects.h"

#include "iostream"
#include "math.h" //This is for mathematical calculation

using namespace std;


const unsigned int N_maximum_divide = 5000; //This is maximum of dividing line between points inside polyhedron
const double mistake = 1e-5; //This is mistake from solving radius of circle. At that way the radius is a little bit more then ecxactly, but thats not important

Polyhedrons_Objects::Polyhedrons_Objects(void)
{
	N_vertex_body = 0;

	u_body = 0;
	v_body = 0;

	p_body = 0;
	T_body = 0;

	is_VelocitySlipBC_body = false;
	F_VelocitySlip = 1.1466;
	w_VelocitySlipBC = 0.1;

	is_TemperatureJumpBC_body = false;
	F_TemperatureJump = 2.1904;
	dTdn_on_wall_0 = 0;


	ToSolveDragCoefficientForThisBody = false;

	point1_inside_x = 0;
	point1_inside_y = 0;

	point2_inside_x = 1;
	point2_inside_y = 1;

	rectangular_around_polyhedron_is_defined = false;
}

Polyhedrons_Objects::~Polyhedrons_Objects(void)
{
}


void Polyhedrons_Objects::define_vector_for_vertexes(void)
//Make vector of vertexes
{
	x_coordinate = new double [N_vertex_body];
	y_coordinate = new double [N_vertex_body];
}


bool Polyhedrons_Objects::is_data_for_polyhedron_correct(void)
/*	Target:	to check already entered data for polyhedron
	Return:	true if all data are correct
			false if data are not correct


	Step for cheking validation of entered data:
	1. Check for coincidence of vertexes of polyhedron
	2. Check for coincidence of points iside
	3. Check for coincident for lines between points inside and one of vertexes
	4. Check if line inside intersect sides of polyhedron

*/
{
	bool continue_checking = true;


	unsigned int counter = 0;
	//1. Check for coincidence of vertexes of polyhedron
	do
	{
		if(counter == N_vertex_body - 1)
		{
			if(y_coordinate[counter] == y_coordinate[0]
				&& x_coordinate[counter] == x_coordinate[0])
			{
				cout << "Two of vertexes coincident and the check will be stopped! You must enter the information again." << endl;

				continue_checking = false;
			}
		}
		else
		{
			if(y_coordinate[counter] == y_coordinate[counter + 1] && x_coordinate[counter] == x_coordinate[counter + 1])
			{
					cout << "Two of vertexes coincident and the check will be stopped! You must enter the information again." << endl;

					continue_checking = false;
			}
		}

		counter++;
	}while(counter < N_vertex_body && continue_checking);


	//2. Check for coincidence of points iside
	if(point1_inside_x == point2_inside_x
		&& point1_inside_y == point2_inside_y)
	{
		cout << "The point inside are coincident! You must enter the information again." << endl;

		continue_checking = false;
	}


	//3. Check for coincident for lines between points inside and one of vertexes
	if(continue_checking)
	{
		counter = 0;
		do
		{
			if(y_coordinate[counter] == (((point2_inside_y - point1_inside_y) * (x_coordinate[counter] - point1_inside_x))
				/ (point2_inside_x - point1_inside_x) - point1_inside_y))
			{
				cout << "The point inside are at the same straight line with one of vertexes and the check will be stopped! You must enter the information again." << endl;

				continue_checking = false;
			}

			counter++;
		}while(counter < N_vertex_body && continue_checking);
	}


	//4. Check if line inside intersect sides of polyhedron
	if(continue_checking)
	{
		counter = 0;
		do
		{
			if(counter < N_vertex_body - 1)
			{
				continue_checking = !is_two_lines_intersections(point1_inside_x, point1_inside_y, point2_inside_x, point2_inside_y,
					x_coordinate[counter], y_coordinate[counter], x_coordinate[counter + 1], y_coordinate[counter + 1]);
			}
			else
			{
				continue_checking = !is_two_lines_intersections(point1_inside_x, point1_inside_y, point2_inside_x, point2_inside_y,
					x_coordinate[counter], y_coordinate[counter], x_coordinate[0], y_coordinate[0]);
			}


			counter++;

		}while(counter < N_vertex_body && continue_checking);


		if(!continue_checking)
			cout << "The line between points inside intersect one of sides of polyhedron! You must enter the information again." << endl;

	}



	return(continue_checking);
}


bool Polyhedrons_Objects::is_point_inside_polyhedron(const double& x, const double& y)
/*	Target:		to return true if point(x, y) is in polyhedron
	Recieve:	coordinates of point which we check - x and y
	Return:		true if point is inside and false if point is outside polyhedron
*/
{
	bool continue_checking = true;

	int counter, dividing_points_inside;

	//Number of intersection between line for check ond side of polyhedron
	//if N_intersections = an even number the point is outside polihedron
	//else point is inside
	double N_intersections;

	//Point inside - this is point whith which we will make checking
	double point_inside_x, point_inside_y;


	//Solving equation for rectangular arround polyhedron
	if(!rectangular_around_polyhedron_is_defined)
		define_rectangular_around_polyhedron();

	//if point is outside rectangular the checking will be stoped
	continue_checking = (x_rect_min <= x && x <= x_rect_max
						&& y_rect_min <= y && y <= y_rect_max);


	//if point(x,y) is not in circle then the point is not in the polihedron
	if(continue_checking)
	{
		//Here is solving min and max coordinate by OX axis for points inside
		x_p_min = minimum(point1_inside_x, point2_inside_x);
		x_p_max = maximum(point1_inside_x, point2_inside_x);

		//Here is solving min and max coordinate by OY axis for points inside
		y_p_min = minimum(point1_inside_y, point2_inside_y);
		y_p_max = maximum(point1_inside_y, point2_inside_y);


		//Solve coefficient for line inside polyhdron
		if((x_p_min != x_p_max) && (y_p_min != y_p_max))
		{
			kpi = (point2_inside_y - point1_inside_y) / (point2_inside_x - point1_inside_x);
			bpi = point1_inside_y - kpi * point1_inside_x;
		}



		//To define inside point of polyhedron for check
		bool continiue_dividing;
		dividing_points_inside = 1;//This is becouse inside point1 may be not good for solving and must be another
		do
		{
			continiue_dividing = false;

			//Solve coordinate for inside point of polyhedron for chech
			if(x_p_min == x_p_max)
			{
				point_inside_x = x_p_min;
				point_inside_y = y_p_min + (y_p_max - y_p_min) / dividing_points_inside;
			}
			else if(y_p_min == y_p_max)
			{
				point_inside_x = x_p_min + (x_p_max - x_p_min) / dividing_points_inside;
				point_inside_y = y_p_min;
			}
			else
			{
				point_inside_x = x_p_min + (x_p_max - x_p_min) / dividing_points_inside;
				point_inside_y = kpi * point_inside_x + bpi;
			}


			//Check if point inside and point for check are coincident
			if(point_inside_x == x && point_inside_y == y)
				continiue_dividing = true;


			//Check if point inside, point what we check and any of vertexes are at the same line
			if(!continiue_dividing)
			{
				counter = 0;
				do
				{
					if(point_inside_x == x && x == x_coordinate[counter])
						continiue_dividing = true;
					else if(point_inside_y == y && y == y_coordinate[counter])
						continiue_dividing = true;
					else if(y_coordinate[counter] == (point_inside_y - y) * (x_coordinate[counter] - point_inside_x) / (point_inside_x - x) + point_inside_y)
						continiue_dividing = true;


					counter++;
				}while(counter < N_vertex_body && !continiue_dividing);

			}



			dividing_points_inside++;
		}while(dividing_points_inside < N_maximum_divide && continiue_dividing);


		//If can not be found point inside polihedron to make check the program will show window for mistake
		if(continiue_dividing)
			cout << "The check with point inside can NOT be done! The divide is made maximum division. You must stop solving and give another points inside polyhedron." << endl;
		else
		{
			//Cheking how many side of polyhedron are intersecsion by line for check
			N_intersections = 0;

			for(counter = 0; counter < N_vertex_body; counter++)
			{
				if(counter < N_vertex_body - 1)
					N_intersections = N_intersections + is_two_lines_intersections(point_inside_x, point_inside_y, x, y,
						x_coordinate[counter], y_coordinate[counter], x_coordinate[counter + 1], y_coordinate[counter + 1]);
				else
					N_intersections = N_intersections + is_two_lines_intersections(point_inside_x, point_inside_y, x, y,
						x_coordinate[counter], y_coordinate[counter], x_coordinate[0], y_coordinate[0]);

			}

			//Deside if the checking point is in polyhedron
			continue_checking = ((N_intersections / 2.0) == floor(N_intersections / 2.0));

		}


	}



	return(continue_checking);
}


bool Polyhedrons_Objects::is_two_lines_intersections(const double& x11, const double& y11, const double& x12, const double& y12,
													const double& x21, const double& y21, const double& x22, const double& y22)
/*	Target:		to return true if two lines wre intersections
	Recieve:	x11, x12 - coordinates point 1 and 2 on OX of first line
				y11, y12 - coordinates point 1 and 2 on OY of first line
				x21, x22 - coordinates point 1 and 2 on OX of second line
				y21, y22 - coordinates point 1 and 2 on OY of second line


	Step for cheking intersection of lines:
	1. Define variables:
		//Coordinate for intersection point
		double x_intersection, y_intersection;
		//This is maximum ond minimum of coordinate of points of fist line
		double x1min, x1max, y1min, y1max;
		//This is maximum ond minimum of coordinate of points of second line
		double x2min, x2max, y2min, y2max;

	2. Solving min and max coordinate by OX axis for first line
	3. solving min and max coordinate by OY axis for first line
	4. Solving min and max coordinate by OX axis for second line
	5. Solving min and max coordinate by OY axis for second line

	6. There is nine combinations between two lines. For each line, line can be horisontal, vertical and slanting
	I'll check all, ecxept when lines are horizontal or vertical, bu they are parallel

	7. Check for coordinate of intersection point
		The check is in rectangular:	x_min <= x <= x_max
										y_min <= y <= y_max
		for two point iside polyhedron and two vertex of side.


*/
{
	bool lines_intersect = false;


	//1. Define variables
	//Coordinate for intersection point
	double x_intersection, y_intersection;

	//This is maximum ond minimum of coordinate of points of fist line
	double x1min, x1max, y1min, y1max;

	//This is maximum ond minimum of coordinate of points of second line
	double x2min, x2max, y2min, y2max;


	//2. Solving min and max coordinate by OX axis for first line
	x1min = minimum(x11, x12);
	x1max = maximum(x11, x12);

	//3. solving min and max coordinate by OY axis for first line
	y1min = minimum(y11, y12);
	y1max = maximum(y11, y12);


	//4. Solving min and max coordinate by OX axis for second line
	x2min = minimum(x21, x22);
	x2max = maximum(x21, x22);

	//5. Solving min and max coordinate by OY axis for second line
	y2min = minimum(y21, y22);
	y2max = maximum(y21, y22);


	//Coefficient for first and second line
	double k1, b1, k2, b2;


	//6. There is nine combinations between two lines. For each line, line can be horisontal, vertical and slanting
	//I'll check all, ecxept when lines are horizontal or vertical, bu they are parallel
	if((x1min != x1max) && (y1min != y1max) && (x2min != x2max) && (y2min != y2max))
	{
		//first is slanting, second is slanting

		//Solving coefficient of first line
		k1 = (y11 - y12) / (x11 - x12);
		b1 = y11 - k1 * x11;

		//Solving coefficient of second line
		k2 = (y21 - y22) / (x21 - x22);
		b2 = y21 - k2 * x21;

		if(k1 != k2)
		{
			x_intersection = (b1 - b2) / (k2 - k1);
			y_intersection = k1 * x_intersection + b1;

			lines_intersect = true;
		}

	}
	else if(y1min == y1max && x2min == x2max)
	{
		//first is horisontal, second is vertical
		x_intersection = x2min;
		y_intersection = y1min;

		lines_intersect = true;
	}
	else if((y1min == y1max) && (x2min != x2max) && (y2min != y2max))
	{
		//first is horizontal, second is slanting
		y_intersection = y1min;

		k2 = (y21 - y22) / (x21 - x22);
		b2 = y21 - k2 * x21;
		x_intersection = (y_intersection - b2) / k2;

		lines_intersect = true;
	}
	else if(x1min == x1max && y2min == y2max)
	{
		//first is verical, second is horisontal
		x_intersection = x1min;
		y_intersection = y2min;

		lines_intersect = true;
	}
	else if((x1min == x1max) && (x2min != x2max) && (y2min != y2max))
	{
		//first is vertical, second is slanting
		x_intersection = x1min;

		k2 = (y21 - y22) / (x21 - x22);
		b2 = y21 - k2 * x21;

		y_intersection = k2 * x_intersection + b2;

		lines_intersect = true;
	}
	else if((x1min != x1max) && (y1min != y1max) && (x2min == x2max))
	{
		//first is slanting, second is vertical
		x_intersection = x2min;

		k1 = (y11 - y12) / (x11 - x12);
		b1 = y11 - k1 * x11;
		y_intersection = k1 * x_intersection + b1;

		lines_intersect = true;
	}
	else if((x1min != x1max) && (y1min != y1max) && (y2min == y2max))
	{
		//first is slanting, second is horisontal
		y_intersection = y2min;

		k1 = (y11 - y12) / (x11 - x12);
		b1 = y11 - k1 * x11;
		x_intersection = (y_intersection - b1) / k1;

		lines_intersect = true;
	}
	else
		lines_intersect = false;


	//7. Check for coordinate of intersection point
	/*The check is in rectangular:	x_min <= x <= x_max
									y_min <= y <= y_max
	for two point iside polyhedron and two vertex of side.
	*/
	if(lines_intersect &&
		!(x1min <= x_intersection && x_intersection <= x1max
		&& x2min <= x_intersection && x_intersection <= x2max
		&& y1min <= y_intersection && y_intersection <= y1max
		&& y2min <= y_intersection && y_intersection <= y2max))
	{
		lines_intersect = false;
	}


	return(lines_intersect);
}




void Polyhedrons_Objects::define_rectangular_around_polyhedron(void)
/*	Target:	to define rectangular around polyhedron
*/
{
	//Given first approximation
	x_rect_min = x_coordinate[0];
	x_rect_max = x_coordinate[0];
	y_rect_min = y_coordinate[0];
	y_rect_max = y_coordinate[0];

	unsigned int counter;
	for(counter = 0; counter < N_vertex_body; counter++)
	{
		x_rect_min = minimum(x_rect_min, x_coordinate[counter]);
		x_rect_max = maximum(x_rect_max, x_coordinate[counter]);

		y_rect_min = minimum(y_rect_min, y_coordinate[counter]);
		y_rect_max = maximum(y_rect_max, y_coordinate[counter]);
	}

	rectangular_around_polyhedron_is_defined = true;

}
