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


//File name: Objects_for_Solving.h
#ifndef Objects_for_Solving_H
#define Objects_for_Solving_H


//#pragma once

#include "Circles_Objects.h"
#include "Polyhedrons_Objects.h"

//#include <cstring>
//#include "atlstr.h"

#include <iostream>
#include <fstream> //This is for writing and reading data from file

using namespace std;


struct /*class*/ Objects_for_Solving
{
//public:
	Objects_for_Solving(void);
	~Objects_for_Solving(void);


	//Numbers of circles for solving
	unsigned int N_circles;
	//Matrix which contain data for circles
	CCircles_Objects * M_circles;
	//Variable which shows are data in matrix for circle
	//are entered - false means that data are not inside, true that data are inside
	bool is_enter_M_circles;

	//This function make new matrix M_circles_b and delete previous if there is
	void make_new_M_circles(void)
	{
		if(is_enter_M_circles)
		{
			delete[] M_circles;
		}

		M_circles = new CCircles_Objects [N_circles];

		is_enter_M_circles = true;
	};


	//Numbers of polyhedrons for solving
	unsigned int N_polyhedrons;
	//Matrix which contain data for polyhedrons
	Polyhedrons_Objects * M_polyhedrons;


	//Variable which shows are data in matrix vor polyhedrons
	//are entered - false means that data are not inside true that data are inside
	bool is_enter_M_polyhedrons;


	//This function make new matrix M_polyhedrons and delete previous if there is
	void make_new_M_polyhedrons(void)
	{
		if(is_enter_M_polyhedrons)
		{
			delete[] M_polyhedrons;
		}

		M_polyhedrons = new Polyhedrons_Objects [N_polyhedrons];

		is_enter_M_polyhedrons = true;
	};


	template <typename T_FileNameEC>
	void Export_data_for_circles_to_file(const T_FileNameEC& FileName)
	/*	Target: to export data for circles to file
				with name FileName
	*/
	{
		if(!is_enter_M_circles)
			cout << "There is no circles to export in file: "
				<< FileName << endl;
		else
		{
			ofstream output_data;
			output_data.open(FileName);

			output_data << N_circles << "	//N_circles" << endl;

			unsigned int counter;

			for(counter = 0; counter < N_circles; counter++)
			{
				output_data
					<< M_circles[counter].x_c << "	"
					<< M_circles[counter].y_c << "	"
					<< M_circles[counter].r_c << "	//x_c	y_c	r_c" << endl

					<< M_circles[counter].u_body << "	"
					<< M_circles[counter].v_body << "	//u_body	v_body" << endl

					<< M_circles[counter].p_body << "	//p_body" << endl
					<< M_circles[counter].T_body << "	//T_body" << endl

					<< M_circles[counter].is_VelocitySlipBC_body << "	//is_VelocitySlipBC_body" << endl
					<< M_circles[counter].F_VelocitySlip << "	//F_VelocitySlip" << endl
					<< M_circles[counter].w_VelocitySlipBC << "	//w_VelocitySlipBC" << endl

					<< M_circles[counter].is_TemperatureJumpBC_body << "	//is_TemperatureJumpBC_body" << endl
					<< M_circles[counter].F_TemperatureJump << "	//F_TemperatureJump" << endl
					<< M_circles[counter].dTdn_on_wall_0 << "	//dTdn_on_wall_0" << endl

					<< M_circles[counter].ToSolveDragCoefficientForThisBody << "	//ToSolveDragCoefficientForThisBody" << endl;
			}

			output_data.close();


			cout << "Data for circles was succesfully exported to file: "
				<< FileName << endl;

		}

	}


	template <typename T_FileNameIC>
	void Import_data_for_circles_from_file(const T_FileNameIC& FileName)
	/*	Target: to import data dor circles from file
				with name file_name_for_circles
	*/
	{
		using namespace std;

		const int MaxCharRead = 5000;
		char buffer[MaxCharRead];


		ifstream input_data;
		input_data.open(FileName);

		if(input_data.is_open())
		{
			input_data >> N_circles;
			input_data.getline(&buffer[0], MaxCharRead, '\n');

			make_new_M_circles();

			unsigned int counter;

			for(counter = 0; counter < N_circles; counter++)
			{
				input_data >> M_circles[counter].x_c;
				input_data >> M_circles[counter].y_c;
				input_data >> M_circles[counter].r_c;
				input_data.getline(&buffer[0], MaxCharRead, '\n');

				input_data >> M_circles[counter].u_body;
				input_data >> M_circles[counter].v_body;
				input_data.getline(&buffer[0], MaxCharRead, '\n');

				input_data >> M_circles[counter].p_body;
				input_data.getline(&buffer[0], MaxCharRead, '\n');
				input_data >> M_circles[counter].T_body;
				input_data.getline(&buffer[0], MaxCharRead, '\n');

				input_data >> M_circles[counter].is_VelocitySlipBC_body;
				input_data.getline(&buffer[0], MaxCharRead, '\n');
				input_data >> M_circles[counter].F_VelocitySlip;
				input_data.getline(&buffer[0], MaxCharRead, '\n');
				input_data >> M_circles[counter].w_VelocitySlipBC;
				input_data.getline(&buffer[0], MaxCharRead, '\n');

				input_data >> M_circles[counter].is_TemperatureJumpBC_body;
				input_data.getline(&buffer[0], MaxCharRead, '\n');
				input_data >> M_circles[counter].F_TemperatureJump;
				input_data.getline(&buffer[0], MaxCharRead, '\n');
				input_data >> M_circles[counter].dTdn_on_wall_0;
				input_data.getline(&buffer[0], MaxCharRead, '\n');

				input_data >> M_circles[counter].ToSolveDragCoefficientForThisBody;
				input_data.getline(&buffer[0], MaxCharRead, '\n');
			}

			input_data.close();


			cout << "Data for circles was succesfully imported from file: "
				<< FileName << endl;

		}
		else
		{
			cout << "The program can't open file : "
				<< FileName << endl;
		}


	}


	template <typename T_FileNameEP>
	void Export_data_for_polyhedrons_to_file(const T_FileNameEP& FileName)
	/*	Target: to export data for polyhedrond to file
				with name file_name_for_polyhedrond
	*/
	{
		if(!is_enter_M_polyhedrons)
			cout << "There is no polyhedrons to export in file: "
				<< FileName << endl;
		else
		{
			ofstream output_data;
			output_data.open(FileName);

			output_data << N_polyhedrons << "	//Number of polyhedrons"  << endl;

			unsigned int counter1, counter2;

			for(counter1 = 0; counter1 < N_polyhedrons; counter1++)
			{
				output_data
					<< M_polyhedrons[counter1].N_vertex_body << "	//Number of vertexes of polyhedron" << endl

					<< M_polyhedrons[counter1].point1_inside_x << "	"
					<< M_polyhedrons[counter1].point1_inside_y << "	//point1_inside_x point1_inside_y" << endl
					<< M_polyhedrons[counter1].point2_inside_x << "	"
					<< M_polyhedrons[counter1].point2_inside_y << "	//point2_inside_x point2_inside_y" << endl

					<< M_polyhedrons[counter1].u_body << "	"
					<< M_polyhedrons[counter1].v_body << "	//u_body v_body" << endl

					<< M_polyhedrons[counter1].p_body << "	//pressure in polyhedron - must be 0" << endl
					<< M_polyhedrons[counter1].T_body << "	//Temperature in polyhedron" << endl

					<< M_polyhedrons[counter1].is_VelocitySlipBC_body << "	//is_VelocitySlipBC_body" << endl
					<< M_polyhedrons[counter1].F_VelocitySlip << "	//F_VelocitySlip" << endl
					<< M_polyhedrons[counter1].w_VelocitySlipBC << "	//w_VelocitySlipBC" << endl

					<< M_polyhedrons[counter1].is_TemperatureJumpBC_body << "	//is_TemperatureJumpBC_body" << endl
					<< M_polyhedrons[counter1].F_TemperatureJump << "	//F_TemperatureJump" << endl
					<< M_polyhedrons[counter1].dTdn_on_wall_0 << "	//dTdn_on_wall_0" << endl

					<< M_polyhedrons[counter1].ToSolveDragCoefficientForThisBody << "	//ToSolveDragCoefficientForThisBody" << endl;


				for(counter2 = 0; counter2 < M_polyhedrons[counter1].N_vertex_body; counter2++)
				{

					output_data <<  M_polyhedrons[counter1].x_coordinate[counter2] << "	"
						<<  M_polyhedrons[counter1].y_coordinate[counter2] << endl;
				}

			}

			output_data.close();


			cout << "Data for polyhedrons was succesfully exported to file: "
				<< FileName << endl;

		}

	}


    template <typename T_FileNameIP>
    void Import_data_for_polyhedrons_from_file(
                                                const T_FileNameIP,
                                                const unsigned int);


};

//contain informatin about body for solving, where pressure is given, where velocity is given.
extern Objects_for_Solving * OforS;


template <typename T_FileNameIP>
void Objects_for_Solving::Import_data_for_polyhedrons_from_file(
                                            const T_FileNameIP FileName,
                                            //const string FileName,
                                            const unsigned int gd)
/*	Target: to import data for polyhedrond to file
            with name file_name_for_polyhedrond
    Recieve:	FileName
                ga - number of given domain - gb - given body (int)
                                            gp - given pressure (int)
                                            gV - given Velocity (int)
*/
{
    using namespace std;

	const int MaxCharRead = 5000;
	char buffer[MaxCharRead];

    ifstream input_data;
    input_data.open(FileName);

    if(input_data.is_open())
    {
        input_data >> N_polyhedrons;
		input_data.getline(&buffer[0], MaxCharRead, '\n');

        make_new_M_polyhedrons();

        unsigned int counter1, counter2;

        for(counter1 = 0; counter1 < N_polyhedrons; counter1++)
        {
            input_data >> M_polyhedrons[counter1].N_vertex_body;
			input_data.getline(&buffer[0], MaxCharRead, '\n');
            M_polyhedrons[counter1].define_vector_for_vertexes();

            input_data >> M_polyhedrons[counter1].point1_inside_x;
            input_data >> M_polyhedrons[counter1].point1_inside_y;
			input_data.getline(&buffer[0], MaxCharRead, '\n');
            input_data >> M_polyhedrons[counter1].point2_inside_x;
            input_data >> M_polyhedrons[counter1].point2_inside_y;
			input_data.getline(&buffer[0], MaxCharRead, '\n');

            input_data >> M_polyhedrons[counter1].u_body;
            input_data >> M_polyhedrons[counter1].v_body;
			input_data.getline(&buffer[0], MaxCharRead, '\n');

            input_data >> M_polyhedrons[counter1].p_body;
			input_data.getline(&buffer[0], MaxCharRead, '\n');
            input_data >> M_polyhedrons[counter1].T_body;
			input_data.getline(&buffer[0], MaxCharRead, '\n');

			input_data >> M_polyhedrons[counter1].is_VelocitySlipBC_body;
			input_data.getline(&buffer[0], MaxCharRead, '\n');
			input_data >> M_polyhedrons[counter1].F_VelocitySlip;
			input_data.getline(&buffer[0], MaxCharRead, '\n');
			input_data >> M_polyhedrons[counter1].w_VelocitySlipBC;
			input_data.getline(&buffer[0], MaxCharRead, '\n');

			input_data >> M_polyhedrons[counter1].is_TemperatureJumpBC_body;
			input_data.getline(&buffer[0], MaxCharRead, '\n');
			input_data >> M_polyhedrons[counter1].F_TemperatureJump;
			input_data.getline(&buffer[0], MaxCharRead, '\n');
			input_data >> M_polyhedrons[counter1].dTdn_on_wall_0;
			input_data.getline(&buffer[0], MaxCharRead, '\n');

			input_data >> M_polyhedrons[counter1].ToSolveDragCoefficientForThisBody;
			input_data.getline(&buffer[0], MaxCharRead, '\n');


            for(counter2 = 0; counter2 < M_polyhedrons[counter1].N_vertex_body; counter2++)
            {

                input_data >> M_polyhedrons[counter1].x_coordinate[counter2];
                input_data >> M_polyhedrons[counter1].y_coordinate[counter2];
            }

        }

        input_data.close();

        bool data_are_correct = true;

        counter1 = 0;
        do
        {
            data_are_correct = OforS[gd].M_polyhedrons[counter1].is_data_for_polyhedron_correct();

            counter1++;
        }while(counter1 < N_polyhedrons && data_are_correct);


        if(data_are_correct)
            cout << "Data for polyhedrons: "
            << FileName << " are correct!" << endl;
        else
            cout << "WARNING: Data for polyhedrons: "
            << FileName << " are NOT correct!" << endl;

    }
    else
    {
        cout << "The program can't open  file : "
            << FileName << endl;
    }

}



#endif //ifndef Objects_for_Solving_H
