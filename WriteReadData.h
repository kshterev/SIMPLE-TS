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


//File name: WriteReadData.h

#ifndef WriteReadData_H
#define WriteReadData_H



#include <fstream> //This is for writing and reading data from file
#include <iomanip> //This is for formatting when writing data to file

#include <iostream>
using namespace std;


//Read massive from binary file.
template <typename T_FileName, typename T_massive>
void ReadMassiveFromBinaryFile(
		const T_FileName& FileName, T_massive massive[],
		unsigned int& size)
/*
Target:		Read massive from binary file.
Receive:	FileName - name of file from where data will be readed
			massive - massive which will be readed from file
			size - number of elements in massive
*/
{
	using namespace std;

	ifstream input_data;
	input_data.open(FileName, ios::in | ios::binary);

	if (input_data.is_open())
	{
		input_data.read((char *)& massive[0], sizeof(T_massive) * size);
		input_data.close();
	}
	else
	{
		cout << "The program can not open file " << FileName << " to read data." << endl;
	}
}


//To read data (massive) from file.
template <typename T_FileName, typename T_massive>
void ReadMassiveFromFile_1D(
	const T_FileName& FileName, T_massive massive[],
	const unsigned int& size)
/*
Target:		To read data (massive) from file.
Receive:	FileName - name of file from where data will be readed
			massive - massive which will be readed from file
			size - number of elements of massive
Node:		This procedure read massive from file like 1D in ASCII code.
*/
{
	using namespace std;

	ifstream input_data;
	input_data.open(FileName);

	if (input_data.is_open())
	{
		unsigned int counter;

		for (counter = 0; counter < size; counter++)
			input_data >> massive[counter];			

		input_data.close();
	}
	else
	{
		cout << "The program can not open file " << FileName << " to read data." << endl;
	}
}


//To write data (massive) to file.
template <typename T_FileName, typename T_massive>
void WriteMassiveToBinaryFile(
	const T_FileName& FileName, const T_massive massive[],
	const unsigned int& size)
/*
Target:		To write data (massive) to file.
Receive:	FileName - name of file where data will be written
			Extension - extension of file
			massive - massive which will be writted to file
			size - number of elements of massive
Node:		This procedure write massive to file in bunary mode.
*/
{
	using namespace std;

	ofstream output_data;
	output_data.open(FileName, ios::out | ios::trunc | ios::binary);

	if(output_data.is_open())
	{
		output_data.write((char *)& massive[0], sizeof(T_massive) * size);
		output_data.close();
	}
	else
	{
		cout << "The program can not open file " << FileName << " to write data." << endl;
	}

}


//To write data (massive) to file.
template <typename T_string, typename T_massive>
void WriteMassiveToFile_1D(
		const T_string& FileName, const T_massive massive[],
		const unsigned int& size, const unsigned int& Ndigits = 15)
/*
Target:		To write data (massive) to file.
Receive:	FileName - name of file where data will be written
			massive - massive which will be writted to file
			size - number of elements of massive
Node:		This procedure write massive from file like 1D in ASCII code.
*/
{
	using namespace std;

	ofstream output_data;
	output_data.open(FileName);

	if(output_data.is_open())
	{
		unsigned int counter;

		for (counter = 0; counter < size; counter++)
		{
			output_data << std::setiosflags(std::ios::scientific)
						<< std::setprecision(Ndigits)
						<< massive[counter]
						<< endl;
		}

		output_data.close();
	}
	else
	{
		cout << "The program can not open file " << FileName << " to write data." << endl;
	}
}


//To add line to and of text file without delete file.
template <typename T_FileName, typename T_Line>
void AddLineToTextFile(
	const T_FileName& FileName, const T_Line& AddLine)
/*
Receive:	FileName - name of file where data will be added
			AddLine - CString which will be added
Node:		This procedure write string to end of text file.
*/
{
	using namespace std;

	ofstream output_data;
	output_data.open(FileName, ios::out | ios::app);

	if(output_data.is_open())
	{
		output_data << AddLine << endl;

		output_data.close();
	}
	else
	{
		cout << "The program can not open file " << FileName << " to write data." << endl;
	}

}


template <typename T_string, typename T_massive>
void WriteMassiveToFile_2D(
		const T_string& FileName, const T_massive massive[],
		const unsigned int& Nx, const unsigned int& Ny,
		const unsigned int& width_between_columns = 30,
		const unsigned int& Ndigits = 15)
/*
Target:		To write data (massive) to file
Receive:	file_name - name of file where data will be written
			Extension - extension of file
			massive - massive which will be writted to file
			Nx - number of nodes on OX axis
			Ny - number of nodes on OY axis
			width_between_columns - width between columns of data in file
Node:		This procedure write massive from floating points data
			like 2D, the massive in program looks diferent,
			like vector - 1D
*/
{
	using namespace std;


	ofstream output_data;
	output_data.open(FileName);

	if(output_data.is_open())
	{
		unsigned int i, j, node_a;
		for (i = 0; i < Nx; i++)
		{
			for (j = 0; j < Ny; j++)
			{
				node_a = i + j * Nx;
				
				if  (j != 0) 
					output_data << std::setiosflags(std::ios::scientific)
						<< std::setw(width_between_columns)
						<< std::setprecision(Ndigits)
						<< massive[node_a];
				else output_data << std::setiosflags(std::ios::scientific)
					<< std::setw(width_between_columns)
					<< std::setprecision(Ndigits)
					<< massive[node_a];
			}
			
			if (i < (Nx - 1)) output_data << endl;
		}

		output_data.close();
	}
	else
	{
		cout << "The program can not open file " << FileName << " to write data." << endl;
	}

}



#endif //ifndef WriteReadData_H
