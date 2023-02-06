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


//: C06:stringConv.h
// From "Thinking in C++, 2nd Edition, Volume 2"
// by Bruce Eckel & Chuck Allison, (c) 2001 MindView, Inc.
// Available at www.BruceEckel.com.
// Chuck Allison's string converter
#ifndef STRINGCONV_H
#define STRINGCONV_H
#include <string>
#include <sstream>

template<typename T>
T fromString(const std::string& s) {
  std::istringstream is(s);
  T t;
  is >> t;
  return t;
}

template<typename T>
std::string toString(const T& t) {
  std::ostringstream s;
  s << t;
  return s.str();
}

////my functions
//template<typename string_T, typename n_T>
//string_T intToString(/*const */n_T n)
///*
//Target:		to convert integer to string
//Receive:	n - ineger number
//Return:		Integer number in string format
//*/
//{
//	const unsigned int size_buffer = 50;
//
//	char buffer[size_buffer];
//	_itoa(n, buffer, 10);
//
//	return buffer;
//}
//
//template<typename n_T>
//string intTo_string(const n_T& n)
//
///*
//Target:		to convert integer to string
//Receive:	n - ineger number
//Return:		Integer number in string format
//*/
//{
//	const unsigned int size_buffer = 50;
//
//	char buffer[size_buffer];
//	_itoa(n, buffer, 10);
//
//	return buffer;
//}


#endif // STRINGCONV_H ///:~
