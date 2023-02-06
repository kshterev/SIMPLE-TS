# SIMPLE-TS
SIMPLE-TS is a finite volume method for the calculation of unsteady, viscous, compressible, and heat-conductive flows.
It can be reduced straightforwardly to calculate unsteady, viscous, incompressible, and isothermal flows.
SIMPLE-TS is part of SIMPLE-like algorithms. SIMPLE-TS calculates numerically equations of the Navier-Stokes.

In SIMPLE-TS for the first time is substituted density in the unsteady term in the pressure equation with pressure
and temperature using the equation of state. In this way, we turned the numerical pressure equation into a stable
that does not need a relaxation coefficient to ensure convergence.
Furthermore, SIMPLE-TS is five times faster than SIMPLE and slightly faster than PISO.

A derivation of the SIMPLE-TS method and considerations are available here: https://kirilshterev.com/index.php/finite-volume-method-for-the-compressible-viscous-gas-flows-simple-ts-derivation/
The web page contains:
 - detailed derivation of the numerical equations and explanations;
 - detailed explanation of differences between SIMPLE-TS and other SIMPLE-like methods considering an example of one-dimensional unsteady, isothermal pressure-driven flow in a duct;
 - derivation of the numerical equations of the algorithm SIMPLE-TS using a first-order upwind scheme for the approximation of convective terms with Mathematica;
 - derivation of the numerical equations of the algorithm SIMPLE-TS using a second-order total variation diminishing (TVD) scheme for the approximation of convective terms with Mathematica;

SIMPLE-TS source codes and examples are placed here: https://kirilshterev.com/index.php/finite-volume-method-for-the-compressible-viscous-gas-flows-source-codes-and-examples/
The web page contains:
 - parallel C++ source code with first-order approximation of convective terms;
 - unsteady supersonic, compressible, viscous, heat-conductive fluid flow past a confined square in a micro-channel – Mach number 2.43 and Knudsen number 0.00283 (Reynolds number 1415);
 - unsteady subsonic, compressible, viscous, heat-conductive fluid flow past a confined square in a micro-channel – Mach number 0.1 and Knudsen number 0.00194 (Reynolds number 85);
 - increasing velocity at the channel inflow from Mach number 2.43 to Mach number 4.86, for Knudsen number 0.05;
 - Rayleigh-Bènard flow of a rarefied gas;

The main work presented in this section is published in the following papers:
    K. Shterev and S. Stefanov, Pressure based finite volume method for calculation of compressible viscous gas flows, Journal of Computational Physics 229 (2010) pp. 461-480,  doi:10.1016/j.jcp.2009.09.042 – the accepted manuscript can be downloaded from here, the paper in its final mode is available here. IF 3.023
    K. S. Shterev and S. K. Stefanov, A Parallelization of Finite Volume Method for Calculation of Gas Microflows by Domain Decomposition Methods, 7th International Conference – Large-ScaleScientific Computat
    K. S. Shterev, S. K. Stefanov, Finite Volume Calculations of Gaseous Slip Flow Past a Confined Square in a Two-Dimensional Microchannel, Proceedings of the 1st European Conference on Microfluidics – Microfluidics 2008 – Bologna, December 10-12, 2008.
    K. S. Shterev, S. K. Stefanov, E.Atanasov, A parallel algorithm with improved performance of Finite Volume Method (SIMPLE-TS), 8th International Conference – Large-ScaleScientific Computations, Sozopol, Bulgaria, June 06-10, 2011.
    K. S. Shterev and S. Ivanovska, Comparison of some approximation schemes for convective terms for solving gas flow past a square in a microchannel, APPLICATION OF MATHEMATICS IN TECHNICAL AND NATURAL SCIENCES: 4th International Conference – AMiTaNS ’12, 11-16 June 2012, St. Constantine and Helena, Bulgaria, AIP Conf. Proc. 1487, pp. 79-87; doi:http://dx.doi.org/10.1063/1.4758944, ISBN 978-0-7354-1099-2, SJR 0.168. The original publication is available here. The authors’ copy is available free from here.

