function [xstar,u,itotal,itotalint,f0val,time] = test_plate80x20_2
%***********************************************************************************************************************
%[xstar,itotal,itotalint,f0val,time] = test_cantilever_beam_100x25
%
%This function calls the main programs (test_structures.m and linseq_struct.m) to solve the topology
%optimization problem to find the optimum topology for a cantilever beam, using 
%a globally convergent version of the Sequential Linear Programming (SLP) algorithm*, 
%proposed by Gomes and Senne.
%
%*Gomes, F. A. M.; Senne, T. A.; An SLP algorithm and its application to
%topology optimization. Computational and Applied Mathematics, Vol. 30, N.
%1, pp. 53-89, 2011.
%
%
%INPUT PARAMETERS:
%-----------------
%str: Structure with the problem data.
%x0: Initial point for the SLP algorithm (initial element densities).
%xmin: Small positive value that the element densities can assume, to avoid numerical instabilities.
%delta: Initial trust region radius for the SLP algorithm.
%frmax: Volume fraction of the rectangular domain that the optimal
%       structure can contain. 
%tol: Minimum value for the step norm obtained by the SLP algorithm (stopping criterion).
%maxit: Maximum number of iterations that the algorithm can perform (stopping criterion).
%prnt: Number of times that the outer iteration data and the topology will be
%      printed on the screen. 
%opfil: Parameter that indicates what filter is used.
%       opfil = 0 (no filter)
%       opfil = 1 (mean density filter)
%       opfil = 2 (dilation filter)
%       opfil = 3 (erosion filter)
%       opfil = 4 (sinh filter)
%rmin: Filter radius.
%
%OUTPUT PARAMETERS:
%------------------
%xstar: Approximated optimal solution of the topology optimization problem.
%itotal: Total number of outer iterations performed by the SLP algorithm to find
%        the approximated optimal solution of the topology optimization
%        problem.
%itotalint: Total number of inner iterations performed by the SLP algorithm to find
%           the approximated optimal solution of the topology optimization
%           problem.
%f0val: Approximated optimal value of the objective function.
%f0val2: Approximated optimal value of the objective function, obtained setting p = 1 and without using any filter
%        at the approximated optimal solution.
%time: Total time taken by the SLP algorithm to find the approximated
%      optimal solution of the topology optimization problem. 
%***********************************************************************************************************************

%profile clear;
%profile on;
load plate80x20_2;
str = plate80x20_2;
x0 = 0.25*ones(2500,1);
xmin = 0.001;
delta = 0.1;
frmax = 0.25;
prnt = 2000;
opfil = 1;
rmin = 2.0;
cd ..
[xstar,u,itotal,itotalint,f0val,time] = test_structures(str,x0,xmin,delta,frmax,prnt,opfil,rmin);
%profsave;
%profile off;
%cd C:\Users\thadeu\Desktop\SLP_Gomes_Senne_Matlab\Tests


