function [xstar,itotal,itotalint,f0val,f0val2,time] = test_force_inverter_60x30
%***************************************************************************************************************************************
%[xstar,itotal,itotalint,f0val,f0val2,time] = test_force_inverter_60x30
%
%This function calls the main programs (test_mechanisms.m and linseq_mech.m) to solve the topology
%optimization problem for compliant mechanisms, following the Nishiwaki et al. formulation* 
%and using a globally convergent version of the Sequential Linear Programming (SLP) algorithm**, 
%proposed by Gomes and Senne. Here, we want to find the optimum topology
%for a force inverter.
%
%*Nishiwaki, S.; Frecker, M. I.; Seungjae, M.; Kikuchi, N.; Topology
%optimization of compliant mechanisms using the homogenization method.
%International Journal for Numerical Methods in Engineering, 42, 1998, p. 535-559.
%
%**Gomes, F. A. M.; Senne, T. A.; An SLP algorithm and its application to
%topology optimization. Computational and Applied Mathematics, Vol. 30, N. 1, pp. 53-89, 2011.
%
%INPUT PARAMETERS:
%-----------------
%str1: Structure with the problem data, considering only the applied traction  
%      and the respective boundary conditions in the first load case (Nishiwaki et al. formulation*).
%str2: Structure with the problem data, considering only the dummy load  
%      and the respective boundary conditions in the first load case (Nishiwaki et al. formulation*).
%str3: Structure with the problem data, considering only the reaction force
%      and the respective boundary conditions in the second load case (Nishiwaki et al. formulation*).
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
%hole: Parameter that indicates if the rectangular domain contains some hole.
%      If hole = 1, the domain has a rectangular hole.
%      If hole = 0, the domain doesn't have any hole. 
%empty: Vector that contains the indexes of the elements that belongs to
%       the rectangular hole. If the domain doesn't have any hole, set
%       empty = [].
%
%OUTPUT PARAMETERS:
%------------------
%xstar: Approximated optimal solution of the topology optimization problem.
%iter: Total number of outer iterations performed by the SLP algorithm to find
%      the approximated optimal solution of the topology optimization
%      problem.
%itint: Total number of inner iterations performed by the SLP algorithm to find
%       the approximated optimal solution of the topology optimization
%       problem.
%f0val: Approximated optimal value of the objective function.
%f0val2: Approximated optimal value of the objective function, obtained setting p = 1 and without using any filter
%        at the approximated optimal solution.
%beta: Parameter used in the dilation and erosion filters. 
%time: Total time taken by the SLP algorithm to find the approximated
%      optimal solution of the topology optimization problem. 
%***************************************************************************************************************************************

load force_inverter_60x30_1;
load force_inverter_60x30_2;
load force_inverter_60x30_3;
str1 = force_inverter_60x30_1;
str2 = force_inverter_60x30_2;
str3 = force_inverter_60x30_3;
x0 = 0.2*ones(1800,1);
xmin = 0.001;
delta = 0.1;
frmax = 0.2;
tol = 1e-3;
maxit = 5000;
prnt = 5000;
opfil = 1;
rmin = 1.5;
hole = 0;
empty = [];
cd ..
[xstar,itotal,itotalint,f0val,f0val2,time] = test_mechanisms(str1,str2,str3,x0,xmin,delta,frmax,tol,maxit,prnt,opfil,rmin,hole,empty);
%cd C:\Users\thadeu\Desktop\SLP_Gomes_Senne\Tests