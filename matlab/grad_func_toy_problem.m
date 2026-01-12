function [df0dx,A] = grad_func_toy_problem(x)

df0dx = [2*x(1); 2*x(2)];
A = [2*x(1); -1];

%df0dx = [-2*x(1); -2*x(2)];
%A = [-1; 1];