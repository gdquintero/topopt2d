function [f0val,b] = func_toy_problem(x)

f0val = x(1)^2 + x(2)^2;
b = x(1)^2 - x(2) + 1;

%f0val = -x(1)^2-x(2)^2;
%b = -x(1)+x(2)-1;