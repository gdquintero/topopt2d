function [lowplin, upplin] = get_plpbounds(n,gamma,x,xmin,breakpoints0,breakpoints1)

sL = max(-gamma, xmin-x);
sU = min(gamma, 1-x);
lowplin = zeros(3*n,1);
upplin = lowplin;

for i = 1:n
    lowplin(i) = sL(i);
    lowplin(i+n) = 0.0;
    lowplin(i+2*n) = 0.0;
    upplin(i) = breakpoints0(i);
    upplin(i+n) = breakpoints1(i)-breakpoints0(i);
    upplin(i+2*n) = sU(i)-breakpoints1(i);
end

