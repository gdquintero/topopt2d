function [A2,c,s0] = vol_constr_normal(A,b,n,sigma)
%[A2,c,s0] = vol_constr_normal(A,b,n,sigma)

sn0 = zeros(n,1);
if (b < 0.0)
   E1 = 1.0;
else
   E1 = -1.0;
end
s0 = [sn0; 0.0; abs(b)];
A2 = [A 0.0 E1];
c = [zeros(n,1); sigma; 1.0];
