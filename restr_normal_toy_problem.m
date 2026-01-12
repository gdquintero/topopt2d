function [A2,c,s0] = restr_normal_toy_problem(A,b)

if (abs(b) > 1e-10)
   if (b < 0.0)
      A2 = [A' 1];
   else
      A2 = [A' -1];
   end
end
c = [0; 0; 1];
s0 = [0; 0; abs(b)];
