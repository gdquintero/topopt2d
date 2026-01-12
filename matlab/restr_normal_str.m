function [A2,c,s0] = restr_normal_str(A,b,n)
%[A2,c,s0,ncol] = restr_normal(A,b,n)

sn0 = zeros(n,1);
if (abs(b) > 1e-10)
    if (b < 0.0)
        E1 = 1;
    else
        E1 = -1;
    end
    s0 = [sn0; abs(b(1))];
else
    E1 = [];
end
A2 = [A E1];
c = [zeros(n,1); 1];