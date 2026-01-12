function [A2,c,s0,ncol] = restr_normal(A,b,n)
%[A2,c,s0,ncol] = restr_normal(A,b,n)

ncol = 0;
s0 = zeros(n+2,1);
if (abs(b(1)) > 1e-10)
    if (b(1) < 0.0)
        E1 = [1;0];
    else
        E1 = [-1;0];
    end
    ncol = ncol+1;
    s0 = [s0; abs(b(1))];
else
    E1 = [];
end
if (abs(b(2)) > 1e-10)
    if (b(2) < 0.0)
        E2 = [0;1];
    else
        E2 = [0;-1];
    end
    ncol = ncol+1;
    s0 = [s0; abs(b(2))];
else
    E2 = [];
end
A2 = [A E1 E2];
c = [zeros(n+2,1); ones(ncol,1)];