function [A2,c,s0,ncol] = restr_normal_eq(A,b,n)
%[A2,c,s0,ncol] = restr_normal_eq(A,b,n)

ncol = 0;
sn0 = zeros(n,1);
if (abs(b(1)) > 1e-10)
    if (b(1) < 0.0)
        E1 = [1;0];
    else
        E1 = [-1;0];
    end
    ncol = ncol+1;
    s0 = [sn0; abs(b(1))];
    sn0 = s0;
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
    s0 = [sn0; abs(b(2))];
else
    E2 = [];
end
A2 = [A E1 E2];
c = [zeros(n,1); ones(ncol,1)];