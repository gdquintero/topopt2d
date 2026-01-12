function bk = hess_diag_reciprocal(n,x,df0dx)

bk = -2*(df0dx./x);
for i = 1:n
    if (bk(i) < 0.0)
       bk(i) = 0.0;
       %bk(i) = bk(i) + 1.1*abs(bk(i));
    end
end




