function bk = hess_diag_exponential(n,xold,x,df0dxold,df0dx)

bk = zeros(n,1);

for i = 1:n
    if ((df0dxold(i)/df0dx(i) > 1e-6) && (abs(x(i)-xold(i)) > 1e-6))
       aik = 1 + (log(df0dxold(i)/df0dx(i)))/(log(xold(i)/x(i)));
    else
       aik = -1;
    end
    bk(i) = ((aik-1)/x(i))*df0dx(i);
end
%bk = ((aik-1)./x).*df0dx;
       
for i = 1:n
    if (bk(i) < 0.0)
       bk(i) = 0.0;
       %bk(i) = bk(i) + 1.1*abs(bk(i));
    end
end
