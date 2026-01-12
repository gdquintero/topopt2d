function x = gauss_seidel(A,b,x0,maxit,tol)

[mA,nA] = size(A);
if (mA ~= nA)
   disp('Matriz A deve ser uma matriz n x n.');
   return;
end

[mx0,nx0] = size(x0);
if ((mx0 ~= nA) || (nx0 ~= 1))
   disp('Vetor x0 deve ser uma matriz n x 1');
   return;
end

[mb,nb] = size(b);
if ((mb ~= nA) || (nb ~= 1))
   disp('Vetor b deve ser uma matriz n x 1');
   return;
end

k = 0;
L = tril(A,-1);
D = diag(diag(A));
U = triu(A,1);
res = b-A*x0;
normres = norm(res,inf);
rel = 1.0;

disp(' k     normres         rel    ');
disp('---  -----------   -----------');
while ((normres > tol) && (rel > tol) && (k < maxit))   
    x = (L+D)\(b-U*x0);
    res = b-A*x;
    normres = norm(res,inf);
    normx0 = norm(x0,inf);
    if (normx0 < tol)
       rel = norm(x-x0,inf);
    else
       rel = norm(x-x0,inf)/norm(x0,inf);
    end
    x0 = x;
    k = k+1;
    fprintf('%2u %13.5E %13.5E\n',k,normres,rel);
end
      


