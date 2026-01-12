function [xstar,u,itotal,itotalint,f0val,time] = test_mechanisms(str,x0,xmin,delta,frmax,prnt,opfil,rmin,uin,ks,ophole,empty)

itotal = 0;
itotalint = 0;
xstar = x0;
tt = 0.0;
%p = 1.0;
p = 3.0;
deltap = 1.0;
tol = 1e-3;

while p <= 3
    if (p < 3)
        maxit = 50;
    else
        maxit = 100;
    end
    [xstar,u,iter,itint,f0val,time] = linseq_ir_mech(str,xstar,xmin,p,delta,frmax,tol,maxit,prnt,opfil,rmin,uin,ks,ophole,empty);
    %[xstar,u,iter,itint,f0val,time] = linseq_ir_mech_quadrest(str,xstar,xmin,p,delta,frmax,tol,maxit,prnt,opfil,rmin,uin,ks,ophole,empty);
    p = p + deltap;
    itotal = itotal + iter;
    itotalint = itotalint + itint;
    tt = tt + time;
end


time = tt;
