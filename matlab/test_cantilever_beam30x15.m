function [xstar,u,itotal,itotalint,f0val,time] = test_cantilever_beam30x15

load cantilever_beam30x15;
str = cantilever_beam30x15;
x0 = 0.4*ones(450,1);
xmin = 0.001;
delta = 0.1;
frmax = 0.4;
prnt = 2000;
opfil = 1;
rmin = 1.5;
[xstar,u,itotal,itotalint,f0val,time] = test_structures(str,x0,xmin,delta,frmax,prnt,opfil,rmin);


