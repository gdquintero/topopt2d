function [xstar,u,itotal,itotalint,f0val,time] = test_cb120x30_F120

load cb120x30_F120;
str = cb120x30_F120;
x0 = 0.5*ones(3600,1);
xmin = 0.001;
delta = 0.01;
frmax = 0.5;
prnt = 2000;
opfil = 1;
rmin = 1.5;
[xstar,u,itotal,itotalint,f0val,time] = test_structures(str,x0,xmin,delta,frmax,prnt,opfil,rmin);

