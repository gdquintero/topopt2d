function [xstar,u,itotal,itotalint,f0val,time] = test_inverter30x15

%profile clear;
%profile on;
load inverter30x15;
str = inverter30x15;
%x0 = 0.4*ones(1800,1);
%x0 = [0.1*ones(900,1);0.7*ones(900,1)];
%x0 = rand(1800,1);
x0 = 0.2*ones(450,1);
xmin = 0.001;
delta = 0.1;
%delta = 0.01;
frmax = 0.2;
prnt = 2000;
opfil = 1;
rmin = 1.5;
uin = 5.0;
%uin = 0.05;
%ks = 5000.0;
ks = 4.0;
ophole = 0;
empty = [];
%cd ..
[xstar,u,itotal,itotalint,f0val,time] = test_mechanisms(str,x0,xmin,delta,frmax,prnt,opfil,rmin,uin,ks,ophole,empty);
%profsave;
%profile off;
%cd C:\Users\thadeu\Desktop\SLP_Gomes_Senne_Matlab\Tests


