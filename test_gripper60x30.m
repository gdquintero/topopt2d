function [xstar,u,itotal,itotalint,f0val,time] = test_gripper60x30

%profile clear;
%profile on;
load gripper60x30;
str = gripper60x30;
%x0 = 0.4*ones(1800,1);
%x0 = [0.1*ones(900,1);0.7*ones(900,1)];
%x0 = rand(1800,1);
x0 = 0.2*ones(1800,1);
xmin = 0.001;
delta = 0.1;
%delta = 0.01;
frmax = 0.2;
prnt = 2000;
opfil = 1;
rmin = 1.5;
uin = 2.0;
%uin = -5.0;
%uin = 0.5;
%ks = 0.004;
ks = 1.0;
ophole = 0;
load vazio60x30;
empty = vazio60x30;
%cd ..
[xstar,u,itotal,itotalint,f0val,time] = test_mechanisms(str,x0,xmin,delta,frmax,prnt,opfil,rmin,uin,ks,ophole,empty);
%profsave;
%profile off;
%cd C:\Users\thadeu\Desktop\SLP_Gomes_Senne_Matlab\Tests


