function [f0val,rdesl,rvol,u,lambdaout,lambdain] = obj_func_mech(xfil,ein2,eout2,f2,K,perm,velem,Vmax,din,dout,uin)

%Getting the nodal displacements vector:
K2 = K(perm,perm);
R2 = chol(K2);
u = get_nodal_disp(R2,f2,perm);

%Objective function value at the current point:
%f0val = -u(dout);
f0val = u(dout);

%Adjoint vectors for the computation of the gradient vector of the
%objective function and of the gradient vector of the displacement
%constraint:
lambdaout = get_nodal_disp(R2,-eout2,perm);
lambdain = get_nodal_disp(R2,-ein2,perm);

%Displacement constraint at the point where the input force is applied:
rdesl = (u(din)-uin)/uin;

%Volume constraint value at the current point:
rvol = (velem'*xfil - Vmax)/Vmax;

   

