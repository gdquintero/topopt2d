function [f0val,vol,xnew,xvol,u] = obj_func_str(n,str,xfil,p,f,f2,K,perm,velem,Vmax,opfil,N,weigh,wi,numviz,nonempty,sinhp)
%***********************************************************************************************************************
%[f0val,vol,xnew,xvol,u] = obj_func_str(n,str,xfil,p,f,f2,K,perm,velem,Vmax,opfil,N,weigh,wi,numviz,nonempty,sinhp)
%
%This function computes the objective function and the volume constraint 
%values, at the current point.
%
%INPUT PARAMETERS:
%-----------------
%n: Number of elements of the rectangular domain.
%str: Structure with the problem data.
%xfil: Filtered element densities vector.
%p: Penalty parameter of the SIMP model for the element densities.
%f: Nodal equivalent forces vector.
%f2: Permuted nodal equivalent forces vector.
%K: Global stiffness matrix, stored as a sparse matrix.
%perm: Permutation vector associated to Cholesky factorization
%      of the global stiffness matrix.
%velem: Vector of element volumes. 
%Vmax: Maximum volume that the structure can contain.
%opfil: Parameter that indicates what filter is used.
%       opfil = 1 (mean density filter)
%       opfil = 2 (dilation filter)
%       opfil = 3 (erosion filter)
%       opfil = 4 (sinh filter)
%N: Matrix where each column contains the neighbors of each finite element.
%weigh: Matrix where each row contains the weights of each finite element,
%       associated to the mean density filter or to the sinh filter.
%       If the dilation or the erosion filter are used or 
%       if no filter is applied, "weigh" is an empty matrix.
%wi: Vector where each element is the sum of the finite element weights,
%    associated to the mean density filter or to the sinh filter.
%    If the dilation of erosion filter are used or 
%    if no filter is applied, "wi" is an empty vector.
%numviz: Vector that contains the number of neighbors of each finite element.
%nonempty: Vector with the indexes of the finite elements that are nonempty.
%sinhp: Constant used when the sinh filter is used.
%
%OUTPUT PARAMETERS:
%------------------
%f0val: Objective function value.
%vol: Volume constraint value.
%xnew: Filtered element densities vector that is used to plot the
%      structure, when the Sinh filter is applied. 
%      If we use another filter or if no filter is used, we have xnew = xfil.
%xvol: Filtered element densities vector, obtained by applying the mean
%      density filter. This vector is used only when the sinh filter is applied.
%      Otherwise, "xvol" is an empty vector.
%u: Nodal displacements vector.
%***********************************************************************************************************************

%Getting the nodal displacements vector:
K2 = K(perm,perm);
R2 = chol(K2);
u = get_nodal_disp(R2,f2,perm);

%Objective function value at the current point:
f0val = f'*u;

%Volume constraint value at the current point:
if opfil ~= 4  
   xvol = [];
   xnew = xfil;
else
   xvol = mean_density_filter(str,xfil,weigh,wi,nonempty,n,numviz,N);
   xnew = 1 - (sinh(p*(1-xvol))/sinhp); 
end
vol = (velem'*xnew - Vmax)/Vmax;
   

