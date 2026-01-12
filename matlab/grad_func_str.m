function [df0dx,A] = grad_func_str(str,xfil,xdens,beta,p,u,kel,velem,Vmax,gradxnew,gradxdil,gradxer,numviz,nonempty,n,N,weigh,wi,opfil,csnhp,lin)
%**********************************************************************************************************************************
%[df0dx,A] = grad_func_str(str,xfil,xvol,p,u,kel,velem,Vmax,gradxnew,gradxdil,gradxer,numviz,nonempty,n,N,opfil,csnhp,lin)
%
%This function computes the objective function gradient vector (only for structures) and the
%volume constraint gradient vector, at the current point.
%
%INPUT PARAMETERS:
%-----------------
%str: Structure with the problem data.
%xfil: Filtered element densities vector.
%xvol: Filtered element densities vector, obtained by applying the mean
%      density filter. This vector is used only when the sinh filter is applied.
%      Otherwise, "xvol" is an empty vector.
%p: Penalty parameter of the SIMP model for the element densities, or 
%   penalty parameter of the sinh filter for the volume constraint.
%u: Nodal displacements vector.
%kel: Element stiffness matrix.
%velem: Element volumes vector. 
%Vmax: Maximum volume that the structure can contain.
%gradxnew: Matrix where each line contains the coefficients associated to
%          the chain rule for the objective function derivatives,
%          when the mean density filter is used. 
%          If the dilation or the erosion filter are used, "gradxnew" is an empty
%          matrix.
%gradxdil: Matrix where each line contains the coefficients associated to
%          the chain rule for the objective function derivatives,
%          when the dilation filter is used. 
%          If another filter is used, "gradxdil" is an empty matrix.
%gradxer: Matrix where each line contains the coefficients associated to
%         the chain rule for the objective function derivatives,
%         when the erosion filter is used. 
%         If another filter is used, "gradxer" is an empty matrix.
%numviz: Vector that contains the number of neighbors of each finite element.
%nonempty: Vector with the indexes of the finite elements that are
%          nonempty.
%n: Number of elements of the rectangular domain.
%N: Matrix where each column contains the neighbors of each finite element.
%opfil: Parameter that indicates what filter is used.
%       opfil = 1 (mean density filter)
%       opfil = 2 (dilation filter)
%       opfil = 3 (erosion filter)
%       opfil = 4 (sinh filter)
%csnhp: Constant used when the sinh filter is used.
%lin: matrix of order 8 x nelem where the j-th column contains
%     the indexes of the degrees of freedom associated to the j-th
%     element.
%
%OUTPUT PARAMETERS:
%------------------
%df0dx: Objective function gradient vector.
%A: Volume constraint gradient vector. 
%***********************************************************************************************************************


if opfil == 1 %Mean density filter (opfil == 1)
   df0dx = grad_compliance_mean_density(str,xfil,p,u,kel,gradxnew,numviz,nonempty,n,N,lin);
   A = [];
elseif opfil == 2 %Heaviside filter (opfil == 2)
   df0dx = grad_compliance_heaviside(str,xfil,xdens,weigh,wi,p,u,kel,gradxnew,numviz,nonempty,beta,n,N,lin);
   A = grad_volume_heaviside(str,xfil,xdens,velem,N,weigh,wi,beta,nonempty,n,Vmax,numviz)';
end

% if opfil == 0 %No filter (opfil == 0)
%    df0dx = grad_compliance(str,n,xfil,p,u,kel,nonempty);
%    A = [];
% elseif opfil == 1 %Mean density filter (opfil == 1)
%    df0dx = grad_compliance_mean_density(str,xfil,p,u,kel,gradxnew,numviz,nonempty,n,N,lin);
%    A = [];
% elseif opfil == 2 %Dilation filter (opfil == 2)
%    df0dx = grad_compliance_dilation(str,xfil,p,u,kel,gradxdil,N,nonempty,numviz,n);
%    A = grad_volume_dilation(str,gradxdil,velem,Vmax,nonempty,numviz,n,N)';
% elseif opfil == 3 %Erosion filter (opfil == 3)
%    df0dx = grad_compliance_erosion(str,xfil,p,u,kel,gradxer,N,nonempty,numviz,n);
%    A = grad_volume_erosion(str,gradxer,velem,Vmax,nonempty,numviz,n,N)';
% elseif opfil == 4 %Sinh filter (opfil == 4)
%    df0dx = grad_compliance(str,n,xfil,1,u,kel,nonempty);
%    A = grad_volume_sinh(str,xvol,p,velem,numviz,N,nonempty,n,csnhp,gradxnew)';
% end

