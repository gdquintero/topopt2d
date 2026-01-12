function u = get_nodal_disp(R,f,perm)
%**************************************************************************
%u = get_nodal_disp(R,f,perm)
%
%This function returns the nodal displacements vector.
%
%INPUT PARAMETERS:
%-----------------
%R: Cholesky factor of the permuted global stiffness matrix.
%f: Nodal equivalent forces vector.
%perm: Permutation vector associated to Cholesky factorization
%      of the global stiffness matrix.
%
%OUTPUT PARAMETERS:
%------------------
%u: Nodal displacements vector.
%**************************************************************************

%Solving the linear system (R'*R)u = f:
w = (R')\f;
y = R\w;
u(perm,:) = y;

