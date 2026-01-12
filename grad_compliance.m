function df0dx = grad_compliance(str,n,x,p,u,kel,nonempty)
%**************************************************************************
%df0dx = grad_compliance(str,n,x,p,u,kel,nonempty)
%
%This function computes the objective function gradient vector for a 
%structure when no filter is applied.
%
%INPUT PARAMETERS:
%-----------------
%str: Structure with the problem data.
%n: Number of elements of a rectangular domain.
%x: Original element densities vector.
%p: Penalty parameter of the SIMP model for the element densities.
%u: Nodal displacements vector.
%kel: Element stiffness matrix.
%nonempty: Vector with the indexes of the finite elements that are nonempty.
%
%OUTPUT PARAMETER:
%-----------------
%df0dx: Objective function gradient vector, when no filter is applied.
%**************************************************************************

df0dx = zeros(n,1);
lin = zeros(8,1);
ny = str.ny;
k = ny-1; 

for i = 1:n
    m1 = nonempty(i);
    lin(2) = 2*((floor((m1-1)/k))*ny+mod(m1-1,k)+1); 
    lin(4) = lin(2) + 2*ny;
    lin(6) = lin(4) + 2;
    lin(8) = lin(2) + 2;
    lin(1) = lin(2) - 1;
    lin(3) = lin(4) - 1;
    lin(5) = lin(6) - 1;
    lin(7) = lin(8) - 1;
    uel = u(lin);
    df0dx(m1) = -p*(x(m1)^(p-1))*(uel'*(kel*uel));
end