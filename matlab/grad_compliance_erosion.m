function df0dx = grad_compliance_erosion(str,xnew,p,u,kel,gradxer,N,nonempty,numviz,kmax)
%***********************************************************************************************************************
%df0dx = grad_compliance_erosion(str,xnew,p,u,kel,gradxer,N,nonempty,numviz,kmax)
%
%This function computes the objective function gradient vector when 
%the erosion filter is applied.
%
%INPUT PARAMETERS:
%-----------------
%str: Structure with the problem data.
%xnew: Filtered element densities vector, using the erosion filter.
%p: Penalty parameter of the SIMP model for the element densities.
%u: Nodal displacements vector.
%kel: Element stiffness matrix.
%gradxer: Matrix where each line contains the coefficients associated to
%         the chain rule for the objective function derivatives,
%         when the erosion filter is used. 
%N: Matrix where each column contains the neighbors of each finite element.
%nonempty: Vector with the indexes of the finite elements that are
%          nonempty.
%numviz: Vector that contains the number of neighbors of each finite
%        element.
%kmax: Number of nonempty elements.
%
%OUTPUT PARAMETER:
%-----------------
%df0dx: Objective function gradient vector, when the erosion filter is applied.
%***********************************************************************************************************************

n = str.nelem;
ny = str.ny;
grdxer = zeros(n,1);
df0dx = zeros(n,1);
k = ny-1; 

%Gradient of the objective function in relation to the filtered densities:
for i = 1:kmax
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
    grdxer(m1) = -p*(xnew(m1)^(p-1))*(uel'*(kel*uel)); 
end

%Gradient of the objective function in relation to the original densities, using the chain rule.
for i = 1:kmax
    m1 = nonempty(i);
    viz = numviz(m1);
    ind = N(1:viz,m1);
    df0dx(m1) = gradxer(m1,ind)*grdxer(ind);
end