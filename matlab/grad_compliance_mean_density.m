function df0dx = grad_compliance_mean_density(str,xnew,p,u,kel,gradxnew,numviz,nonempty,kmax,N,lin)
%***********************************************************************************************************************
%df0dx = grad_compliance_mean_density(str,xnew,p,u,kel,gradxnew,numviz,nonempty,kmax,N,lin) 
%
%This function computes the objective function gradient vector when 
%the mean density filter is applied.
%
%INPUT PARAMETERS:
%-----------------
%str: Structure with the problem data.
%xnew: Filtered element densities vector, using the mean density filter.
%p: Penalty parameter of the SIMP model for the element densities.
%u: Nodal displacements vector.
%kel: Element stiffness matrix.
%gradxnew: Matrix where each line contains the coefficients associated to
%          the chain rule for the objective function derivatives,
%          when the mean density filter is used. 
%numviz: Vector that contains the number of neighbors of each finite
%        element.
%nonempty: Vector with the indexes of the finite elements that are
%          nonempty.
%kmax: Number of nonempty elements.
%N: Matrix where each column contains the neighbors of each finite element.
%lin: matrix of order 8 x nelem where the j-th column contains
%     the indexes of the degrees of freedom associated to the j-th
%     element.
%
%OUTPUT PARAMETER:
%-----------------
%df0dx: Objective function gradiente vector, when the mean density filter is applied.
%***********************************************************************************************************************

n = str.nelem;
gradxdens = zeros(n,1);

%Gradient of the objective function in relation to the filtered densities:
for i = 1:kmax
    m1 = nonempty(i);
    uel = u(lin(:,m1));
    gradxdens(m1) = -p*(xnew(m1)^(p-1))*(uel'*(kel*uel));
end
%kel*uel

%Gradient of the objective function in relation to the original densities, using the chain rule.
df0dx = zeros(n,1);
for i = 1:kmax
    m1 = nonempty(i);
    ind = numviz(m1);
    viz = N(1:ind,m1);
    df0dx(m1) = gradxnew(m1,viz)*gradxdens(viz);
end

