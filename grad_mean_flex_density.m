function ucKucdx = grad_mean_flex_density(str,xnew,p,uc,kel,gradxnew,numviz,nonempty,kmax,N)
%*************************************************************************************
%ucKucdx = grad_mean_flex_density(str,xnew,p,uc,kel,gradxnew,numviz,nonempty,kmax,N)
%
%This function computes the mean flexibility gradient vector for 
%compliant mechanisms when the mean density filter is applied.
%
%*Nishiwaki, S.; Frecker, M. I.; Seungjae, M.; Kikuchi, N.; Topology
%optimization of compliant mechanisms using the homogenization method.
%International Journal for Numerical Methods in Engineering, 42, 1998, p. 535-559.
%
%INPUT PARAMETERS:
%-----------------
%str: Structure with the problem data.
%xnew: Filtered element densities vector, using the mean density filter.
%p: Penalty parameter of the SIMP model for the element densities.
%uc: Nodal displacements vector associated to the reaction force in the
%    second load case (Nishiwaki et al. formulation*).
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
%
%OUTPUT PARAMETER:
%-----------------
%ucKucdx: Mean flexibility gradient vector for compliant mechanisms, 
%         when the mean density filter is applied.
%*************************************************************************************

n = str.nelem;
ucKucdxdens = zeros(n,1);
lin = zeros(8,1);
ny = str.ny;
k = str.ny-1;

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
    ucel = uc(lin);
    ucKucdxdens(m1) = -p*(xnew(m1)^(p-1))*(ucel'*(kel*ucel));
end

%Gradient of the objective function in relation to the original densities, using the chain rule.
ucKucdx = zeros(n,1);
for k = 1:kmax
    m1 = nonempty(k);
    ind = numviz(m1);
    viz = N(1:ind,m1);
    ucKucdx(m1) = gradxnew(m1,viz)*ucKucdxdens(viz);
end
