function ubKuadx = grad_mutual_energy_density(str,xfil,p,ua,ub,kel,gradxnew,numviz,nonempty,kmax,N)
%********************************************************************************************
%ubKuadx = grad_mutual_energy_density(str,xfil,p,ua,ub,kel,gradxnew,numviz,nonempty,kmax,N)
%
%This function computes the gradient of mutual energy of 
%compliant mechanisms when the mean density filter is applied.
%
%*Nishiwaki, S.; Frecker, M. I.; Seungjae, M.; Kikuchi, N.; Topology
%optimization of compliant mechanisms using the homogenization method.
%International Journal for Numerical Methods in Engineering, 42, 1998, p. 535-559.
%
%INPUT PARAMETERS:
%-----------------
%str: Structure with the problem data.
%xfil: Filtered element densities vector, using the mean density filter.
%p: Penalty parameter of the SIMP model for the element densities.
%ua: Nodal displacements vector associated to the applied traction in 
%    the first load case (Nishiwaki et al. formulation*).
%ub: Nodal displacements vector associated to the dummy load in the
%    first load case (Nishiwaki et al. formulation*).
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
%ubKuadx: Mutual energy gradient vector compliant mechanisms,
%         when the mean density filter is applied.
%********************************************************************************************

n = str.nelem;
ubKuadxdens = zeros(n,1);
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
    uael = ua(lin);
    ubel = ub(lin);
    ubKuadxdens(m1) = -p*(xfil(m1)^(p-1))*(ubel'*(kel*uael));
end

%Gradient of the objective function in relation to the original densities, using the chain rule.
ubKuadx = zeros(n,1);
for i = 1:kmax
    m1 = nonempty(i);
    ind = numviz(m1);
    viz = N(1:ind,m1);
    ubKuadx(m1) = gradxnew(m1,viz)*ubKuadxdens(viz);
end

