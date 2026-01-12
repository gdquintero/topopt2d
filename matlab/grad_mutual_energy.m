function ubKuadx = grad_mutual_energy(str,n,x,p,ua,ub,kel,nonempty,kmax)
%**************************************************************************
%ubKuadx = grad_mutual_energy(str,n,x,p,ua,ub,kel,nonempty,kmax)
%
%This function computes the mutual energy gradient vector for 
%compliant mechanisms when no filter is applied.
%
%*Nishiwaki, S.; Frecker, M. I.; Seungjae, M.; Kikuchi, N.; Topology
%optimization of compliant mechanisms using the homogenization method.
%International Journal for Numerical Methods in Engineering, 42, 1998, p.
%535-559.
%
%INPUT PARAMETERS:
%-----------------
%str: Structure with the problem data.
%n: Number of elements of a rectangular domain.
%x: Original element densities vector.
%p: Penalty parameter of the SIMP model for the element densities.
%ua: Nodal displacements vector associated to the applied traction in 
%    the first load case (Nishiwaki et al. formulation*).
%ub: Nodal displacements vector associated to the dummy load in the
%    first load case (Nishiwaki et al. formulation*).
%kel: Element stiffness matrix.
%nonempty: Vector with the indexes of the finite elements that are
%          nonempty.
%kmax: Number of nonempty elements.
%
%OUTPUT PARAMETER:
%-----------------
%ubKuadx: Mutual energy gradient vector for compliant mechanisms,
%         when no filter is applied.
%**************************************************************************

ubKuadx = zeros(n,1);
lin = zeros(8,1);
ny = str.ny;
k = ny-1; 

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
    ubKuadx(m1) = -p*(x(m1)^(p-1))*(ubel'*(kel*uael));
end