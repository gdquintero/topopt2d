function gradvol = grad_volume_sinh(str,xvol,p,velem,numviz,N,nonempty,kmax,csnhp,gradxnew)
%***********************************************************************************************************************
%gradvol = grad_volume_sinh(str,xvol,p,velem,numviz,N,nonempty,kmax,csnhp,gradxnew)
%
%This function computes the volume constraint gradient vector, when 
%the sinh filter is applied.
%
%INPUT PARAMETERS:
%-----------------
%str: Structure with the problem data.
%xvol: Filtered element densities vector, obtained by applying the mean
%      density filter. 
%p: Penalty parameter of the sinh filter for the volume constraint.
%velem: Element volumes vector. 
%numviz: Vector that contains the number of neighbors of each finite
%        element.
%N: Matrix where each column contains the neighbors of each finite element.
%nonempty: Vector with the indexes of the finite elements that are
%          nonempty.
%kmax: Number of nonempty elements.
%csnhp: Constant used when the sinh filter is used.
%gradxnew: Matrix where each line contains the coefficients associated to
%          the chain rule for the volume constraint derivatives,
%          when the mean density filter is used. 
%
%OUTPUT PARAMETER:
%-----------------
%gradvol: Volume constraint gradient vector, when the sinh filter is applied.
%***********************************************************************************************************************

gradvol = zeros(str.nelem,1);
for k = 1:kmax
    m1 = nonempty(k);
    viz = numviz(m1);
    ind = N(1:viz,m1);
    gradvol(m1) = csnhp*(velem(ind).*cosh(p*(1-xvol(ind))))'*(p*gradxnew(m1,ind)'); 
end
