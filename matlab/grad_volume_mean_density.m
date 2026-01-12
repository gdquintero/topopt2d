function gradvol = grad_volume_mean_density(str,velem,N,weigh,wi,nonempty,kmax,Vmax,numviz)
%***********************************************************************************************************************
%gradvol = grad_volume_mean_density(str,velem,N,weigh,wi,nonempty,kmax,Vmax,numviz)
%
%This function computes the volume constraint gradient vector, when 
%the mean density filter is applied.
%
%INPUT PARAMETERS:
%-----------------
%str: Structure with the problem data.
%velem: Element volumes vector. 
%N: Matrix where each column contains the neighbors of each finite element.
%weigh: Matrix where each row contains the weights of each finite element,
%       associated to the mean density filter or to the sinh filter.
%wi: Vector where each element is the sum of the finite element weights,
%    associated to the mean density filter or to the sinh filter.
%nonempty: Vector with the indexes of the finite elements that are
%          nonempty.
%kmax: Number of nonempty elements.
%Vmax: Maximum volume that the structure can contain.
%numviz: Vector that contains the number of neighbors of each finite
%        element.
%
%OUTPUT PARAMETER:
%-----------------
%gradvol: Volume constraint gradient vector, when the mean density filter is applied.
%***********************************************************************************************************************

gradvol = zeros(str.nelem,1);
for i = 1:kmax
    m1 = nonempty(i);
    viz = numviz(m1);
    ind = N(1:viz,m1);
    gradvol(m1) = sum(weigh(m1,1:viz)'.*velem(ind)./wi(ind))/Vmax;
end
        

