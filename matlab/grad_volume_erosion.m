function gradvol = grad_volume_erosion(str,gradxer,velem,Vmax,nonempty,numviz,kmax,N)
%***********************************************************************************************************************
%gradvol = grad_volume_erosion(str,gradxer,velem,Vmax,nonempty,numviz,kmax,N)
%
%This function computes the volume constraint gradient vector, when 
%the erosion filter is applied.
%
%INPUT PARAMETERS:
%-----------------
%str: Structure with the problem data.
%gradxer: Matrix where each line contains the coefficients associated to
%         the chain rule for the volume constraint derivatives,
%         when the erosion filter is used.
%velem: Element volumes vector. 
%Vmax: Maximum volume that the structure can contain.
%nonempty: Vector with the indexes of the finite elements that are
%          nonempty.
%numviz: Vector that contains the number of neighbors of each finite
%        element.
%kmax: Number of nonempty elements.
%N: Matrix where each column contains the neighbors of each finite element.
%
%OUTPUT PARAMETER:
%-----------------
%gradvol: Volume constraint gradient vector, when the erosion filter is applied.
%***********************************************************************************************************************

gradvol = zeros(str.nelem,1);
for i = 1:kmax
    m1 = nonempty(i);
    viz = numviz(m1);
    ind = N(1:viz,m1);
    gradvol(m1) = gradxer(m1,ind)*velem(ind)/Vmax;
end