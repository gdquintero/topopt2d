function xvol = mean_density_filter(str,x,weigh,wi,nonempty,kmax,numviz,N)
%************************************************************************************************************
%xvol = mean_density_filter(str,x,weigh,wi,nonempty,kmax,numviz,N)
%
%This function applies the mean density filter to the original element
%densities.
%
%INPUT PARAMETERS:
%-----------------
%str: Structure with the problem data.
%x: Original element densities vector.
%weigh: Matrix where each row contains the weights of each finite element,
%       associated to the mean density filter or to the sinh filter.
%wi: Vector where each element is the sum of the finite element weights,
%    associated to the mean density filter or to the sinh filter.
%nonempty: Vector with the indexes of the finite elements that are nonempty.
%kmax: Number of nonempty elements.
%numviz: Vector that contains the number of neighbors of each finite element.
%N: Matrix where each column contains the neighbors of each finite element.
%
%OUTPUT PARAMETERS:
%------------------
%x: Filtered element densities vector. 
%************************************************************************************************************

xvol = zeros(str.nelem,1); 
for k = 1:kmax
    m1 = nonempty(k);
    viz = numviz(m1);
    ind = N(1:viz,m1);
    xvol(m1) = weigh(m1,1:viz)*x(ind)/wi(m1);
end
       
