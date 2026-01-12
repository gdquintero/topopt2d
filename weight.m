function [weigh,wi] = weight(str,rmin,nonempty,kmax,numviz,distviz)
%************************************************************************************************************
%[weigh,wi] = weight(str,rmin,nonempty,kmax,numviz,distviz)
%
%This function computes the finite elements weights associated to
%the application of the mean density filter.
%
%INPUT PARAMETERS:
%-----------------
%str: Structure with the problem data.
%rmin: Filter radius.
%nonempty: Vector with the indexes of the finite elements that are nonempty.
%kmax: Number of nonempty elements.
%numviz: Vector that contains the number of neighbors of each finite element.
%distviz: Matrix where each column contains the Euclidean distance between
%         the centers of two elements of the neighborhood.
%
%OUTPUT PARAMETERS:
%------------------
%weigh: Matrix where each row contains the weights of each finite element,
%       associated to the mean density filter or to the sinh filter.
%wi: Vector where each element is the sum of the finite element weights,
%    associated to the mean density filter or to the sinh filter.
%************************************************************************************************************

frmin = floor(rmin);
ncol = (2*frmin + 1)^2;
weigh = zeros(str.nelem,ncol);
wi = zeros(str.nelem,1);  
rmin32 = 2*((rmin/3)^2);
twopirmin3 = 2*pi*rmin/3;

for k = 1:kmax
    m1 = nonempty(k);
    viz = numviz(m1);
    for j = 1:viz
        weigh(m1,j) = exp(-((distviz(j,m1))^2)/rmin32)/twopirmin3;
    end
    wi(m1) = sum(weigh(m1,1:viz));
end
        
    
